"""EMIT V001 scene screening. Run --help; no credentials are embedded."""
import argparse
import csv
import hashlib
import json
from pathlib import Path
from urllib.parse import urlencode
from urllib.request import urlopen

import numpy as np
from scipy import ndimage

HERE = Path(__file__).resolve().parent


def discover(cfg, out):
    """Public CMR inventory; downloading is a separate, explicit step."""
    out.mkdir(parents=True, exist_ok=True)
    inventory = {}
    for product in ("EMITL2BMIN", "EMITL2ARFL"):
        entries, page = [], 1
        while True:
            query = urlencode({"short_name": product, "version": cfg["version"],
                "bounding_box": ",".join(map(str, cfg["bbox"])),
                "temporal": ",".join(cfg["temporal"]), "page_size": 200,
                "page_num": page, "sort_key": "-start_date"})
            with urlopen("https://cmr.earthdata.nasa.gov/search/granules.json?" + query,
                         timeout=90) as response:
                batch = json.load(response)["feed"]["entry"]
            entries.extend(batch)
            if len(batch) < 200:
                break
            page += 1
        inventory[product] = entries
        print(f"{product}: {len(entries)} intersecting granules")
    (out / "inventory.json").write_text(json.dumps(inventory, indent=2), encoding="utf-8")


def download(cfg, out, granule_id, product):
    import earthaccess
    earthaccess.login()
    results = earthaccess.search_data(short_name=product, version=cfg["version"],
                                     granule_name=granule_id)
    if len(results) != 1:
        raise ValueError("Specify exactly one producer_granule_id from inventory.json")
    earthaccess.download(results, str(out))


def array(variable):
    return np.ma.filled(variable[:].astype(float), np.nan)


def read_scene(mineral, uncertainty, quality):
    """Read native V001 geometry; reject unpaired acquisitions or missing quality."""
    from netCDF4 import Dataset, chartostring
    from rasterio import Affine
    key = Path(mineral).stem.split("_001_")[-1]
    for path in (uncertainty, quality):
        if Path(path).stem.split("_001_")[-1] != key:
            raise ValueError("MIN, MINUNCERT and MASK must share version/acquisition/orbit/scene")
    with Dataset(mineral) as ds, Dataset(uncertainty) as unc, Dataset(quality) as mask:
        loc = ds.groups["location"]
        gx, gy = array(loc["glt_x"]), array(loc["glt_y"])
        if gx.shape != gy.shape:
            raise ValueError("GLT shapes differ")
        transform = Affine.from_gdal(*ds.getncattr("geotransform"))
        meta = ds.groups["mineral_metadata"]
        names = meta["name"][:]
        if names.dtype.kind == "S" and names.ndim > 1:
            names = chartostring(names)
        names = [n.decode() if isinstance(n, bytes) else str(n) for n in names]
        lookup = dict(zip(np.asarray(meta["index"][:], dtype=int), names))
        raw_shape = ds["group_1_mineral_id"].shape
        valid_glt = (np.isfinite(gx) & np.isfinite(gy) & (gx >= 1) & (gy >= 1)
                     & (gx <= raw_shape[1]) & (gy <= raw_shape[0]))
        xx = np.where(valid_glt, gx - 1, 0).astype(int)
        yy = np.where(valid_glt, gy - 1, 0).astype(int)

        def ortho(a):
            if a.shape[:2] != raw_shape:
                raise ValueError("Raw arrays have different shapes")
            result = a[yy, xx].copy()
            result[~valid_glt] = np.nan
            return result

        # Select binary masks by metadata; never interpret atmospheric quantities as flags.
        raw_mask = array(mask["mask"])
        if raw_mask.ndim != 3 or raw_mask.shape[-1] < 3:
            raise ValueError("Expected V001 mask (downtrack, crosstrack, bands)")
        band_names = mask.groups["sensor_band_parameters"]["mask_bands"][:]
        if band_names.dtype.kind == "S" and band_names.ndim > 1:
            band_names = chartostring(band_names)
        band_names = [n.decode() if isinstance(n, bytes) else str(n) for n in band_names]
        selected = [i for i, name in enumerate(band_names)
                    if any(t in name.lower() for t in ("cloud", "cirrus", "water flag", "spacecraft", "aggregate"))]
        if not all(any(term in n.lower() for n in band_names) for term in ("cloud", "cirrus", "water flag")):
            raise ValueError("Missing required cloud/cirrus/water mask labels")
        flags = raw_mask[..., selected]
        clear = np.all(np.isfinite(flags) & (flags == 0), axis=-1)
        groups = []
        for g in (1, 2):
            prefix = f"group_{g}_"
            groups.append({"id": ortho(array(ds[prefix + "mineral_id"])),
                "depth": ortho(array(ds[prefix + "band_depth"])),
                "uncertainty": ortho(array(unc[prefix + "band_depth_unc"])),
                "fit": ortho(array(unc[prefix + "fit"]))})
        return groups, lookup, ortho(clear.astype(float)) == 1, transform


def classify(groups, lookup, clear, cfg):
    target = np.zeros(clear.shape, bool)
    support = target.copy()
    best_snr = np.zeros(clear.shape)
    matched = set()
    for g in groups:
        depth, uncertainty = g["depth"], g["uncertainty"]
        snr = np.divide(depth, uncertainty, out=np.full(depth.shape, np.nan),
                        where=np.isfinite(uncertainty) & (uncertainty > 0))
        good = clear & np.isfinite(depth) & (depth >= cfg["min_depth"]) & (snr >= cfg["min_snr"])
        if cfg["min_fit"] is not None:
            if cfg["fit_multiplier"] is None:
                raise ValueError("Set a verified product-specific fit_multiplier before min_fit")
            good &= np.isfinite(g["fit"]) & (g["fit"] * cfg["fit_multiplier"] >= cfg["min_fit"])
        for mid, name in lookup.items():
            is_target = any(term in name.lower() for term in cfg["target_terms"])
            is_support = any(term in name.lower() for term in cfg["support_terms"])
            pixels = good & (g["id"] == mid)
            if is_target:
                target |= pixels
                best_snr[pixels] = np.maximum(best_snr[pixels], snr[pixels])
                if pixels.any():
                    matched.add(name)
            if is_support:
                support |= pixels
    return target, support, best_snr, sorted(matched)


def sha256(path):
    h = hashlib.sha256()
    with open(path, "rb") as f:
        for chunk in iter(lambda: f.read(1024 * 1024), b""):
            h.update(chunk)
    return h.hexdigest()


def analyze(groups, lookup, clear, transform, cfg, out, provenance):
    import rasterio
    from rasterio.features import shapes
    from rasterio.transform import xy
    from pyproj import Geod
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    out.mkdir(parents=True, exist_ok=True)
    # Crop eligibility to AOI while preserving the native geographic grid.
    rr, cc = np.indices(clear.shape)
    longitude = transform.c + (cc + .5) * transform.a + (rr + .5) * transform.b
    latitude = transform.f + (cc + .5) * transform.d + (rr + .5) * transform.e
    west, south, east, north = cfg["bbox"]
    clear &= (longitude >= west) & (longitude <= east) & (latitude >= south) & (latitude <= north)
    target, support, snr, matched = classify(groups, lookup, clear, cfg)
    labels, n = ndimage.label(target)  # 4-connectivity agrees with polygonization
    sizes = np.bincount(labels.ravel())
    keep = sizes >= cfg["min_pixels"]
    keep[0] = False
    labels = np.where(keep[labels], labels, 0).astype("int32")
    records = []
    geod = Geod(ellps="WGS84")
    features = []
    for geometry, value in shapes(labels, mask=labels > 0, transform=transform):
        label = int(value)
        pixels = labels == label
        neighborhood = ndimage.binary_dilation(pixels, iterations=3) & clear
        support_fraction = float(support[neighborhood].mean()) if neighborhood.any() else 0.
        median_snr = float(np.median(snr[pixels]))
        # Transparent screening heuristic, NOT a calibrated deposit probability.
        score = 100 * (.6 * min(median_snr / 10, 1) + .4 * support_fraction)
        rings = geometry["coordinates"]
        area = abs(geod.polygon_area_perimeter(*zip(*rings[0]))[0])
        area -= sum(abs(geod.polygon_area_perimeter(*zip(*ring))[0]) for ring in rings[1:])
        row = {"target_id": label, "pixels": int(pixels.sum()), "area_ha": area / 10000,
            "longitude": float(longitude[pixels].mean()), "latitude": float(latitude[pixels].mean()),
            "median_snr": median_snr, "nearby_iron_fraction": support_fraction,
            "screening_score": score}
        records.append(row)
        features.append({"type": "Feature", "geometry": geometry, "properties": row})
    records.sort(key=lambda r: r["screening_score"], reverse=True)
    fields = ["target_id", "pixels", "area_ha", "longitude", "latitude", "median_snr",
              "nearby_iron_fraction", "screening_score"]
    with open(out / "targets.csv", "w", newline="", encoding="utf-8") as f:
        writer = csv.DictWriter(f, fieldnames=fields)
        writer.writeheader()
        writer.writerows(records)
    (out / "targets.geojson").write_text(json.dumps({"type": "FeatureCollection", "features": features}, indent=2))
    layers = {"valid_surface": clear.astype(float), "target_detection": np.where(clear, target, np.nan),
              "target_snr": np.where(target, snr, np.nan),
              "target_regions": np.where(clear, labels, np.nan)}
    for name, data in layers.items():
        with rasterio.open(out / f"{name}.tif", "w", driver="GTiff", height=data.shape[0],
            width=data.shape[1], count=1, dtype="float32", crs="EPSG:4326", transform=transform,
            nodata=-9999, compress="deflate") as dst:
            dst.write(np.where(np.isfinite(data), data, -9999).astype("float32"), 1)
            dst.update_tags(description=name)
    fig, axes = plt.subplots(1, 3, figsize=(12, 4), constrained_layout=True)
    for ax, data, title in zip(axes, [clear, np.where(clear, target, np.nan), np.where(target, snr, np.nan)],
                               ["Eligible surface", "Alteration candidates", "Detection depth / uncertainty"]):
        im = ax.imshow(data, interpolation="nearest")
        ax.set_title(title)
        ax.set_xlabel("Native grid column")
        ax.set_ylabel("Native grid row")
        fig.colorbar(im, ax=ax, shrink=.7)
    fig.suptitle(cfg["name"])
    fig.savefig(out / "overview.png", dpi=150)
    plt.close(fig)
    provenance.update({"config": cfg, "matched_library_entries": matched,
        "eligible_pixels": int(clear.sum()), "candidate_pixels": int(target.sum()),
        "target_count": len(records), "ranking": "100*(0.6*min(median_snr/10,1)+0.4*nearby_iron_fraction)"})
    (out / "run.json").write_text(json.dumps(provenance, indent=2), encoding="utf-8")
    print(f"Saved {len(records)} regions to {out.resolve()}")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("command", choices=["discover", "download", "analyze"])
    parser.add_argument("--config", type=Path, default=HERE/"config.json")
    parser.add_argument("--out", type=Path, default=HERE/"results")
    parser.add_argument("--granule-id")
    parser.add_argument("--product", choices=["EMITL2BMIN", "EMITL2ARFL"], default="EMITL2BMIN")
    parser.add_argument("--mineral", type=Path)
    parser.add_argument("--uncertainty", type=Path)
    parser.add_argument("--quality", type=Path)
    args = parser.parse_args()
    cfg = json.loads(args.config.read_text(encoding="utf-8"))
    if cfg["version"] != "001":
        parser.error("This adapter supports V001 only; verify schema and masks before other versions")
    if args.command == "discover":
        discover(cfg, args.out)
    elif args.command == "download":
        if not args.granule_id:
            parser.error("download requires --granule-id")
        download(cfg, args.out, args.granule_id, args.product)
    else:
        if not all((args.mineral, args.uncertainty, args.quality)):
            parser.error("analyze requires --mineral, --uncertainty and --quality")
        groups, lookup, clear, transform = read_scene(args.mineral, args.uncertainty, args.quality)
        provenance = {"inputs": {str(p.resolve()): sha256(p)
                      for p in (args.mineral, args.uncertainty, args.quality)}}
        analyze(groups, lookup, clear, transform, cfg, args.out, provenance)


if __name__ == "__main__":
    main()