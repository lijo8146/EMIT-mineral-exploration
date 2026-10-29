# EMIT mineral exploration pilot

Python workflow for screening advanced argillic alteration around **Cuprite, Nevada**.
This is an alteration mapping benchmark.
The initial AOI is longitude -117.35 to -117.10, latitude 37.40 to 37.65.

## Quick start

Use Python 3.11 or newer. From this folder:

```powershell
python -m venv .venv
.venv\Scripts\Activate.ps1
python -m pip install -r requirements.txt
python workflow.py discover --out inventory
python -m unittest discover -s tests -v
```

## Workflow

1. Run `python workflow.py discover --out inventory`. This queries public NASA CMR for
   V001 L2B mineral and L2A reflectance collections intersecting the AOI. Inspect
   `inventory/inventory.json`: `producer_granule_id`, `time_start`, footprint and browse links.
   Intersection does not guarantee full AOI coverage or clear ground. Choose a scene with
   exposed rock and a clear browse image. If empty, expand the date window/AOI or inspect
   current collections; this adapter deliberately supports V001 only.
2. Download exactly one selected mineral granule:
   `python workflow.py download --granule-id EXACT_MINERAL_PRODUCER_ID --out data`.
   The command uses `earthaccess.login()` and downloads the granule assets. NASA Earthdata
   authentication is required. Never put credentials into source or config files.
3. Find its matching acquisition/orbit/scene in the reflectance inventory, then run
   `python workflow.py download --product EMITL2ARFL --granule-id EXACT_REFLECTANCE_PRODUCER_ID --out data`.
   Obtain the accompanying MASK asset from that granule along with reflectance and uncertainty.
   Keep MIN, MINUNCERT, MASK and RFL from the same acquisition and version. Downloads can
   be large; review the inventory before choosing a granule.
4. Run:

```powershell
python workflow.py analyze --mineral data/EMIT_L2B_MIN_001_TIMESTAMP_ORBIT_SCENE.nc --uncertainty data/EMIT_L2B_MINUNCERT_001_TIMESTAMP_ORBIT_SCENE.nc --quality data/EMIT_L2A_MASK_001_TIMESTAMP_ORBIT_SCENE.nc --out results
```

Replace the filenames with actual paths. Each invocation analyzes one scene independently.
Do not mix overlapping scenes before checking temporal consistency.

## Analysis implemented

- Reads the native NetCDF mineral library `index` and `name` fields, both spectral groups,
  depths, uncertainty and fit scores. It does not assume fixed IDs for individual minerals.
- Applies one-based GLT indices to the native geographic grid; invalid GLT locations fail closed.
- Uses named cloud, cirrus, water, spacecraft, dilated-cloud and aggregate flags when present;
  flags must be finite and zero. Limits eligible pixels to the configured bounding box.
- Screens candidate library entries using explicit name substrings for kaolinite, alunite,
  dickite and pyrophyllite. Mixture names can match; this is screening, not pure-mineral
  identification or quantitative unmixing. Review `matched_library_entries` in `run.json`.
- Requires depth >= 0.01 and depth/standard uncertainty >= 3 by default. These are tunable
  pilot thresholds, not validated mineral-specific detection limits. Invalid/zero uncertainty
  is rejected. Hematite/goethite/jarosite provide secondary spatial context.
- Groups candidates with 4-connectivity; keeps regions of at least nine native grid cells.
  Geographic GLT grid cells are not necessarily independent 60 m observations, so the cell
  count is not a count of independent detections. Polygon areas use WGS84 geodesic calculations.
- Ranks each region with `100*(0.6*min(median_snr/10,1)+0.4*nearby_iron_fraction)`.
  The iron fraction uses the region plus a three-cell dilation. This transparent heuristic
  is neither a calibrated confidence nor a probability of a deposit. Fit filtering is off
  until you explicitly set `min_fit` and a verified `fit_multiplier`: NASA's older tutorial
  describes a factor of two, which should not be applied blindly to every processing revision.
- Exports GeoTIFFs, WGS84 GeoJSON polygons, ranked CSV, overview image and config/input SHA256
  provenance. No detection and excluded/nodata pixels are distinct in the detection raster.

## Outputs

| File | Meaning |
|---|---|
| `valid_surface.tif` | 1 eligible, 0 excluded or outside AOI; no distinction among exclusion causes |
| `target_detection.tif` | 1 candidate, 0 eligible without candidate, -9999 excluded |
| `target_snr.tif` | Depth/uncertainty for candidate pixels; -9999 elsewhere |
| `target_regions.tif` | Region ID, 0 eligible without retained region, -9999 excluded |
| `targets.csv`, `targets.geojson` | Region metrics and heuristic rankings |
| `overview.png` | Quick QA image in native map-grid row/column coordinates |
| `run.json` | Configuration, matched entries, input hashes |

## Geological validation before target selection

The implemented code ends at spectral screening. Complete these interpretation steps before
calling a region an exploration target:

1. Examine matched library entries, including mixtures and experimental mineral classes.
2. Inspect associated L2A observed spectra and uncertainty against resampled USGS reference
   spectra; check bad/interpolated bands, low illumination, vegetation and atmospheric residuals.
   Full reflectance spectral fitting, vegetation/shadow filtering and spectral plots are **not
   yet automated** in this version. Mask eligibility alone does not ensure high-quality exposed rock.
3. Check persistence in a second independent observation, where available. Consider sensitivity
   at SNR thresholds 2, 3 and 5 and minimum depths 0.005, 0.01 and 0.02; retain robust regions.
4. Overlay published lithology, alteration, faults and known occurrences. Distinguish transported
   cover and weathering from bedrock alteration; assess steam-heated versus other origins.
5. Compare against independent AVIRIS/ASTER maps and field spectra; sample appropriate exposures
   for XRD and geochemistry. Record access constraints and validation status in a target register.

Mineral spectral evidence alone cannot establish metal content,
subsurface geometry or economic viability. The 16 June 2026 Cuprite scene has been
processed; see `results/REPORT.md` for the observed
spectral anomalies and limitations. No economic deposit findings are claimed.

## Sources and technical references

- [NASA L2B V001 product schema](https://www.earthdata.nasa.gov/data/catalog/lpcloud-emitl2bmin-001)
- [NASA mineralogy tutorial: GLT, library lookup, uncertainty and fit caveats](https://nasa.github.io/LPDAAC-Data-Resources/external/Working_with_EMIT_L2B_Mineralogy.html)
- [NASA reference Python tools](https://github.com/nasa/EMIT-Data-Resources)
- [USGS Cuprite advanced argillic alteration study](https://www.usgs.gov/publications/mapping-advanced-argillic-alteration-cuprite-nevada-using-imaging-spectroscopy)

## Verification scope

Automated tests use small, isolated fixtures to verify screening, invalid uncertainty,
quality exclusion, scene pairing, GLT indexing and GIS exports. A small fabricated V001 NetCDF
fixture tests the file adapter. Live CMR inventory and protected downloads subsequently succeeded,
and full-size V001 MIN, MINUNCERT and MASK files were processed. Nine threshold combinations
were checked on the AOI. Reflectance fitting and field validation remain outstanding.
Requirements are minimum versions, not a locked environment.

The CLI requires scene files for analysis and has no demo mode yet.
The selected acquisition metadata is stored in `selection.json`.
