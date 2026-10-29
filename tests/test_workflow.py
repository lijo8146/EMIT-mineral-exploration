import json
import sys
import tempfile
import unittest
from pathlib import Path

import numpy as np
from netCDF4 import Dataset

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
import workflow as w


class WorkflowTests(unittest.TestCase):
    def setUp(self):
        self.cfg = json.loads((w.HERE / "config.json").read_text(encoding="utf-8"))
        self.scratch = Path.cwd() / "work"
        self.scratch.mkdir(exist_ok=True)

    def test_excluded_and_invalid_uncertainty_never_detected(self):
        g = {"id": np.ones((2, 3)), "depth": np.full((2, 3), .05),
             "uncertainty": np.array([[.01, 0, np.nan], [.01, -.01, .1]]),
             "fit": np.ones((2, 3))}
        clear = np.ones((2, 3), bool)
        clear[1, 0] = False
        target, _, _, _ = w.classify([g], {1: "Kaolinite"}, clear, self.cfg)
        np.testing.assert_array_equal(target, [[True, False, False], [False, False, False]])

    def test_unverified_fit_scaling_rejected(self):
        self.cfg["min_fit"] = .7
        a = np.ones((2, 2))
        with self.assertRaises(ValueError):
            w.classify([dict(id=a, depth=a, uncertainty=a, fit=a)],
                       {1: "Kaolinite"}, a.astype(bool), self.cfg)

    def test_native_netcdf_adapter_and_glt(self):
        with tempfile.TemporaryDirectory(dir=self.scratch) as tmp:
            tmp = Path(tmp)
            paths = [tmp / f"EMIT_{prefix}_001_20230101T000000_2300101_001.nc"
                     for prefix in ("L2B_MIN", "L2B_MINUNCERT", "L2A_MASK")]
            for idx, path in enumerate(paths):
                with Dataset(path, "w") as ds:
                    ds.createDimension("downtrack", 2)
                    ds.createDimension("crosstrack", 2)
                    if idx < 2:
                        for g in (1, 2):
                            variables = [("mineral_id", [[1, 2], [2, 1]]), ("band_depth", np.full((2, 2), .05))] if idx == 0 else [
                                ("band_depth_unc", np.full((2, 2), .01)), ("fit", np.full((2, 2), .4))]
                            for name, values in variables:
                                ds.createVariable(f"group_{g}_{name}", "f4", ("downtrack", "crosstrack"))[:] = values
                    if idx == 0:
                        ds.geotransform = [-117.3, .001, 0, 37.5, 0, -.001]
                        loc = ds.createGroup("location")
                        loc.createDimension("ortho_y", 2)
                        loc.createDimension("ortho_x", 2)
                        loc.createVariable("glt_x", "i4", ("ortho_y", "ortho_x"))[:] = [[2, 1], [0, 1]]
                        loc.createVariable("glt_y", "i4", ("ortho_y", "ortho_x"))[:] = [[1, 1], [0, 2]]
                        meta = ds.createGroup("mineral_metadata")
                        meta.createDimension("minerals", 2)
                        meta.createVariable("index", "i4", ("minerals",))[:] = [1, 2]
                        meta.createVariable("name", str, ("minerals",))[:] = np.array(["Kaolinite", "Alunite"], dtype=object)
                    if idx == 2:
                        ds.createDimension("bands", 5)
                        bands = ds.createGroup("sensor_band_parameters")
                        bands.createVariable("mask_bands", str, ("bands",))[:] = np.array([
                            "Cloud flag", "Cirrus flag", "Water flag", "Spacecraft Flag", "AOD550"], dtype=object)
                        mask = np.zeros((2, 2, 5))
                        mask[..., 4] = .2  # Atmospheric quantities must not exclude pixels.
                        mask[0, 1, 3] = 1  # Spacecraft must exclude this pixel.
                        ds.createVariable("mask", "f4", ("downtrack", "crosstrack", "bands"))[:] = mask
            groups, lookup, clear, transform = w.read_scene(*paths)
            self.assertEqual(lookup[2], "Alunite")
            self.assertEqual(groups[0]["id"][0, 0], 2)
            self.assertTrue(np.isnan(groups[0]["id"][1, 0]))
            np.testing.assert_array_equal(clear, [[False, True], [False, True]])
            self.assertAlmostEqual(transform.c, -117.3)
            with self.assertRaises(ValueError):
                w.read_scene(paths[0], paths[1], tmp / "wrong_scene.nc")

    def test_exports_and_geographic_polygons(self):
        import rasterio
        with tempfile.TemporaryDirectory(dir=self.scratch) as tmp:
            out = Path(tmp)
            from rasterio import Affine
            # Minimal fixture exercises eligible, excluded and detected export values.
            a = np.ones((5, 5))
            clear = a.astype(bool)
            clear[0, 0] = False
            ids = np.zeros((5, 5))
            ids[1:4, 1:4] = 1
            groups = [dict(id=ids, depth=a * .05, uncertainty=a * .01, fit=a)]
            w.analyze(groups, {1: "Kaolinite"}, clear,
                      Affine(.001, 0, -117.3, 0, -.001, 37.5),
                      self.cfg, out, {"inputs": {"fixture": "test-only"}})
            run = json.loads((out / "run.json").read_text())
            self.assertNotIn("synthetic", run)
            self.assertEqual(run["target_count"], 1)
            features = json.loads((out / "targets.geojson").read_text())["features"]
            for f in features:
                self.assertNotIn("synthetic", f["properties"])
                self.assertGreater(f["properties"]["area_ha"], 0)
                self.assertGreaterEqual(f["properties"]["pixels"], self.cfg["min_pixels"])
                self.assertGreater(f["properties"]["longitude"], -118)
            with rasterio.open(out / "target_detection.tif") as ds:
                self.assertNotIn("synthetic", ds.tags())
                self.assertEqual(ds.crs.to_epsg(), 4326)
                self.assertEqual(ds.nodata, -9999)
                self.assertEqual(set(np.unique(ds.read(1))), {-9999., 0., 1.})


if __name__ == "__main__":
    unittest.main()
