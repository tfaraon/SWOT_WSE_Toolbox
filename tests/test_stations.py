import json
import tempfile
import unittest
from pathlib import Path
import pandas as pd
from examples.make_demo import create
from SWOT_toolbox.run import run
from SWOT_toolbox.stations import station
from SWOT_toolbox.SWOT_tools import extract_wse_timeseries, extract_wse_timeseries_parallel
from SWOT_toolbox.comparison import compare


class StationsTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.temp = tempfile.TemporaryDirectory()
        cls.root = Path(cls.temp.name) / "demo"
        create(cls.root)
        cls.args = (cls.root / "rasters", 137.028098, -28.893022)
        cls.tables = run(cls.root / "config.json")

    @classmethod
    def tearDownClass(cls):
        cls.temp.cleanup()

    def test_gauge_and_lakesp(self):
        rows = self.tables["summary"].query("method == 'pixel' and resolution == '100'").set_index("reference")
        self.assertAlmostEqual(rows.loc["gauge","bias_m"], .2)
        self.assertAlmostEqual(rows.loc["lakesp","bias_m"], -.1)
        self.assertEqual(rows.loc["gauge","matched_passes"], 6)
        self.assertAlmostEqual(rows.loc["coverage","usable_fraction"], 6/7)

    def test_support_and_peak(self):
        audit = self.tables["audit"].query("method == 'pixel'")
        self.assertEqual(set(audit.width_m), {1250, 1300})
        peak = audit[audit.wse_raw.eq(-11)]
        self.assertEqual(len(peak), 2)
        self.assertTrue(peak.iqr_flag.all() and peak.status.eq("accepted").all())
        self.assertEqual(sum(self.tables["audit"].status.eq("iqr")), 2)

    def test_rejections_and_parallel(self):
        kw = dict(width_m=1300, resolution=100, bounds=(-15,-12), min_fraction=.999)
        a = station(*self.args, **kw)
        pd.testing.assert_frame_equal(a, station(*self.args, **kw, workers=2))
        self.assertIn("coverage", a.status.values)
        self.assertIn("no_wse", a.status.values)
        b = station(*self.args, width_m=1300, bounds=(-15,-12))
        self.assertIn("bounds", b.status.values)
        outside = station(self.args[0], 0, 0)
        self.assertTrue(outside.status.eq("outside").all())

    def test_notebook_entry_points(self):
        kw = dict(filter_resolution="100", wse_qual_filter=[0,1])
        a = extract_wse_timeseries(*self.args, **kw)
        b = extract_wse_timeseries_parallel(*self.args, **kw, n_workers=1)
        self.assertEqual(len(a), 5)
        self.assertEqual(a.wse.tolist(), b.wse.tolist())
        new = extract_wse_timeseries_parallel(*self.args, **kw, width_m=1300, quality_mode="pixel")
        self.assertEqual(len(new), 6)

    def test_datum_and_time(self):
        group = self.tables["audit"].query("method == 'pixel' and resolution == '100'")
        spec = dict(kind="gauge", datum="other", tolerance="5min")
        ref = pd.DataFrame(dict(reference_date=pd.to_datetime(["2026-01-01"], utc=True), reference_m=[1.]))
        with self.assertRaises(ValueError):
            compare(group, ref, spec, "demo")
        _, stats = compare(group, ref, dict(spec, datum="demo"), "demo")
        self.assertEqual(stats["matched_passes"], 0)
        empty, stats = compare(group.assign(status="outside"), ref, dict(spec, datum="demo"), "demo")
        self.assertTrue(empty.empty)
        self.assertEqual(stats["eligible_passes"], 0)

    def test_period_boundary(self):
        cfg = json.loads((self.root / "config.json").read_text())
        cfg["sites"][0]["periods"] = [{"name": "peak", "start": "2025-03-05", "end": "2025-03-26"}]
        cfg["methods"] = [{"name": "window", "width_m": 1300}]
        path = self.root / "periods.json"
        path.write_text(json.dumps(cfg))
        res = run(path)
        row = res["summary"].query("phase == 'peak' and resolution == '100' and reference == 'gauge'").iloc[0]
        self.assertEqual(row.matched_passes, 1)
        self.assertAlmostEqual(row.bias_m, .2)

    def test_receding_shore(self):
        row = self.tables["audit"].query(
            "method == 'pixel' and resolution == '100' and cycle == '005'").iloc[0]
        self.assertEqual(row.status, "accepted")
        self.assertFalse(row.centre_valid)
        self.assertEqual(row.optical_centre_water, 0)
        self.assertEqual(row.optical_lag_h, 2)
        self.assertAlmostEqual(row.nearest_valid_m, 200, delta=1)
        self.assertAlmostEqual(row.centroid_distance_m, 400, delta=1)
        self.assertAlmostEqual(row.optical_wet_fraction, 5/13)
        self.assertEqual(row.wse_mad_m, 0)
        pixels = self.tables["pixels"].query(
            "method == 'pixel' and resolution == '100' and cycle == '005'")
        self.assertEqual(len(pixels), 169)
        self.assertEqual(pixels.used.sum(), 65)
        self.assertTrue(pixels.loc[pixels.used, "optical_water"].eq(1).all())

    def test_clipped_window(self):
        row = station(*self.args, radius=10, resolution=100).iloc[0]
        self.assertEqual(row.n_in_tile, 225)
        self.assertAlmostEqual(row.tile_fraction, 225/441)
        self.assertEqual(row.available_fraction_in_tile, 1)
        self.assertAlmostEqual(row.quality_retained_fraction, 224/225)
        self.assertAlmostEqual(row.valid_fraction, 224/441)

    def test_estimator_and_missing_optical_date(self):
        result = station(*self.args, radius=1, levels=None, estimator="median")
        row = result.iloc[0]
        self.assertAlmostEqual(row.wse_mean_m, -13)
        self.assertAlmostEqual(row.wse, -14)
        self.assertAlmostEqual(row.wse_mad_m, 0)
        row = self.tables["audit"].query("method == 'pixel' and cycle == '002'").iloc[0]
        self.assertTrue(pd.isna(row.optical_centre_water))

    def test_common_sample(self):
        rows = self.tables["common_passes"].query("resolution == '100' and reference == 'gauge'")
        self.assertEqual(len(rows), 4)
        self.assertTrue(rows.common_passes.eq(4).all())
        pairs = self.tables["matches"].query("common_pass and resolution == '100' and reference == 'gauge'")
        self.assertTrue(pairs.groupby(["cycle","pass","reference_date"]).size().eq(4).all())

    def test_shared_engine_and_cache_adapter(self):
        from SWOT_toolbox.SWOT_tools import _process_file
        reference = self.tables["audit"].query("method == 'pixel' and resolution == '100' and status == 'accepted'")
        options = dict(width_m=1300, wse_qual_filter=[0,1], quality_mode="pixel", filter_resolution="100", iqr="flag")
        for table in (extract_wse_timeseries(*self.args, **options),
                      extract_wse_timeseries_parallel(*self.args, **options, n_workers=2)):
            self.assertEqual(table.date.tolist(), reference.date.tolist())
            self.assertEqual(table.wse.tolist(), reference.wse.tolist())
        row = _process_file(reference.filename.iloc[0], self.args[1], self.args[2], 0, [0,1], False)
        self.assertEqual(row["wse"], -14)

    def test_holdout_iqr(self):
        # An insufficient calibration period must never borrow evaluation dates.
        frame = station(*self.args, radius=1, iqr_fit_period=["2025-01-01","2025-02-01"])
        self.assertFalse(frame.iqr_flag.any())

    def test_real_filename_duplicate_selection(self):
        import shutil
        from netCDF4 import Dataset
        from SWOT_toolbox.stations import select_acquisitions
        src = sorted(self.args[0].glob("*100m*"))[0]
        stem = "SWOT_L2_HR_Raster_100m_UTM53J_N_x_x_x_007_116_102F_20231127T131422_20231127T131443"
        with tempfile.TemporaryDirectory() as tmp:
            directory = Path(tmp)
            for ending in ["PIC0_05","PID0_01","PID0_10","PID0_2","PID0_10 (1)"]:
                shutil.copy(src, directory / f"{stem}_{ending}.nc")
            chosen = directory / f"{stem}_PID0_10.nc"
            with Dataset(chosen, "r+") as ds:
                ds["wse"][:] = -12
            sub = directory / "z_cache"
            sub.mkdir()
            shutil.copy(chosen, sub / chosen.name)
            result = extract_wse_timeseries_parallel(directory, self.args[1], self.args[2],
                buffer_size=6, filter_outliers=True, filter_bound=True,
                wse_qual_filter=[0,1,2,3], mnt_corr=None)
            self.assertEqual(len(result), 1)
            row = result.iloc[0]
            self.assertEqual(row.filename, str(chosen))
            self.assertEqual(row.duplicate_count, 5)
            self.assertEqual(len(json.loads(row.skipped_files)), 5)
            self.assertEqual(row.wse, -12)
            with self.assertRaisesRegex(ValueError, "Duplicate acquisition"):
                station(directory, *self.args[1:], duplicate_policy="error")
            other = [Path(str(chosen).replace(old,new)) for old,new in
                     [("53J","54J"),("131443","131444"),("102F","103F"),("100m","250m")]]
            selected, _ = select_acquisitions([(p,) for p in [chosen, *other]])
            self.assertEqual(len(selected), 5)

    def test_unreadable_granules(self):
        import shutil
        import warnings
        src = sorted(self.args[0].glob("*100m*"))[0]
        with tempfile.TemporaryDirectory() as tmp, warnings.catch_warnings():
            warnings.simplefilter("ignore", RuntimeWarning)
            root = Path(tmp)
            stem = src.name.removesuffix("demo.nc")
            good = root / (stem + "PIC0_01.nc")
            bad = root / (stem + "PID0_01.nc")
            bad.write_bytes(b"truncated")
            shutil.copy(src, good)
            other = root / bad.name.replace("_001_394_", "_002_394_")
            other.write_bytes(b"not a NetCDF file")
            for workers in (1, 2):
                result = station(root, *self.args[1:], workers=workers)
                recovered = result[result.status.eq("accepted")].iloc[0]
                self.assertEqual(recovered.filename, str(good))
                self.assertEqual(json.loads(recovered.read_errors)[0]["filename"], str(bad))
                self.assertEqual(result.status.eq("unreadable").sum(), 1)
                out = extract_wse_timeseries_parallel(root, *self.args[1:], buffer_size=6,
                    filter_outliers=True, filter_bound=True, wse_qual_filter=[0,1,2,3],
                    mnt_corr=None, n_workers=workers)
                self.assertEqual(len(out), 1)
                self.assertEqual({err["filename"] for err in out.attrs["read_errors"]}, {str(bad), str(other)})
            good.unlink()
            out = extract_wse_timeseries_parallel(root, *self.args[1:], filter_bound=True, filter_outliers=True)
            self.assertTrue(out.empty)
            self.assertEqual(len(out.attrs["read_errors"]), 2)
