import hashlib
import sys
import unittest
from pathlib import Path
from tempfile import TemporaryDirectory
from unittest.mock import Mock, patch

import pandas as pd

from examples.make_demo import create
from SWOT_toolbox.catalog import COLUMNS, download_selected, extract_selection


class CatalogueTest(unittest.TestCase):
    def test_download_retry_integrity_and_selection(self):
        with TemporaryDirectory() as temp:
            root = Path(temp)
            create(root / "demo")
            source = next((root / "demo/rasters").rglob("*100m*.nc"))
            data = source.read_bytes()
            dest = root / source.name
            row = dict(filename=source.name, local_path=str(dest), size_bytes=len(data),
                       checksum=hashlib.md5(data).hexdigest(), algorithm="MD5", url="https://example.test/file")
            selection = pd.DataFrame([row], columns=COLUMNS)
            api = Mock()
            with patch.dict(sys.modules, earthaccess=api):
                failed = download_selected(selection)
                self.assertEqual(failed.local_state.iloc[0], "download_error")
                self.assertFalse(dest.exists())
                def transfer(urls, folder, provider):
                    (Path(folder) / source.name).write_bytes(data)
                api.download.side_effect = transfer
                result = download_selected(selection)
                download_selected(result)
                self.assertEqual(api.download.call_count, 2)
            audit = extract_selection(result, 137.028098, -28.893022, root)
            self.assertEqual(len(audit), 1)
            self.assertEqual(audit.local_state.iloc[0], "available")
            dest.write_bytes(bytes(len(data)))
            audit = extract_selection(result, 137.028098, -28.893022, root)
            self.assertEqual(audit.local_state.iloc[0], "invalid_checksum")
            self.assertEqual(audit.status.iloc[0], "not_processed")

    def test_empty_search(self):
        empty = pd.DataFrame(columns=COLUMNS)
        self.assertTrue(extract_selection(empty, 0, 0, ".").empty)
