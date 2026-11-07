import tempfile
import unittest
from pathlib import Path

from blast_search import parse_blast


class AnnotationOutputTests(unittest.TestCase):
    def test_empty_blast_xml_is_rejected(self):
        with tempfile.TemporaryDirectory() as directory:
            xml_path = Path(directory) / "empty.xml"
            csv_path = Path(directory) / "annotation.csv"
            xml_path.write_text("", encoding="utf-8")
            with self.assertRaises(ValueError):
                parse_blast(xml_path, csv_path, 100)


if __name__ == "__main__":
    unittest.main()
