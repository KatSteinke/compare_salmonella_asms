import re
import unittest


import pandas as pd
import pytest

import snake_helpers

class TestGetReadsFromOverview(unittest.TestCase):
    read_overview = pd.DataFrame(data={"sample": ["1199123456-1", "1199123456-2"],
                                       "r1": ["/path/to/r1", pd.NA],
                                       "r2": ["/path/to/r2", pd.NA],
                                       "extra": [pd.NA, "/path/to/long_reads"]})

    def test_success(self):
        """Successfully get the reads from the overview."""
        expected_r1 = "/path/to/r1"
        test_r1 = snake_helpers.get_reads_from_overview("1199123456-1", self.read_overview, "r1")
        assert expected_r1 == test_r1

        expected_r2 = "/path/to/r2"
        test_r2 = snake_helpers.get_reads_from_overview("1199123456-1", self.read_overview, "r2")
        assert expected_r2 == test_r2

        expected_long = "/path/to/long_reads"
        test_long = snake_helpers.get_reads_from_overview("1199123456-2", self.read_overview,
                                                          "extra")
        assert expected_long == test_long

    def test_invalid_type(self):
        """Handle a bad read type."""
        error_msg = "long is not a valid read type. Valid read types are ['r1', 'r2', 'extra']."
        with pytest.raises(ValueError, match=re.escape(error_msg)):
            snake_helpers.get_reads_from_overview("1199123456-2", self.read_overview,
                                                  "long")

    def test_no_reads(self):
        """Fail on missing reads."""
        error_msg = "No extra reads found for 1199123456-1."
        with pytest.raises(KeyError, match = re.escape(error_msg)):
            snake_helpers.get_reads_from_overview("1199123456-1", self.read_overview,
                                                  "extra")

    def test_no_such_number(self):
        """Fail on missing samples."""
        error_msg = "Isolate 1199123456-3 not found in read overview."
        with pytest.raises(KeyError, match = re.escape(error_msg)):
            snake_helpers.get_reads_from_overview("1199123456-3", self.read_overview,
                                                  "extra")
