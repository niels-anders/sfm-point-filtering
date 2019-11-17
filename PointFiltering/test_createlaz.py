from unittest import TestCase
import LazData as LazData
import Grid as Grid


FN_LAS = "python/PointFiltering/tests/testdata/las.las"


class TestCreatelaz(TestCase):

    def test_createlaz(self):
        las = LazData.createlaz(FN_LAS)
        self.assertIsNone(las.get_header().records_count)

