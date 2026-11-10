import os
from unittest import TestCase
import pytest
import pyrtlib

THIS_DIR = os.path.dirname(os.path.abspath(__file__))


class Test(TestCase):
    def test_version_attribute_is_available(self):
        """pyrtlib.__version__ should be accessible lazily."""
        version = pyrtlib.__version__
        assert isinstance(version, str)
        assert version  # non vuota


    def test_unknown_attribute_raises_attribute_error(self):
        """Accessing an unknown attribute should raise AttributeError."""
        with pytest.raises(AttributeError) as excinfo:
            pyrtlib.__dummy__
        assert "__dummy__" in str(excinfo.value)
