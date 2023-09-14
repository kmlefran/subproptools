"""Tests to ensure the subproptools module works as intended"""
# import math  # sqrt
# import os  # file system stuff

# import numpy as np  # arrays
# import pandas as pd  # data frames
import pkg_resources

# from subproptools import qtaim_extract as qt


def load_sumfile(filename):
    """Test loading a file"""
    stream = pkg_resources.resource_stream(__name__, "test_data/" + filename + ".sum")
    return stream.readlines()
