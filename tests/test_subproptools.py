import pkg_resources
import os # file system stuff
import pandas as pd #data frames
import math #sqrt
import numpy as np #arrays
from subproptools import qtaimExtract as qt

def load_sumfile(filename):
    stream = pkg_resources.resource_stream(__name__,'test_data/'+filename+'.sum')
    return stream.readlines()