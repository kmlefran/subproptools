"""Tests to ensure the subproptools module works as intended"""
# import math  # sqrt
# import os  # file system stuff

import os

import numpy as np  # arrays
import pandas as pd  # data frames
import pkg_resources

from subproptools import qtaim_extract as qt

# from subproptools import qtaim_extract as qt


def load_sumfile(filename):
    """Test loading a file"""
    stream = pkg_resources.resource_stream(__name__, "test_data/" + filename + ".sum")
    return stream.readlines()


def test_atomic_properties():
    """Test building atomic properties"""
    data = load_sumfile("SubH_NHOCH3_wb97xd_aug-cc-pvtz_reor")
    atomic_properties = qt.get_atomic_props(data)
    assert isinstance(atomic_properties, dict)
    assert len(atomic_properties.keys()) == 8
    assert (
        "N1" in atomic_properties.keys()  # pylint:disable=consider-iterating-dictionary
    )  # pylint:disable=consider-iterating-dictionary
    assert (
        "q"
        in atomic_properties[
            "N1"
        ].keys()  # pylint:disable=consider-iterating-dictionary
    )


def test_bcp_properties():
    """Test extraction of bcp properties"""
    data = load_sumfile("SubH_NHOCH3_wb97xd_aug-cc-pvtz_reor")
    bcp_properties = qt.get_bcp_properties(data, ["N1", "H2"])
    assert isinstance(bcp_properties, dict)


def test_ldm():
    """Test constuction of the LDM matrix"""
    data = load_sumfile("SubH_NHOCH3_wb97xd_aug-cc-pvtz_reor")
    ldm = qt.get_ldm(data)
    a_props = {  # dummy atomic property data
        "N1": {"q": -0.41769231890},  # actual charge
        "H2": {},
        "H3": {},
        "O4": {},
        "C5": {},
        "H6": {},
        "H7": {},
        "H8": {},
    }
    n_atoms = len(a_props.keys())
    assert isinstance(ldm, pd.DataFrame)
    assert len(ldm) == n_atoms
    assert len(ldm.columns) == n_atoms
    assert np.abs(7 - np.sum(ldm["N1"]) - a_props["N1"]["q"]) < 0.01


def test_get_sub_di():
    """Test that getting substituent di returns a float"""
    data = load_sumfile("SubH_NHOCH3_wb97xd_aug-cc-pvtz_reor")
    sub_di = qt.get_sub_di(data, ["N1", "H3", "O4", "C5"])
    assert isinstance(sub_di, float)
    assert sub_di > 0


def test_get_cc_props():
    """Test that getting charge concentration properties returns a dict"""
    file_with_relpath = (
        os.getcwd() + "/" + "test_data" + "/" + "SubH_NHOCH3_wb97xd_aug-cc-pvtz_reor"
    )
    cc_props = qt.get_cc_props(file_with_relpath, "N1")
    assert isinstance(cc_props, dict)
