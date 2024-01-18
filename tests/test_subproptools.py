"""Tests to ensure the subproptools module works as intended"""
# import math  # sqrt
# import os  # file system stuff


import numpy as np  # arrays
import pandas as pd  # data frames

from subproptools import qtaim_extract as qt

# import pkg_resources


# from subproptools import qtaim_extract as qt


def load_sumfile(filepath_tests, filename):
    """Test loading a file"""
    fname = filepath_tests / "test_data" / filename
    with open(fname, encoding="utf-8") as f:
        data = f.readlines()
    return data


def test_atomic_properties(filepath_tests):
    """Test building atomic properties"""
    data = load_sumfile(
        filepath_tests, filename="SubH_NHOCH3_wb97xd_aug-cc-pvtz_reor.sum"
    )
    atomic_properties = qt.get_atomic_props(data)
    assert isinstance(atomic_properties, dict)
    assert len(atomic_properties.keys()) == 8
    assert "N1" in atomic_properties
    assert (
        "q" in atomic_properties["N1"]
    )  # pylint:disable=consider-iterating-dictionary


def test_bcp_properties(filepath_tests):
    """Test extraction of bcp properties"""
    data = load_sumfile(
        filepath_tests, filename="SubH_NHOCH3_wb97xd_aug-cc-pvtz_reor.sum"
    )
    bcp_properties = qt.get_bcp_properties(data, ["N1", "H2"])
    assert isinstance(bcp_properties, dict)


def test_ldm(filepath_tests):
    """Test constuction of the LDM matrix"""
    data = load_sumfile(
        filepath_tests, filename="SubH_NHOCH3_wb97xd_aug-cc-pvtz_reor.sum"
    )
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


def test_get_sub_di(filepath_tests):
    """Test that getting substituent di returns a float"""
    data = load_sumfile(
        filepath_tests, filename="SubH_NHOCH3_wb97xd_aug-cc-pvtz_reor.sum"
    )
    sub_di = qt.get_sub_di(data, ["N1", "H3", "O4", "C5"])
    assert isinstance(sub_di, float)
    assert sub_di > 0


def test_no_vscc_in_ammonium(filepath_tests):
    """Test filtering the charge concentrations for VSCC"""
    file_with_relpath = (
        filepath_tests
        / "test_data"
        / "SubH_NHOCH3_wb97xd_aug-cc-pvtz_reor_atomicfiles"
        / "n1.agpviz"
    )
    with open(file_with_relpath, encoding="utf-8") as f:
        data = f.readlines()
    cc_props = qt.get_cc_props(data, "N1", is_lines_data=True)
    data = load_sumfile(
        filepath_tests, filename="SubH_NHOCH3_wb97xd_aug-cc-pvtz_reor.sum"
    )
    atomic_properties = qt.get_atomic_props(data)
    vscc_props = qt.identify_vscc(cc_props, atomic_properties, "N1")
    assert isinstance(vscc_props, dict)
    assert len(vscc_props) == 0


def test_get_atom_vscc(filepath_tests):
    """Test get_atom_vscc"""
    file_with_relpath = (
        filepath_tests
        / "test_data"
        / "SubH_NHOCH3_wb97xd_aug-cc-pvtz_reor_atomicfiles"
        / "o4.agpviz"
    )
    with open(file_with_relpath, encoding="utf-8") as f:
        data = f.readlines()
    data2 = load_sumfile(
        filepath_tests, filename="SubH_NHOCH3_wb97xd_aug-cc-pvtz_reor.sum"
    )
    atomic_properties = qt.get_atomic_props(data2)
    vscc = qt.get_atom_vscc(data, "O4", atomic_properties, True)
    assert isinstance(vscc, dict)
    assert len(vscc) == 3


def test_get_sub_props(filepath_tests):
    """Test get_sub_props"""
    data = load_sumfile(
        filepath_tests, filename="SubH_NHOCH3_wb97xd_aug-cc-pvtz_reor.sum"
    )
    atomic_properties = qt.get_atomic_props(data)
    sub_props = qt.get_sub_props(
        atomic_properties,
        [1, 3, 4, 5, 6, 7, 8],
        ["N1", "H2", "H3", "O4", "C5", "H6", "H7", "H8"],
    )
    assert isinstance(sub_props, dict)


def test_extract_sub_props(filepath_tests):
    """Text extract_sub_props"""
    data = load_sumfile(
        filepath_tests, filename="SubH_NHOCH3_wb97xd_aug-cc-pvtz_reor.sum"
    )
    sub_props = qt.extract_sub_props(
        data,
        [1, 3, 4, 5, 6, 7, 8],
        str(filepath_tests)
        + "/"
        + "test_data"
        + "/"
        + "SubH_NHOCH3_wb97xd_aug-cc-pvtz_reor",
        lapRhoCpAtoms=[4],
    )
    assert isinstance(sub_props, dict)
    assert "Group" in sub_props
    assert "Atomic" in sub_props
    assert "BCP" in sub_props
    assert "VSCC" in sub_props


def test_extract_requested_cc_props(filepath_tests):
    """Test extract_requested_cc_props"""
    data = load_sumfile(
        filepath_tests, filename="SubH_NHOCH3_wb97xd_aug-cc-pvtz_reor.sum"
    )
    atomic_properties = qt.get_atomic_props(data)
    vscc_props = qt.extract_requested_cc_props(
        [1, 4],
        str(filepath_tests)
        + "/"
        + "test_data"
        + "/"
        + "SubH_NHOCH3_wb97xd_aug-cc-pvtz_reor",
        ["N1", "H2", "H3", "O4", "C5", "H6", "H7", "H8"],
        atomic_properties,
    )
    assert len(vscc_props) == 2
    assert len(vscc_props["N1"]) == 0
    assert len(vscc_props["O4"]) == 3


def test_get_selected_bcps(filepath_tests):
    """Test get_selected_bcps"""
    data = load_sumfile(
        filepath_tests, filename="SubH_NHOCH3_wb97xd_aug-cc-pvtz_reor.sum"
    )
    bcp_dict = qt.get_selected_bcps(data, [["N1", "H2"], ["N1", "H3"]])
    assert isinstance(bcp_dict, dict)
    assert len(bcp_dict) == 2
    assert "N1-H2" in bcp_dict
    assert "N1-H3" in bcp_dict


def test_find_connected(filepath_tests):
    """Test find_connected"""
    data = load_sumfile(
        filepath_tests, filename="SubH_NHOCH3_wb97xd_aug-cc-pvtz_reor.sum"
    )
    con = qt.find_connected(data, "H2", "N1")
    assert isinstance(con, list)
    assert len(con) == 3


def test_find_all_connected(filepath_tests):
    """Test find_all_connected"""
    data = load_sumfile(
        filepath_tests, filename="SubH_NHOCH3_wb97xd_aug-cc-pvtz_reor.sum"
    )
    bcps = qt.find_all_connections(data)
    assert isinstance(bcps, list)
    assert len(bcps) == 7
    assert len(bcps[0]) == 2
    assert isinstance(bcps[0][0], str)


def test_get_xyz(filepath_tests):
    """Test get_xyz"""
    data = load_sumfile(
        filepath_tests, filename="SubH_NHOCH3_wb97xd_aug-cc-pvtz_reor.sum"
    )
    xyz = qt.get_xyz(data)
    assert isinstance(xyz, dict)
    assert isinstance(xyz["xyz"], (np.ndarray, np.generic))
    assert isinstance(xyz["Atoms"], list)
    assert len(xyz["Atoms"]) == 8
    assert len(xyz["xyz"]) == 8


def test_sub_prop_frame(filepath_tests):
    """test sub_prop_frame"""
    df = pd.DataFrame(
        {
            "Substituent": ["NHOCH3"],
            "subAtoms": [1, 3, 4, 5, 6, 7, 8],
            "SubH": [
                str(filepath_tests)
                + "/"
                + "test_data"
                + "/"
                + "SubH_NHOCH3_wb97xd_aug-cc-pvtz_reor"
            ],
        }
    )
    df.to_csv(str(filepath_tests) + "/" + "test_data" + "/" + "testcsv.csv")
    fr = qt.sub_prop_frame(
        str(filepath_tests) + "/" + "test_data" + "/" + "testcsv.csv"
    )
    assert isinstance(fr, dict)
    assert "SubH" in fr
    assert "Group" in fr["SubH"]
    assert isinstance(fr["SubH"]["Group"], pd.DataFrame)
    assert "BCP" in fr["SubH"]
    assert isinstance(fr["SubH"]["BCP"], pd.DataFrame)


def test_extract_requested_bcp_props(filepath_tests):
    """Test extract_requested_bcp_props"""
    data = load_sumfile(
        filepath_tests, filename="SubH_NHOCH3_wb97xd_aug-cc-pvtz_reor.sum"
    )
    atomic_properties = qt.get_atomic_props(data)
    bcp_props = qt.extract_requested_bcp_props(
        data,
        ["N1", "H2", "H3", "O4", "C5", "H6", "H7", "H8"],
        [[1, 2], [1, 3]],
        ["N1", "H3", "O4", "C5", "H6", "H7", "H8"],
        atomic_properties,
    )
    assert isinstance(bcp_props, dict)
    assert "N1-H2" in bcp_props
    assert "N1-H3" in bcp_props


def test_get_cc_props(filepath_tests):
    """Test that getting charge concentration properties returns a dict"""
    file_with_relpath = (
        filepath_tests
        / "test_data"
        / "SubH_NHOCH3_wb97xd_aug-cc-pvtz_reor_atomicfiles"
        / "n1.agpviz"
    )
    with open(file_with_relpath, encoding="utf-8") as f:
        data = f.readlines()
    cc_props = qt.get_cc_props(data, "N1", is_lines_data=True)
    assert isinstance(cc_props, dict)
    assert len(cc_props) == 5
    assert "xyz" in cc_props[1]
    assert "rho" in cc_props[1]
    assert "delsqrho" in cc_props[1]
    assert "distFromNuc" in cc_props[1]
