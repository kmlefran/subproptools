"""Tests for the sub_reor module of subproptools"""
import pandas as pd

from subproptools import sub_reor as sr


def test_rotate_cyclopropyl(filepath_tests):
    """Test rotation of cyclopropyl"""
    rel_name = str(filepath_tests) + "/" "test_data" + "/" + "cyclopropyl"
    rot_frame = sr.rotate_substituent(rel_name, 1, 2)
    assert isinstance(rot_frame, pd.DataFrame)
    assert len(rot_frame["Atom"]) == 9
    assert abs(list(rot_frame["x"])[0]) < 0.0001
    assert abs(list(rot_frame["y"])[0]) < 0.0001
    assert abs(list(rot_frame["z"])[0]) < 0.0001
    assert list(rot_frame["x"])[1] < 0.1
    assert abs(list(rot_frame["y"])[1]) < 0.0001
    assert abs(list(rot_frame["z"])[1]) < 0.0001


def test_rotate_specified_atom(filepath_tests):
    """Test that specifying an atom to be on +y works"""
    rel_name = str(filepath_tests) + "/" "test_data" + "/" + "CFHCH3"
    rot_frame = sr.rotate_substituent(rel_name, 1, 2, 4)
    assert isinstance(rot_frame, pd.DataFrame)
    assert len(rot_frame["Atom"]) == 8
    assert abs(list(rot_frame["x"])[0]) < 0.0001
    assert abs(list(rot_frame["y"])[0]) < 0.0001
    assert abs(list(rot_frame["z"])[0]) < 0.0001
    assert list(rot_frame["x"])[1] < 0.1
    assert abs(list(rot_frame["y"])[1]) < 0.0001
    assert abs(list(rot_frame["z"])[1]) < 0.0001
    assert list(rot_frame["x"])[3] > 0.1
    print(rot_frame)
    assert list(rot_frame["y"])[3] > 0.1

    assert abs(list(rot_frame["z"])[3]) < 0.01


def test_logical_match(filepath_tests):
    """Test that the match found is logical. Here that means that F lines up with F, H with H and CH3 with CH3

    Test relates by analyzing the reference_maps.py module geomtry.

    e.g. in the reference the Fluorine is on this axis, so that should be here too
    """
    rel_name = str(filepath_tests) + "/" "test_data" + "/" + "CFHCH3"
    rot_frame = sr.rotate_substituent(rel_name, 1, 3)
    assert isinstance(rot_frame, pd.DataFrame)
    assert abs(list(rot_frame["x"])[0]) < 0.0001
    assert abs(list(rot_frame["y"])[0]) < 0.0001
    assert abs(list(rot_frame["z"])[0]) < 0.0001
    assert list(rot_frame["x"])[2] < 0.1
    assert abs(list(rot_frame["y"])[2]) < 0.0001
    assert abs(list(rot_frame["z"])[2]) < 0.0001
    assert list(rot_frame["y"])[3] < -0.1
    assert list(rot_frame["z"])[3] < -0.1
    assert list(rot_frame["y"])[4] < -0.1
    assert list(rot_frame["z"])[4] > 0.1
    assert list(rot_frame["y"])[1] > 0.1
    assert abs(list(rot_frame["z"])[1]) < 0.01


def test_rotate_hydroxyl(filepath_tests):
    """Test the rotation of a hydroxyl group"""
    rel_name = (
        str(filepath_tests) + "/" "test_data" + "/" + "water_wb97xd_augccpvtz_qtaim"
    )
    rot_frame = sr.rotate_substituent(rel_name, 1, 2)
    # Test successful rotation - data frame and right number of atoms
    assert isinstance(rot_frame, pd.DataFrame)
    assert len(list(rot_frame["x"])) == 3
    # Test 1st atom at 0,0,0
    assert abs(rot_frame["x"][0]) < 0.01
    assert abs(rot_frame["y"][0]) < 0.01
    assert abs(rot_frame["z"][0]) < 0.01
    # Test 2nd atom on -x axis
    assert abs(rot_frame["y"][1]) < 0.01
    assert rot_frame["x"][1] < -0.1
    assert abs(rot_frame["z"][1]) < 0.01
    # Test 3rd atom in xy plane with +x and -y  (average LPs oriented to +y)
    assert rot_frame["x"][2] > 0.1
    assert rot_frame["y"][2] < -0.1
    assert abs(rot_frame["z"][2]) < 0.01
