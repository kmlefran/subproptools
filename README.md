
[![ci](https://github.com/kmlefran/subproptools/actions/workflows/ci.yml/badge.svg)](https://github.com/kmlefran/subproptools/actions/workflows/ci.yml)
[![Documentation Status](https://readthedocs.org/projects/subproptools/badge/?version=latest)](https://subproptools.readthedocs.io/en/latest/?badge=latest)
[![Coverage Status](https://coveralls.io/repos/github/kmlefran/subproptools/badge.svg?branch=main)](https://coveralls.io/github/kmlefran/subproptools?branch=main)
[![PyPI version](https://badge.fury.io/py/subproptools.svg)](https://badge.fury.io/py/subproptools)
# Molecule Reorientation and .sum file handling

This package provides tools to extract data from .sum files (in the qtaimExtract module), and tools to reorient the molecule to defined coordinate systems (in the subreor module). The subreor module rotates to a defined coordinate system. First an atom is positioned at the origin. By convention, this is the atom of the group that is directly bonded to the rest of the molecule. The rest of the molecule is placed along the -x axis. The remaining axes are defined as follows:
* if there is one lone pair(VSCC), that point lies on the +y
* if there are two lone pairs, the average position of them lies on the +y
* if there are no lone pairs, map the BCPs of the atom at the origin to a reference. Identify the closest match to the reference to determine a BCP to set as +y

# Authors
Kevin Lefrancois-Gagnon
Robert C. Mawhinney

# Installation prior to distribution
```
pip install git+https://github.com/kmlefran/subproptools
```

# User Facing Functions:
## subreor
* rotate_sheet - performs rotations to the desired coordinate system for many molecules as defined in a csv file
* output_to_gjf - writes the reoreinted geometry to a .gjf file
* rotate_substituent - performs the defined rotation for an individual substituent

## qtaimExtract
* get_sub_di - gets delocalization index between a substituent and the rest of the molecule
* get bcp_properties - returns a dictionary of bcp properties for a given bcp
* get_atomic_props - returns dictionary of atomic properties for all atoms in molecule
* get_cc_props - returns dictionary of properties for all charge concentrations
* identify_vscc - returns only the VSCC from a set of charge concentrations
* get_sub_props - returns dictionary of one substituent's properties
* extract_sub_props - returns dictionary of atomic, group, bcp and vscc properties for a substituent
* sub_prop_frame - returns dictioanry of frames of group properties
* get_xyz - extracts a molecule's xyz geometry from sum file
