Property Extraction Usage
=========================

Get xyz of atoms in sumfile
---------------------------

.. autofunction:: subproptools.qtaim_extract.get_xyz

Extracting Atomic Properties
----------------------------

To retrieve a dictionary of the atomic properties, you can use the ``subproptools.qtaim_extract.get_atomic_props()`` function

.. autofunction:: subproptools.qtaim_extract.get_atomic_props


Convert Atomic Properties to Group Properties
---------------------------------------------

.. autofunction:: subproptools.qtaim_extract.get_sub_props

Extracting BCP Properties
-------------------------
To retrieve a dictionary of BCP properties for a specific BCP:

.. autofunction:: subproptools.qtaim_extract.get_bcp_properties

To get a dictionary of BCP properties for a set of requested BCPs:

.. autofunction:: subproptools.qtaim_extract.extract_requested_bcp_props

To find a list of all BCPs in file in a format to use for the above functions:

.. autofunction:: subproptools.qtaim_extract.find_all_connections

Get DI between a substituent and the rest of the molecule
---------------------------------------------------------
.. autofunction:: subproptools.qtaim_extract.get_sub_di

Extracting Properties for a set of molecules
--------------------------------------------
Gets BCP and group properties

.. autofunction:: subproptools.qtaim_extract.sub_prop_frame

Charge Concentration Manipulation
---------------------------------

Find all CCs on an atom

.. autofunction:: subproptools.qtaim_extract.get_cc_props


Filter all CCs for VSCCs
.. autofunction:: subproptools.qtaim_extract.identify_vscc

Combine those two into one step

.. autofunction:: subproptools.qtaim_extract.get_atom_vscc

Get VSCC for a set of atoms

.. autofunction:: subproptools.qtaim_extract.extract_requested_cc_props

Get CC, BCP, Group and Atomic propeties for a substituent
---------------------------------------------------------
.. autofunction:: subproptools.qtaim_extract.extract_sub_props
