"""qtaimExtract
Tools for extracting data from sum files(and CC from agpviz), and calculating substituent properties specifically
"""
import ast
import math  # sqrt
import os  # file system stuff

import numpy as np  # arrays
import pandas as pd  # data frames


def _search_str(linesObj, word, searchStart=0, ignore=0):
    """Given lines of file, return line that word is on.

    Args:
        lines obj - list, each containing a string. e.g. the output of file.readlines() method
        word - string that we look to match
        searchStart - integer, if we don't need to start searching at the beginning, choose \
        line number to start at
        ignore - integer on how many instances of the string to ignore
    Returns:
        index of line that word occurs, or -1 if not found
    """
    wordLine = -1  # will return -1 if string not found
    for ln_num, line in enumerate(linesObj):  # iterate over lines
        if line.find(word) > -1 and linesObj.index(line) >= searchStart:
            if ignore == 0:
                wordLine = ln_num
                break
            ignore = ignore - 1
    return wordLine


def _get_di_table(
    data,
):
    """Given lines of sum file, returns the table containing DI information.
    Args:
        data: list[str] - the lines of a sum file
    Returns:
        pandas DataFrame with columns:
        Atom A, Atom B, 2*D2(A,B), DI(A,B), %Deloc(A,B), %Deloc(B,A)
    """
    tableStart = (
        _search_str(
            data,
            "2*D2(A,B)          DI(A,B)          %Deloc(A,B)       %Deloc(B,A)",
            searchStart=0,
        )
        + 2
    )  # data starts 2 lines after header
    tableEnd = (
        _search_str(data, "Areas of Interatomic Surfaces:", tableStart) - 2
    )  # -2 since there is a 2 line gap
    tableData = data[tableStart:tableEnd]
    headerNames = ["Atom A", "Atom B"]
    for (
        colName
    ) in "2*D2(A,B)          DI(A,B)          %Deloc(A,B)       %Deloc(B,A)".split():
        headerNames.append(colName)
    table = pd.DataFrame(columns=headerNames)
    # Add each row of table to Data Frame
    for i, line in enumerate(tableData):
        table.loc[i] = line.split()[0:6]
    return table


def _get_table(data, tableHeader, ignored=0, endString="Total"):
    """Given lines of sum file and table header, return pandas dataframe of requested table

    Args:
        data:list[str] - lines of sum file
        tableHeader - string containing non-atom columns of table to be found
        ignored - int - find ignored+1th occurence
        endString - the string found 2 lines after the table data ends
    Returns:
        pandas dataframe of requested table
    """
    tableStart = (
        _search_str(data, tableHeader, searchStart=0, ignore=ignored) + 2
    )  # data starts 2 lines after header
    tableEnd = (
        _search_str(data, endString, tableStart) - 1
    )  # only deduct 1, because next line will get up to but not including that line
    tableData = data[tableStart:tableEnd]  # get data
    # Include atom here, as sumfile header is Atom A and the whitespace there will mess with how code is writting
    headerNames = ["Atom"]
    for colName in tableHeader.split():
        headerNames.append(colName)
    table = pd.DataFrame(columns=headerNames)
    n_col = len(headerNames)
    for i, line in enumerate(tableData):
        split_line = line.split()
        if "Vol" in tableHeader and len(split_line) < n_col:
            split_line.insert(1, 0.0)
        table.loc[i] = split_line

    return table


def _get_bcp_block(data, atPair=["C1", "H2"]):  # pylint:disable=dangerous-default-value
    """Given sum file lines and pair of atoms, return lines in file corresponding to that BCP."""
    bcpStart = _search_str(data, word=f"Type = (3,-1) BCP {atPair[0]} {atPair[1]}") - 1
    if bcpStart == -2:
        bcpStart = (
            _search_str(data, word=f"Type = (3,-1) BCP {atPair[1]} {atPair[0]}") - 1
        )
    if bcpStart == -2:
        raise ValueError(f"BCP between {atPair[0]} and {atPair[1]} not found in file")
    bcpEnd = bcpStart + 34
    bcpBlock = data[bcpStart:bcpEnd]
    return bcpBlock  # return the lines of the BCP data


def get_sub_di(data, subAtomLabels):
    """Given lines of sum file and labels of atoms in substituent, return DI between substituent and rest of molecule"""
    diTable = _get_di_table(data)
    diTable = diTable.drop(["2*D2(A,B)", "%Deloc(A,B)", "%Deloc(B,A)"], axis=1)
    diClass = []

    diTable = diTable.astype({"DI(A,B)": float})
    for ind in diTable.index:
        if any(x == diTable.at[ind, "Atom A"] for x in subAtomLabels) and any(
            x == diTable.at[ind, "Atom B"] for x in subAtomLabels
        ):
            diClass.append("Substituent-Substituent")
        elif any(x == diTable.at[ind, "Atom A"] for x in subAtomLabels) or any(
            x == diTable.at[ind, "Atom B"] for x in subAtomLabels
        ):
            diClass.append("Substituent-Substrate")
        else:
            diClass.append("Substrate-Substrate")
    diTable.insert(2, "Interaction", diClass, True)
    return sum(
        diTable.loc[diTable["Interaction"] == "Substituent-Substrate", "DI(A,B)"]
    )


def get_bcp_properties(
    data, atPair=["C1", "H2"]
):  # pylint:disable=dangerous-default-value
    """Given lines of sum file and atom pair, return bcp properties for pair

    Args:
        data: list[str] - lines of a .sum file
        atPair: list of atom labels for which to find BCP, e.g. ['C1','H2']
    Returns:
        dictionary containing properties of atPair bcp
        keys are:
        Coords (np array of xyz coords), 'Rho', lambda1,lambda2,lambda3,DelSqRho,Ellipticity,V,G,H
    """
    bcpBlock = _get_bcp_block(data, atPair)
    bcpDict = {}
    for line in bcpBlock:  # iterate over lines
        splitLine = (
            line.split()
        )  # split the line into individual words based on whitespace
        if "Coords" in splitLine:  # if line contains the coordinates, get coordinates
            bcpDict.update(
                {
                    "xyz": np.array(
                        [float(splitLine[4]), float(splitLine[5]), float(splitLine[6])]
                    )
                }
            )
        elif "Rho" in splitLine:  # get density at BCP
            bcpDict.update({"rho": [float(splitLine[2])]})
        elif "HessRho_EigVals" in splitLine:  # get lambda_i
            bcpDict.update({"lambda1": [float(splitLine[2])]})
            bcpDict.update({"lambda2": [float(splitLine[3])]})
            bcpDict.update({"lambda3": [float(splitLine[4])]})
        elif "DelSqRho" in splitLine:  # get DelSqRho
            bcpDict.update({"DelSqRho": [float(splitLine[2])]})
        elif "Ellipticity" in splitLine:  # get ellipticity
            bcpDict.update({"Ellipticity": [float(splitLine[3])]})
        elif "V" in splitLine:  # get V
            bcpDict.update({"V": [float(splitLine[2])]})
        elif "G" in splitLine:  # get G
            bcpDict.update({"G": [float(splitLine[2])]})
    bcpDict.update(
        {"H": [bcpDict["V"][0] + bcpDict["G"][0]]}
    )  # get total energy density
    return bcpDict


def get_atomic_props(data):
    """Returns a dictionary of atomic properties

    :param data: The string lines of a .sum file
    :type data: list[str]
    :return: Dicionary with a key for each atom in the molecule
    :rtype: dict

    """
    xyzFrame = _get_table(
        data,
        "Charge                X                  Y                  Z",
        endString="Some Atomic Properties:",
    )
    eFrame = _get_table(
        data,
        "q(A)              L(A)              K(A)          K_Scaled(A)      |Mu_Intra(A)|",
    )
    muXFrame = _get_table(data, "Mu_Intra_X(A)     Mu_Bond_X(A)        Mu_X(A)")
    muYFrame = _get_table(data, "Mu_Intra_Y(A)     Mu_Bond_Y(A)        Mu_Y(A)")
    muZFrame = _get_table(data, "Mu_Intra_Z(A)     Mu_Bond_Z(A)        Mu_Z(A)")
    muMagFrame = _get_table(data, "|Mu_Intra(A)|     |Mu_Bond(A)|        |Mu(A)|")
    quadFrame = _get_table(
        data,
        "Q_XX(A)           Q_XY(A)           Q_XZ(A)           Q_YY(A)           Q_YZ(A)           Q_ZZ(A)",
        endString="Eigenvalues and Eigenvectors of Atomic Traceless Quadrupole Moment Tensors",
    )
    radFrame = _get_table(
        data,
        "R-2(A)            R-1(A)            R0(A)             R+1(A)            R+2(A)",
        endString="Atomic Radial Distortion Moments",
    )
    volFrame = _get_table(
        data,
        "Area(A)           Vol(A)          N(Vol(A))      N(Vol(A))/Vol(A)     %N(Vol(A))",
        ignored=1,
    )
    allatomDict = {}  # initialize empty dictionary to store dictionaries of the atoms
    for atom in eFrame["Atom"]:
        atomDict = {}  # create a dictionary for each atom
        atomDict.update(
            {
                "xyz": np.array(
                    [
                        float(xyzFrame[xyzFrame["Atom"] == atom]["X"].iloc[0]),
                        float(xyzFrame[xyzFrame["Atom"] == atom]["Y"].iloc[0]),
                        float(xyzFrame[xyzFrame["Atom"] == atom]["Z"].iloc[0]),
                    ]
                )
            }
        )
        atomDict.update({"q": [float(eFrame[eFrame["Atom"] == atom]["q(A)"].iloc[0])]})
        atomDict.update({"K": [float(eFrame[eFrame["Atom"] == atom]["K(A)"].iloc[0])]})
        atomDict.update(
            {"K_Scaled": [float(eFrame[eFrame["Atom"] == atom]["K_Scaled(A)"].iloc[0])]}
        )
        atomDict.update(
            {
                "Mu_Intra_X": [
                    float(muXFrame[muXFrame["Atom"] == atom]["Mu_Intra_X(A)"].iloc[0])
                ]
            }
        )
        atomDict.update(
            {
                "Mu_Intra_Y": [
                    float(muYFrame[muYFrame["Atom"] == atom]["Mu_Intra_Y(A)"].iloc[0])
                ]
            }
        )
        atomDict.update(
            {
                "Mu_Intra_Z": [
                    float(muZFrame[muZFrame["Atom"] == atom]["Mu_Intra_Z(A)"].iloc[0])
                ]
            }
        )
        atomDict.update(
            {
                "Mu_Bond_X": [
                    float(muXFrame[muXFrame["Atom"] == atom]["Mu_Bond_X(A)"].iloc[0])
                ]
            }
        )
        atomDict.update(
            {
                "Mu_Bond_Y": [
                    float(muYFrame[muYFrame["Atom"] == atom]["Mu_Bond_Y(A)"].iloc[0])
                ]
            }
        )
        atomDict.update(
            {
                "Mu_Bond_Z": [
                    float(muZFrame[muZFrame["Atom"] == atom]["Mu_Bond_Z(A)"].iloc[0])
                ]
            }
        )
        atomDict.update(
            {"Mu_X": [float(muXFrame[muXFrame["Atom"] == atom]["Mu_X(A)"].iloc[0])]}
        )
        atomDict.update(
            {"Mu_Y": [float(muYFrame[muYFrame["Atom"] == atom]["Mu_Y(A)"].iloc[0])]}
        )
        atomDict.update(
            {"Mu_Z": [float(muZFrame[muZFrame["Atom"] == atom]["Mu_Z(A)"].iloc[0])]}
        )
        atomDict.update(
            {
                "|Mu_Intra|": [
                    float(
                        muMagFrame[muMagFrame["Atom"] == atom]["|Mu_Intra(A)|"].iloc[0]
                    )
                ]
            }
        )
        atomDict.update(
            {
                "|Mu_Bond|": [
                    float(
                        muMagFrame[muMagFrame["Atom"] == atom]["|Mu_Bond(A)|"].iloc[0]
                    )
                ]
            }
        )
        atomDict.update(
            {"|Mu|": [float(muMagFrame[muMagFrame["Atom"] == atom]["|Mu(A)|"].iloc[0])]}
        )
        atomDict.update(
            {"Q_XX": [float(quadFrame[quadFrame["Atom"] == atom]["Q_XX(A)"].iloc[0])]}
        )
        atomDict.update(
            {"Q_XY": [float(quadFrame[quadFrame["Atom"] == atom]["Q_XY(A)"].iloc[0])]}
        )
        atomDict.update(
            {"Q_XZ": [float(quadFrame[quadFrame["Atom"] == atom]["Q_XZ(A)"].iloc[0])]}
        )
        atomDict.update(
            {"Q_YY": [float(quadFrame[quadFrame["Atom"] == atom]["Q_YY(A)"].iloc[0])]}
        )
        atomDict.update(
            {"Q_YZ": [float(quadFrame[quadFrame["Atom"] == atom]["Q_YZ(A)"].iloc[0])]}
        )
        atomDict.update(
            {"Q_ZZ": [float(quadFrame[quadFrame["Atom"] == atom]["Q_ZZ(A)"].iloc[0])]}
        )
        atomDict.update(
            {"R+2": [float(radFrame[radFrame["Atom"] == atom]["R+2(A)"].iloc[0])]}
        )
        atomDict.update(
            {"R+1": [float(radFrame[radFrame["Atom"] == atom]["R+1(A)"].iloc[0])]}
        )
        atomDict.update(
            {"Vol": [float(volFrame[volFrame["Atom"] == atom]["Vol(A)"].iloc[0])]}
        )
        atomDict.update({"quadContrib": _get_atomic_quad_contrib(atomDict)})
        allatomDict.update(
            {atom: atomDict}
        )  # add the atomic dictionary to the dictionary of dictionaries
    return allatomDict


def _get_dist(coordListA, coordListB):
    """Return cartesian distance between two 3 dimensional np arrays."""
    dist = math.sqrt(
        (coordListA[0] - coordListB[0]) ** 2
        + (coordListA[1] - coordListB[1]) ** 2
        + (coordListA[2] - coordListB[2]) ** 2
    )
    return dist


def _is_bonded_cc(ccXYZ, atomDict, originAtom):
    """Checks if a charge concentration is a bonded charge concentration"""
    atom_labels = atomDict.keys()
    on_a_line = False
    for a_1 in atom_labels:
        if a_1 != originAtom:
            if _is_on_line(atomDict[a_1]["xyz"], atomDict[originAtom]["xyz"], ccXYZ):
                on_a_line = True
                break

    return on_a_line


def identify_vscc(multiccDict, atomDict, originAtom, thresh=0.7):
    """Given dictionary of charge concentraion properties and atomic properties, identify vscc

    Args:
        multiccDict - dictionary of cc properties for all ccs of an atom
        atomDict - get_atomic_properties object
        threshold - distance between inner shell and outer shell ccs
        #(eg inner shell charge concentration is more than 0.7 au closer to nuclei than VSCC)
    Returns:
        subset of multiccDict correspondng to valence shell charge concentrations

    Note:
        VSCCs identified by: not being on nuclei, not being on line between atoms, and being the outermost CCs
    """

    vsccDict = {}  # initialize empty dictionary for cc
    nucDistList = []  # empty list to store distances
    potentialccList = (
        []
    )  # empty list to store keys for potential VSCC after eliminating some criteria
    for cc in multiccDict:  # for each cc
        if (
            multiccDict[cc]["distFromNuc"] > 0.1
        ):  # this will filter out the innermost cc that are on the nuclei
            ccXYZ = multiccDict[cc]["xyz"]  # get the xyz coords of the cc
            # nucPair = _find_closest_nuclei(ccXYZ,atomDict) #find the two closest atoms to the cc
            # if it is a bonded cc, it will be the nuclei to which it is bonded
            # check if the cc is on the bond between those two nuclei
            isBonded = _is_bonded_cc(ccXYZ, atomDict, originAtom)
            if (
                not isBonded
            ):  # if it is not on the line, it is potentially a VSCC, store it, and its distance
                nucDistList.append(multiccDict[cc]["distFromNuc"])
                potentialccList.append(cc)
    # at this point, non-bonded, non-core CC remain, so just need to check which CC are the outermost
    # do this by comparing to maximum distance from nuclei
    if len(nucDistList) > 0:
        outerShellDist = max(nucDistList)
    for cc in potentialccList:
        if (
            abs(
                multiccDict[cc]["distFromNuc"]
                - outerShellDist  # pylint:disable=used-before-assignment
            )
            < thresh
        ):
            # for given cc, if close (within thresh) to outermost, it is a VSCC, store it as such
            # So, going up to 3rd row for now. May work for 4th but I haven't checked
            # C/O/F easy - any not core not bonded should be VSCC
            # P/S/Si harder
            # in example data inner~0.25-0.3, outer ~ 1.2-1.3
            # arbitrary default threshold seems like 0.7 will work
            vsccDict.update({cc: multiccDict[cc]})
    return vsccDict


def get_atom_vscc(filename, atomLabel, atomicProps, is_lines_data=False):
    """Returns a dicitonary of the properties of VSCC for the atom"""
    all_cc_dict = get_cc_props(filename, atomLabel, is_lines_data)
    vscc_dict = identify_vscc(all_cc_dict, atomicProps, atomLabel)
    return vscc_dict


def _is_on_line(lineStartList, lineEndList, pointToCheckList, epsilon=0.1):
    """Given 3d coords for start of line, end of line and point to check, sees if pointToCheck is on line"""
    # line connecting atoms: (x1 + t(x2-x1),y1 + t(y2-y1),z1 + t(z2-z1))
    # to reconstruct - create equations of lines. (pointToCheck-pointWitht) and (lineStart-pointWitht)
    # dot product of those vectors should be 0(closest point is. 90 degrees)
    # rearrange to solve for t
    t = (
        (pointToCheckList[0] - lineStartList[0]) * (lineEndList[0] - lineStartList[0])
        + (pointToCheckList[1] - lineStartList[1]) * (lineEndList[1] - lineStartList[1])
        + (pointToCheckList[2] - lineStartList[2]) * (lineEndList[2] - lineStartList[2])
    ) / (
        ((lineEndList[2] - lineStartList[2]) ** 2)
        + ((lineEndList[1] - lineStartList[1]) ** 2)
        + ((lineEndList[0] - lineStartList[0]) ** 2)
    )
    perpPoint = np.array(
        [
            lineStartList[0] + t * (lineEndList[0] - lineStartList[0]),
            lineStartList[1] + t * (lineEndList[1] - lineStartList[1]),
            lineStartList[2] + t * (lineEndList[2] - lineStartList[2]),
        ]
    )
    # get distance between test point and closest line point
    lowX = min(lineStartList[0], lineEndList[0])
    highX = max(lineStartList[0], lineEndList[0])
    lowY = min(lineStartList[1], lineEndList[1])
    highY = max(lineStartList[1], lineEndList[1])
    lowZ = min(lineStartList[2], lineEndList[2])
    highZ = max(lineStartList[2], lineEndList[2])
    # if t is between the end points of the line, and the pointToCheck is within epsilon of that point, it is a bonded cc
    in_x = perpPoint[0] >= lowX and perpPoint[0] <= highX
    in_y = perpPoint[1] >= lowY and perpPoint[1] <= highY
    in_z = perpPoint[2] >= lowZ and perpPoint[2] <= highZ
    if in_x and in_y and in_z and _get_dist(pointToCheckList, perpPoint) < epsilon:
        return True
    return False


def get_cc_props(filename, atomLabel, is_lines_data=False):
    """takes a sumfilename with no extension and an atom label that we want cc for (eg F4)
    returns dictionary of all cc for atomLabelwith all the (3,+3) cc and properties

    Args:
        filename: sum file name with no extension
        atomLabel: label of atom that we wish to find VSCC for, e.g. 'N3'

    Returns:
        Dictionary with one nested dictinoary for each VSCC found
        {1: {VSCC1 props},
        2::{VSCC2 props}
        ,...}
        Sub-dictionary keys include: xyz, rho, delsqrho, distFromNuc
    """
    if not is_lines_data:
        # example dir name SubH_CCF-ReorPosY-B3LYP-def2-TZVPPD-Field_atomicfiles
        pathToSubdir = os.getcwd() + "/" + filename + "_atomicfiles" + "/"
        with open(
            pathToSubdir + atomLabel.lower() + ".agpviz", encoding="utf-8"
        ) as atFile:
            atData = atFile.readlines()

    else:
        atData = filename
    allccDict = {}
    alldictcounter = 1  # counter since the key will be numeric

    for lnNum, line in enumerate(atData):
        if "Type = (3,+3)" in line:  # if we're looking at CC
            oneccDict = {}  # create empty dict for storing properties for one cc
            xyzSplit = atData[
                lnNum + 1
            ].split()  # line after label is xyz, split it and store xyz
            oneccDict.update(
                {
                    "xyz": np.array(
                        [float(xyzSplit[2]), float(xyzSplit[3]), float(xyzSplit[4])]
                    )
                }
            )
            # next line is distance from nuc, split and store
            distFromNucSplit = atData[lnNum + 2].split()
            oneccDict.update({"distFromNuc": float(distFromNucSplit[2])})
            # next is rho, then delsqrho - split and store those
            rhoSplit = atData[lnNum + 3].split()
            oneccDict.update({"rho": float(rhoSplit[2])})
            delsqrhoSplit = atData[lnNum + 4].split()
            oneccDict.update({"delsqrho": float(delsqrhoSplit[2])})
            allccDict.update({alldictcounter: oneccDict})
            alldictcounter += 1
    return allccDict


def _get_atomic_quad_contrib(atomDict):
    """Given atomDict from get_atomic_properties, create dictionary containing atomic quadupole contributions

    Note: formula used from paper by Laidig in 1991
    https://www.sciencedirect.com/science/article/abs/pii/000926149180247U
    """
    secondMomentDict = {}
    secondMomentDict.update(
        {
            "Qxx": atomDict["q"][0] * (atomDict["xyz"][0] ** 2)
            + (atomDict["Q_XX"][0] + atomDict["R+2"][0]) / 3
            + atomDict["xyz"][0] * atomDict["Mu_Intra_X"][0]
            + atomDict["xyz"][0] * atomDict["Mu_Intra_X"][0]
        }
    )
    secondMomentDict.update(
        {
            "Qxy": atomDict["q"][0] * atomDict["xyz"][0] * atomDict["xyz"][1]
            + atomDict["Q_XY"][0] / 3
            + atomDict["xyz"][1] * atomDict["Mu_Intra_X"][0]
            + atomDict["xyz"][0] * atomDict["Mu_Intra_Y"][0]
        }
    )
    secondMomentDict.update(
        {
            "Qxz": atomDict["q"][0] * atomDict["xyz"][0] * atomDict["xyz"][2]
            + atomDict["Q_XZ"][0] / 3
            + atomDict["xyz"][2] * atomDict["Mu_Intra_X"][0]
            + atomDict["xyz"][0] * atomDict["Mu_Intra_Z"][0]
        }
    )
    secondMomentDict.update(
        {
            "Qyy": atomDict["q"][0] * (atomDict["xyz"][1] ** 2)
            + (atomDict["Q_YY"][0] + atomDict["R+2"][0]) / 3
            + atomDict["xyz"][1] * atomDict["Mu_Intra_Y"][0]
            + atomDict["xyz"][1] * atomDict["Mu_Intra_Y"][0]
        }
    )
    secondMomentDict.update(
        {
            "Qyz": atomDict["q"][0] * atomDict["xyz"][1] * atomDict["xyz"][2]
            + atomDict["Q_YZ"][0] / 3
            + atomDict["xyz"][2] * atomDict["Mu_Intra_Y"][0]
            + atomDict["xyz"][1] * atomDict["Mu_Intra_Z"][0]
        }
    )
    secondMomentDict.update(
        {
            "Qzz": atomDict["q"][0] * (atomDict["xyz"][2] ** 2)
            + (atomDict["Q_ZZ"][0] + atomDict["R+2"][0]) / 3
            + atomDict["xyz"][2] * atomDict["Mu_Intra_Z"][0]
            + atomDict["xyz"][2] * atomDict["Mu_Intra_Z"][0]
        }
    )
    secondMomentDict.update(
        {
            "trace": secondMomentDict["Qxx"]
            + secondMomentDict["Qyy"]
            + secondMomentDict["Qzz"]
        }
    )
    atomicQuadrupoleDict = {}
    atomicQuadrupoleDict.update(
        {"Q_xx": [0.5 * (3 * secondMomentDict["Qxx"] - secondMomentDict["trace"])]}
    )
    atomicQuadrupoleDict.update({"Q_xy": [0.5 * (3 * secondMomentDict["Qxy"])]})
    atomicQuadrupoleDict.update({"Q_xz": [0.5 * (3 * secondMomentDict["Qxz"])]})
    atomicQuadrupoleDict.update(
        {"Q_yy": [0.5 * (3 * secondMomentDict["Qyy"] - secondMomentDict["trace"])]}
    )
    atomicQuadrupoleDict.update({"Q_yz": [0.5 * (3 * secondMomentDict["Qyz"])]})
    atomicQuadrupoleDict.update(
        {"Q_zz": [0.5 * (3 * secondMomentDict["Qzz"] - secondMomentDict["trace"])]}
    )
    return atomicQuadrupoleDict  # return dictionary


def get_sub_props(atomDict, subAtoms, atomList):
    """Given atomic property dictionary and atoms to use, return dictionary of group properties.

    Args:
        atomDict: output from get_atomic_properties
        subAtoms: list[int] integer labels of atoms in substituent [1, 2,...]
        atomList: list[str] string list of atom labels, ['C1', 'H2'...]

    Returns:
        dictionary of group properties with the following keys: q,K, K_Scaled, Mu_Intra_X, Mu_Intra_Y,
        Mu_Intra_Z, Mu_Bond_X, Mu_Bond_Y,Mu_Bond_Z, Mu_X, Mu_Y, Mu_Z, Q_xx, Q_xy, Q_xz,
        Q_yy,Q_yz,Q_zz,Vol, |Mu_Intra|,|Mu_Bond|,|Mu|

    Note:
        in output dictionary, each property is stored as a one element list(to enable pandas later)
        Access properties as outDict['key'][0]
    """
    groupDict = {}  # create empty dict
    first = True  # flag for first iteration through loop
    for atom in subAtoms:
        if first:  # if first time, create all elements in dicionary
            first = False  # don't come here again
            if atomList[atom - 1] not in list(atomDict.keys()):
                raise ValueError(f"{atomList[atom-1]} not in atoms in file")
            groupDict = {
                "q": [atomDict[atomList[atom - 1]]["q"][0]],
                "K": [atomDict[atomList[atom - 1]]["K"][0]],
                "K_Scaled": [atomDict[atomList[atom - 1]]["K_Scaled"][0]],
                "Mu_Intra_X": [atomDict[atomList[atom - 1]]["Mu_Intra_X"][0]],
                "Mu_Intra_Y": [atomDict[atomList[atom - 1]]["Mu_Intra_Y"][0]],
                "Mu_Intra_Z": [atomDict[atomList[atom - 1]]["Mu_Intra_Z"][0]],
                "Mu_Bond_X": [atomDict[atomList[atom - 1]]["Mu_Bond_X"][0]],
                "Mu_Bond_Y": [atomDict[atomList[atom - 1]]["Mu_Bond_Y"][0]],
                "Mu_Bond_Z": [atomDict[atomList[atom - 1]]["Mu_Bond_Z"][0]],
                "Mu_X": [atomDict[atomList[atom - 1]]["Mu_X"][0]],
                "Mu_Y": [atomDict[atomList[atom - 1]]["Mu_Y"][0]],
                "Mu_Z": [atomDict[atomList[atom - 1]]["Mu_Z"][0]],
                "Q_xx": [atomDict[atomList[atom - 1]]["quadContrib"]["Q_xx"][0]],
                "Q_xy": [atomDict[atomList[atom - 1]]["quadContrib"]["Q_xy"][0]],
                "Q_xz": [atomDict[atomList[atom - 1]]["quadContrib"]["Q_xz"][0]],
                "Q_yy": [atomDict[atomList[atom - 1]]["quadContrib"]["Q_yy"][0]],
                "Q_yz": [atomDict[atomList[atom - 1]]["quadContrib"]["Q_yz"][0]],
                "Q_zz": [atomDict[atomList[atom - 1]]["quadContrib"]["Q_zz"][0]],
                "Vol": [atomDict[atomList[atom - 1]]["Vol"][0]],
            }
        else:  # for the rest, add the atomic property to the group dictionary element
            for prop in groupDict:  # pylint:disable=consider-using-dict-items
                if "Q" not in prop:
                    groupDict[prop][0] += atomDict[atomList[atom - 1]][prop][0]
                else:
                    groupDict[prop][0] += atomDict[atomList[atom - 1]]["quadContrib"][
                        prop
                    ][0]
    groupDict.update(
        {
            "|Mu_Intra|": [
                math.sqrt(
                    groupDict["Mu_Intra_X"][0] ** 2
                    + groupDict["Mu_Intra_Y"][0] ** 2
                    + groupDict["Mu_Intra_Z"][0] ** 2
                )
            ]
        }
    )
    groupDict.update(
        {
            "|Mu_Bond|": [
                math.sqrt(
                    groupDict["Mu_Bond_X"][0] ** 2
                    + groupDict["Mu_Bond_Y"][0] ** 2
                    + groupDict["Mu_Bond_Z"][0] ** 2
                )
            ]
        }
    )
    groupDict.update(
        {
            "|Mu_Bond|": [
                math.sqrt(
                    groupDict["Mu_X"][0] ** 2
                    + groupDict["Mu_Y"][0] ** 2
                    + groupDict["Mu_Z"][0] ** 2
                )
            ]
        }
    )
    return groupDict


def _check_num_atoms(atom_label_list: list[str], atom_int_list: list[str]) -> None:
    num_atoms = len(atom_label_list)
    max_int = max(atom_int_list)
    if max_int > num_atoms:
        raise ValueError(
            f"Largest group atom index {max_int} greater than number of atoms in file"
        )


def extract_sub_props(
    data: list[str],
    subAtoms: list[int],
    sumFileNoExt: str,
    groupProps=True,
    bcpId=[[1, 2]],
    lapRhoCpAtoms=[],
):  # pylint:disable=dangerous-default-value
    # pylint:disable=too-many-arguments
    """returns a dictionary of all group properties - bcp, group, atomic, and vscc
    Args:
        data: list[str] - lines of a .sum file
        subAtoms - indices of atoms in the molecule comprising group - starts at 1
        groupProps - boolean - do you want to compute group properties
        bcpId - list of 2 length lists. Pass empty list if no bcp properties wanted.
        2 length lists are indices of atoms that you want BCPs for
        lapRhoCpAtoms = list of atom indices that you want to find laprhoCPs for. Defaults to empty
    vol surface is 0.001 au isodensity

    Returns:
        nested dictionaries with keys 'Group', 'Atomic', 'BCP', 'VSCC'
    """
    atomList = list(
        _get_table(
            data,
            "q(A)              L(A)              K(A)          K_Scaled(A)      |Mu_Intra(A)|",
        )["Atom"]
    )
    _check_num_atoms(atomList, subAtoms)
    subatomLabels = [
        atomList[i]
        for i in range(0, len(atomList))
        if any(x == i + 1 for x in subAtoms)
    ]

    # if bcpIdx not empty, get bcp properties for each listed bcp

    atomicProps = get_atomic_props(data)  # get atomic dictionary
    bcpProperties = extract_requested_bcp_props(
        data, atomList, bcpId, subatomLabels, atomicProps
    )
    if groupProps:  # if you want group properties
        subDict = get_sub_props(
            atomicProps, subAtoms, atomList
        )  # get substituent properties
    else:
        subDict = "Group Properties not requested"

    if len(lapRhoCpAtoms) > 0:  # if we want laprhocps for at least one atom
        vsccProps = extract_requested_cc_props(
            lapRhoCpAtoms, sumFileNoExt, atomList, atomicProps
        )
    else:
        vsccProps = {"VSCC Props not requested"}
    # create output dictionary to return all requested properties, if a property not requested return a string stating that
    outDict = {
        "Group": subDict,
        "Atomic": atomicProps,
        "BCP": bcpProperties,
        "VSCC": vsccProps,
    }
    return outDict


def extract_requested_cc_props(lapRhoCpAtoms, sumFileNoExt, atomList, atomicProps):
    """Get VSCC dict for requested atoms"""
    vsccProps = {}
    for atom in lapRhoCpAtoms:  # for each atom requested, get lapRhoCps
        allCC = get_cc_props(sumFileNoExt, atomList[atom - 1])
        vsccProps.update(
            {atomList[atom - 1]: identify_vscc(allCC, atomicProps, atomList[atom - 1])}
        )
    return vsccProps


def extract_requested_bcp_props(data, atomList, bcpId, subatomLabels, atomicProps):
    """Get BCP dict for requested bcps"""
    bcpIdx = []
    for bcp in bcpId:
        atList = []
        for at in bcp:
            atList.append(at)
        bcpIdx.append(atList)
    bcpProperties = {}  # initialize empty dictionary
    if len(bcpIdx) > 0:
        for bcpPair in bcpIdx:
            bcpPair[0] -= 1  # lists start from 0, indices would be passed from molecule
            bcpPair[1] -= 1  # which would start at 1, so adjust index
            # add the bcp properties to the dictionary under the atom labels
            prop = get_bcp_properties(data, [atomList[i] for i in bcpPair])
            prop.update({"DI(R,G)": [get_sub_di(data, subatomLabels)]})
            for key in atomicProps:
                if key in (atomList[bcpPair[0]], atomList[bcpPair[1]]):
                    if "1" in key:
                        keynum = 1
                    elif "2" in key:
                        keynum = 2
                    prop.update(
                        {
                            f"r(BCP-{keynum})": [
                                _get_dist(atomicProps[key]["xyz"], prop["xyz"])
                            ]
                        }
                    )
            bcpProperties.update(
                {f"{atomList[bcpPair[0]]}-{atomList[bcpPair[1]]}": prop}
            )
    else:  #
        bcpProperties = "BCP Properties not requested"
    return bcpProperties


def get_selected_bcps(data, bcp_list):
    """Given a list of bcps, return dictionary of bcp properties."""
    bcpProperties = {}  # initialize empty dictionary
    for bcpPair in bcp_list:
        prop = get_bcp_properties(data, bcpPair)
        bcpProperties.update({f"{bcpPair[0]}-{bcpPair[1]}": prop})
    return bcpProperties


def find_connected(data, negXAtomLabel, originAtomLabel):
    """Given lines of sumfile, atom on -x label, and atom on origin Label, find atoms bonded to origin.
    Args:
        data: list[str]: lines of sumfile
        negXAtomLabel: str: eg. 'H2'
        originAtomLabel: str: eg. 'C1'
    Returns:
        List of List of BCPs connected to origin atom
        e.g.
            [['C1',H3'],['C1','H4'],['C1','H5']]
    """
    # find all atoms to which bonded
    bcpLines = []
    for line in data:
        split_line = line.split()
        if (
            "(3,-1)" in split_line
            and negXAtomLabel not in split_line
            and originAtomLabel in split_line
        ):
            bcpLines.append(line)
    bcpList = []
    for bcp in bcpLines:
        splitbcp = bcp.split()
        bcpList.append([splitbcp[4], splitbcp[5]])
    return bcpList


def find_all_connections(data):
    """find all pairs of atoms connected by a BCP"""
    bcpLines = []
    for line in data:
        split_line = line.split()
        if "(3,-1)" in split_line:
            bcpLines.append(line)
    bcpList = []
    for bcp in bcpLines:
        splitbcp = bcp.split()
        bcpList.append([splitbcp[4], splitbcp[5]])
    return bcpList


def sub_prop_frame(csvFile: str) -> dict:  # pylint:disable=too-many-locals
    """Given csv file, extract group properties for all files included and store properties
    Args:
        csvFile: string containing csv file (WITH extension)
        example csvFile structure:
        Substituent, subAtoms,    label1,   label2,...
        CH3          1 3 4 5      SubH_CH3  SubC6H5_CH3
            Substituent: string of substituent
            subAtoms: string of space separated substituent atoms eg '1 3 4'
            label1 - contains .sum file with no extension
            labeli could be e.g. a model chemistry, or substrate
            that is the column "SubH" would have sum files for the substituents attached to H
            SubC6H5 would have sum files for substituents attached to C6H5 etc

    Returns:
        dictionary of dicitonary of data frames containing group properties
        {'label1': {
            'Group': Pandas Data Frame
            'BCP': Pandas Data Frame containing BCP properties for bcp between atoms 1 and 2
        },
        'label2': same as label1 but for second label
        }
    Notes: Group frame has columns: Substituent, q,K, K_Scaled, Mu_Intra_X, Mu_Intra_Y, Mu_Intra_Z,
                Mu_Bond_X, Mu_Bond_Y,Mu_Bond_Z, Mu_X, Mu_Y, Mu_Z, Q_xx, Q_xy, Q_xz, Q_yy,Q_yz,Q_zz,Vol,
                |Mu_Intra|,|Mu_Bond|,|Mu|
            BCP frame has columns: Substituent, rho, delsqrho, lambda1, lambda2, lambda3,V,G,H,DI
    """
    csvFrame = pd.read_csv(csvFile)

    all_label_dict = {}

    ncolumns = csvFrame.shape[1]
    nrow = csvFrame.shape[0]
    subAtoms = []
    for sub in range(0, nrow):  # get subAtoms from csv in list format
        subAtomsString = csvFrame.loc[sub]["subAtoms"].split()
        subAtomsInt = [ast.literal_eval(i) for i in subAtomsString]
        subAtoms.append(subAtomsInt)
    for col in range(
        2, ncolumns
    ):  # 2 is the start of the label columns, so for all labels:  # create empty dictionary for this label
        count = 0  # For first iteration for each label will create the data frame
        for row_num, sub in enumerate(
            csvFrame["Substituent"]
        ):  # iterate over substituents
            sumFileName = csvFrame[csvFrame["Substituent"] == sub][
                csvFrame.columns[col]
            ].iloc[0]
            # get properties
            with open(sumFileName + ".sum", encoding="utf-8") as sumFile:
                data = sumFile.readlines()
            extracted_props = extract_sub_props(data, subAtoms[row_num], sumFileName)
            excludeKeys = ["xyz"]
            extracted_props["Group"].update({"Substituent": sub})
            if count == 0:
                count = 1  # don't come here again for this label (don't need to make data frame again)
                groupFrame = pd.DataFrame.from_dict(
                    extracted_props["Group"], orient="columns"
                )
                # for each bcp properties gotten for - currently only do this for one
                # should work generally soon
                for bnum, bcp in enumerate(
                    extracted_props["BCP"]
                ):  # currently only use for one bcp query
                    tempbcpDict = {
                        k: extracted_props["BCP"][bcp][k]
                        for k in set(list(extracted_props["BCP"][bcp].keys()))
                        - set(excludeKeys)
                    }
                    tempbcpDict.update({"x": [extracted_props["BCP"][bcp]["xyz"][0]]})
                    tempbcpDict.update({"y": [extracted_props["BCP"][bcp]["xyz"][1]]})
                    tempbcpDict.update({"z": [extracted_props["BCP"][bcp]["xyz"][2]]})
                    tempbcpDict.update({"Substituent": sub})
                    tempbcpDict.update({"BCP": bcp})
                    if bnum == 0:  # create data frame for first bcp in list
                        bcpFrame = pd.DataFrame.from_dict(tempbcpDict, orient="columns")
                    else:  # else add to data frame
                        bcpFrame = pd.concat(
                            [bcpFrame, pd.DataFrame(tempbcpDict)], ignore_index=True
                        )
            else:  # add to data frame after first iteration
                groupFrame = pd.concat(
                    [groupFrame, pd.DataFrame(extracted_props["Group"])],
                    ignore_index=True,
                )
                for bcp in extracted_props["BCP"]:
                    tempbcpDict = {
                        k: extracted_props["BCP"][bcp][k]
                        for k in set(list(extracted_props["BCP"][bcp].keys()))
                        - set(excludeKeys)
                    }
                    tempbcpDict.update({"x": [extracted_props["BCP"][bcp]["xyz"][0]]})
                    tempbcpDict.update({"y": [extracted_props["BCP"][bcp]["xyz"][1]]})
                    tempbcpDict.update({"z": [extracted_props["BCP"][bcp]["xyz"][2]]})
                    tempbcpDict.update({"Substituent": sub})
                    tempbcpDict.update({"BCP": bcp})
                    bcpFrame = pd.concat(
                        [bcpFrame, pd.DataFrame(tempbcpDict)], ignore_index=True
                    )
        all_label_dict.update(
            {csvFrame.columns[col]: {"Group": groupFrame, "BCP": bcpFrame}}
        )
    return all_label_dict


def get_xyz(data) -> dict:
    """Given sumfile, return dicitonary containing xyzcoordinates and atom labels

    Args:
        sumfile: string with sumfile to be used, without .sum extenstion

    Returns:
        Dictionary
            'xyz': dataframe of xyz coordinates of nuclei in sum file
            'Atoms': list of atom labels

    """

    xyzTable = _get_table(
        data,
        "Charge                X                  Y                  Z",
        endString="Some Atomic Properties:",
    )
    xyzTable["X"] = pd.to_numeric(xyzTable["X"], downcast="float")
    xyzTable["Y"] = pd.to_numeric(xyzTable["Y"], downcast="float")
    xyzTable["Z"] = pd.to_numeric(xyzTable["Z"], downcast="float")
    return {"xyz": xyzTable[["X", "Y", "Z"]].to_numpy(), "Atoms": xyzTable["Atom"]}