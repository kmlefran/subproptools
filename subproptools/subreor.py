import sys
import os
sys.path.append(sys.path[0].replace('subproptools','')+'/'+'referenceMaps')
import pandas as pd #data frames
import math #sqrt
import numpy as np
from subproptools import qtaimExtract as qt #sum file manipulation and property extraction

#define stat dicitnoary for scaling bcp properties
_DEFAULT_STAT_DICT = {'rho': {'mean':0.290686,'sd':0.077290},
                 'lambda1':{'mean':-0.725552,'sd':0.299756},
                 'lambda2':{'mean':-0.678830,'sd':0.291123},
                 'lambda3':{'mean':0.583261,'sd':0.449474},
                 'DelSqRho':{'mean':-0.821120,'sd':0.570553},
                 'Ellipticity':{'mean':0.077722,'sd':0.137890},
                 'V':{'mean':-0.475963,'sd':0.332327},
                 'G':{'mean':0.135342,'sd':0.184923},
                 'H':{'mean':-0.340622,'sd':0.176627},
                 'DI(R,G)':{'mean':1.081894,'sd':0.369324},}

#define reference maps for rotation
_SP3_CARBON_BCP_DICT = {
    'C1-H3': {'xyz': np.array([ 0.        ,  1.25015821, -0.00235558]),
  'rho': [0.29632696936],
  'lambda1': [-0.85302912761],
  'lambda2': [-0.81303659844],
  'lambda3': [0.47413487677],
  'DelSqRho': [-1.1919308493],
  'Ellipticity': [0.049189088474],
  'V': [-0.37184352147],
  'G': [0.036930404577],
  'H': [-0.33491311689299996]},
 'C1-F8': {'xyz': np.array([ 0.        , -0.35624989, -0.74223654]),
  'rho': [0.23829406815],
  'lambda1': [-0.4245753785],
  'lambda2': [-0.41016717535],
  'lambda3': [0.89122326411],
  'DelSqRho': [0.056480710264],
  'Ellipticity': [0.035127635803],
  'V': [-0.68209769455],
  'G': [0.34810893606],
  'H': [-0.33398875849000004]},
 'C1-C4': {'xyz': np.array([ 0.        , -0.7598954 ,  1.12297281]),
  'rho': [0.26194645384],
  'lambda1': [-0.53899298565],
  'lambda2': [-0.51656172708],
  'lambda3': [0.3428304022],
  'DelSqRho': [-0.71272431054],
  'Ellipticity': [0.043424159007],
  'V': [-0.29901033382],
  'G': [0.060414628094],
  'H': [-0.23859570572599997]}
}

_SP2_CARBON_BCP_DICT = {'C1-H3': {'xyz': np.array([0.00000000e+00, 1.16430959e+00, 3.38675717e-17]),
  'rho': [0.29183770628],
  'lambda1': [-0.80520264462],
  'lambda2': [-0.79941295281],
  'lambda3': [0.44085193209],
  'DelSqRho': [-1.1637636653],
  'Ellipticity': [0.0072424293269],
  'V': [-0.37436126977],
  'G': [0.041710176716],
  'H': [-0.332651093054]},
 'C1-C4': {'xyz': np.array([ 0.00000000e+00, -1.06357498e+00, -5.56714996e-16]),
  'rho': [0.36288408448],
  'lambda1': [-0.83468307102],
  'lambda2': [-0.62545003471],
  'lambda3': [0.15660468687],
  'DelSqRho': [-1.3035284189],
  'Ellipticity': [0.3345319765],
  'V': [-0.62700069759],
  'G': [0.15055929644],
  'H': [-0.47644140115]}}

_SP3_BORON_BCP_DICT = {'B1-H3': {'xyz': np.array([ 0.        ,  0.92742229, -0.18477357]),
  'rho': [0.12300243683],
  'lambda1': [-0.17214219218],
  'lambda2': [-0.14265466718],
  'lambda3': [0.42509227636],
  'DelSqRho': [0.110295417],
  'Ellipticity': [0.20670564506],
  'V': [-0.23579025344],
  'G': [0.13168205384],
  'H': [-0.1041081996]},
 'B1-C5': {'xyz': np.array([ 0.        , -0.22289873, -0.93254826]),
  'rho': [0.12047147227],
  'lambda1': [-0.14773895803],
  'lambda2': [-0.11298066129],
  'lambda3': [0.33173067015],
  'DelSqRho': [0.071011050836],
  'Ellipticity': [0.30764819702],
  'V': [-0.22525108014],
  'G': [0.12150192142],
  'H': [-0.10374915872]},
 'B1-H4': {'xyz': np.array([ 0.        , -0.43919783,  0.70581714]),
  'rho': [0.18145419362],
  'lambda1': [-0.38722702492],
  'lambda2': [-0.37404958811],
  'lambda3': [0.45979611296],
  'DelSqRho': [-0.30148050007],
  'Ellipticity': [0.035229117263],
  'V': [-0.31882116911],
  'G': [0.12172552204],
  'H': [-0.19709564707000002]}}

_SP2_BORON_BCP_DICT = {'B1-H3': {'xyz': np.array([ 0.00000000e+00,  8.11870811e-01, -4.43740716e-16]),
  'rho': [0.19110739213],
  'lambda1': [-0.46141987776],
  'lambda2': [-0.40921992633],
  'lambda3': [0.44346099627],
  'DelSqRho': [-0.42717880782],
  'Ellipticity': [0.12755965211],
  'V': [-0.32443350486],
  'G': [0.10881940145],
  'H': [-0.21561410340999998]},
 'B1-F4': {'xyz': np.array([ 0.00000000e+00, -7.37862359e-01,  1.93163647e-16]),
  'rho': [0.20662382934],
  'lambda1': [-0.67411870986],
  'lambda2': [-0.56240110796],
  'lambda3': [2.5243450739],
  'DelSqRho': [1.2878252561],
  'Ellipticity': [0.19864399328],
  'V': [-0.6323324474],
  'G': [0.47714438071],
  'H': [-0.15518806669000007]}}

_SP3_NITROGEN_BCP_DICT = {'N1-H3': {'xyz': np.array([0.        , 1.39422451, 0.00543887]),
  'rho': [0.33833704323],
  'lambda1': [-1.2918559059],
  'lambda2': [-1.2716835],
  'lambda3': [0.85793011206],
  'DelSqRho': [-1.7056092939],
  'Ellipticity': [0.015862756642],
  'V': [-0.52200823266],
  'G': [0.047802954594],
  'H': [-0.474205278066]},
 'N1-O4': {'xyz': np.array([ 0.        , -0.63122426, -1.0353918 ]),
  'rho': [0.35088013303],
  'lambda1': [-0.78650875363],
  'lambda2': [-0.78506939426],
  'lambda3': [1.3832612581],
  'DelSqRho': [-0.18831688977],
  'Ellipticity': [0.0018334167423],
  'V': [-0.53352667926],
  'G': [0.24322372841],
  'H': [-0.29030295084999996]},
 'N1-C5': {'xyz': np.array([ 0.        , -0.65245651,  1.47278505]),
  'rho': [0.24674849376],
  'lambda1': [-0.48766358511],
  'lambda2': [-0.46565621172],
  'lambda3': [0.32464807311],
  'DelSqRho': [-0.62867172372],
  'Ellipticity': [0.047260989639],
  'V': [-0.35743098181],
  'G': [0.10013152544],
  'H': [-0.25729945637]}}

_SP3_PHOSPHOROUS_BCP_DICT = {'P1-H3': {'xyz': np.array([0.        , 1.26385682, 0.00494592]),
  'rho': [0.17523148631],
  'lambda1': [-0.30194810012],
  'lambda2': [-0.28264362027],
  'lambda3': [0.44213586293],
  'DelSqRho': [-0.14245585747],
  'Ellipticity': [0.068299719025],
  'V': [-0.33080503772],
  'G': [0.14759553668],
  'H': [-0.18320950103999997]},
 'P1-O4': {'xyz': np.array([ 0.        , -0.58969576, -0.83336938]),
  'rho': [0.2366336491],
  'lambda1': [-0.39519084069],
  'lambda2': [-0.39303248004],
  'lambda3': [2.1516401064],
  'DelSqRho': [1.3634167857],
  'Ellipticity': [0.005491557963],
  'V': [-0.74971032231],
  'G': [0.54528225936],
  'H': [-0.20442806295000004]},
 'P1-C5': {'xyz': np.array([ 0.        , -0.35462054,  1.23291434]),
  'rho': [0.17350204252],
  'lambda1': [-0.26162785795],
  'lambda2': [-0.24116223985],
  'lambda3': [0.29034324591],
  'DelSqRho': [-0.21244685189],
  'Ellipticity': [0.08486244827],
  'V': [-0.31020609825],
  'G': [0.12854719264],
  'H': [-0.18165890561]}}

_SP3_SI_BCP_DICT = {'Si1-H3': {'xyz': np.array([0.        , 1.27973328, 0.00517033]),
  'rho': [0.12414909335],
  'lambda1': [-0.20581129439],
  'lambda2': [-0.19632482872],
  'lambda3': [0.58847810247],
  'DelSqRho': [0.18634197935],
  'Ellipticity': [0.048320254397],
  'V': [-0.21031525709],
  'G': [0.12845037596],
  'H': [-0.08186488113000001]},
 'Si1-F8': {'xyz': np.array([ 0.        , -0.53009648, -1.08479626]),
  'rho': [0.13244199788],
  'lambda1': [-0.27302642228],
  'lambda2': [-0.27272649439],
  'lambda3': [1.5898643458],
  'DelSqRho': [1.0441114291],
  'Ellipticity': [0.0010997387276],
  'V': [-0.33781650703],
  'G': [0.29942218215],
  'H': [-0.03839432488]},
 'Si1-C4': {'xyz': np.array([ 0.        , -0.71188794,  1.04126176]),
  'rho': [0.12665648789],
  'lambda1': [-0.1843126043],
  'lambda2': [-0.1757363122],
  'lambda3': [0.56650837256],
  'DelSqRho': [0.20645945606],
  'Ellipticity': [0.048802048906],
  'V': [-0.21962871024],
  'G': [0.13562178713],
  'H': [-0.08400692310999999]}}

_SP2_SI_BCP_DICT = {'Si1-H3': {'xyz': np.array([ 0.        ,  1.23757956, -0.00766368]),
  'rho': [0.12096117641],
  'lambda1': [-0.1905393755],
  'lambda2': [-0.17657152668],
  'lambda3': [0.57387809704],
  'DelSqRho': [0.20676719486],
  'Ellipticity': [0.079105895977],
  'V': [-0.20375573354],
  'G': [0.12772376613],
  'H': [-0.07603196741000001]},
 'Si1-NNA7': {'xyz': np.array([ 0.        , -1.26148618,  0.24490124]),
  'rho': [0.10836999179],
  'lambda1': [-0.11714251117],
  'lambda2': [-0.057140471201],
  'lambda3': [0.092875642396],
  'DelSqRho': [-0.081407339975],
  'Ellipticity': [1.0500795444],
  'V': [-0.13782730849],
  'G': [0.058737736748],
  'H': [-0.07908957174199999]}}

_REFERENCE_MAP = {
    'C':{'sp3':_SP3_CARBON_BCP_DICT,'sp2':_SP2_CARBON_BCP_DICT},
    'B':{'sp3':_SP3_BORON_BCP_DICT,'sp2':_SP2_BORON_BCP_DICT},
    'N':{'sp3':_SP3_NITROGEN_BCP_DICT},
    'P':{'sp3':_SP3_PHOSPHOROUS_BCP_DICT},
    'Si':{'sp3':_SP3_SI_BCP_DICT,'sp2':_SP2_SI_BCP_DICT},
}

def _get_bcp_reference(originAtom,numBonds):
    """Takes atom and number of bonds and chooses the right map to use."""
    if originAtom =='C' and numBonds==4: #sp3 carbon
        retDict = _REFERENCE_MAP['C']['sp3']
    elif originAtom =='C' and numBonds==3: #sp2 carbon
        retDict = _REFERENCE_MAP['C']['sp2']
        #note linear shouldn't get to this point
    elif originAtom =='B' and numBonds==3: #planar boron
        # print('planar boron')
        retDict = _REFERENCE_MAP['B']['sp2']
    elif originAtom =='B' and numBonds==4: #sp3 boron
        # print('sp3 boron')
        retDict = _REFERENCE_MAP['B']['sp3']
    elif originAtom == 'N' and numBonds==4:
        # print('ammonium')
        retDict = _REFERENCE_MAP['N']['sp3']
    # elif originAtom =='Al' and numBonds==3: #planar boron
    #     print('planar aluminum')
    # elif originAtom =='Al' and numBonds==4: #sp3 boron
    #     print('sp3 aluminum') 
    elif originAtom =='Si' and numBonds==4: #sp3 carbon
        # print('sp3 silicon')
        retDict = _REFERENCE_MAP['Si']['sp3']
    elif originAtom =='Si' and numBonds==3: #sp2 carbon
        # print('sp2 silicon')
        retDict = _REFERENCE_MAP['Si']['sp2']
    elif originAtom == 'P' and numBonds==4:
        # print('phosphonium')
        retDict = _REFERENCE_MAP['P']['sp3']
    return retDict 

# def _default_stat_dict():
#     """defines mean and sd for scaling properties if no other available
    
#     Scaling based on mean & sd of all bcps of files used in transferability paper(not just Substrate-Substituent)
#     See:notebooks/bcpStatExtraction.ipynb
#     Returns:
#       a dictionary for each property, with names matching the names returned in a get_bcp_properties dictionary
#     """
#     stat_dict = {'rho': {'mean':0.290686,'sd':0.077290},
#                  'lambda1':{'mean':-0.725552,'sd':0.299756},
#                  'lambda2':{'mean':-0.678830,'sd':0.291123},
#                  'lambda3':{'mean':0.583261,'sd':0.449474},
#                  'DelSqRho':{'mean':-0.821120,'sd':0.570553},
#                  'Ellipticity':{'mean':0.077722,'sd':0.137890},
#                  'V':{'mean':-0.475963,'sd':0.332327},
#                  'G':{'mean':0.135342,'sd':0.184923},
#                  'H':{'mean':-0.340622,'sd':0.176627},
#                  'DI(R,G)':{'mean':1.081894,'sd':0.369324},}
#     return stat_dict

def _find_bcp_match(data,originAtomXYZ,negXAtomLabel, originAtomLabel):
    """
    Finds the atoms connected to the origin atom, and arranges the BCPs in a clockwise sense
    (assuming -x atom going into page)

    Args:
        data: the lines of a .sum file stored as a list
        originAtomXYZ - np.array of xyz coordinates of atom that was to be set to origin
        negXAtomLabel - the atom connected to origin that is positioned along -x axis. e.g. "H2"
        originAtomLabel - the origin atom of substituent bonded to substrate. Format e.g. "C1"

    Returns:
        Dictionary of BCP properties for atoms (non-neg-x atom) that are bonded to the origin atom.
        These are ordered in a clockwise rotational sense        
    """

    #find bcps connected to origin atom that are not the -x axis atom
    originBCPs = qt._find_connected(data,negXAtomLabel,originAtomLabel)
    # print(originBCPs)
    bcpPropDict = {}
    #get the bcp properties
    for bcp in originBCPs:
        #bcpBlock = qt.lock()
        bcpPropDict.update({bcp[0]+'-'+bcp[1]: qt.get_bcp_properties(data,atPair=bcp)})
    # print(bcpPropDict)
    if len(bcpPropDict) == 2:
        clockwiseKeys = []
    else:
        clockwiseKeys = _find_clockwise_rot(bcpPropDict,originAtomLabel,originAtomXYZ)
    #at this point have bcpDictionary ordered from 1st to last with clockwise bcp
    return clockwiseKeys #this is a dictionary of bcps

# def _angle_btw_bcp(xyzA,xyzB,atomXYZ=np.array([0.,0.,0.])):
#     """find angle between two bcp orientation vectors defined by x,y,z np.array, after flattening to yz plane."""
#     xyzA[0], xyzA[1], xyzA[2] = [xyzA[0]-atomXYZ[0], xyzA[1]-atomXYZ[1], xyzA[2]-atomXYZ[2]]
#     xyzB[0], xyzB[1], xyzB[2] = [xyzB[0]-atomXYZ[0], xyzB[1]-atomXYZ[1], xyzB[2]-atomXYZ[2]]
#     angle= math.acos((xyzA[0]*xyzB[0]+xyzA[1]*xyzB[1]+xyzA[2]*xyzB[2])/(math.sqrt(xyzA[0]**2+xyzA[1]**2+xyzA[2]**2)*math.sqrt(xyzB[0]**2+xyzB[1]**2+xyzB[2]**2)))
#     return angle

def _find_clockwise_rot(bcpPropDict,originAtomLabel,originAtomXYZ=np.array([0.0,0.0,0.0])):
    """given dictionary of bcp properties, find which ones are in a clockwise rotation"""
    #return list of dictionary keys ordered for clockwise rotation
    # print(bcpPropDict)
    # print(originAtomLabel)
    crossDict = {}
    for key1 in bcpPropDict:
        print(key1)
        for key2 in bcpPropDict:
            print(key2)
            if key1 != key2:
                xyz1 = bcpPropDict[key1]['xyz']
                xyz1[0] = 0.#project to yz plane
                xyz2=bcpPropDict[key2]['xyz']
                xyz2[0]=0.
                cross=np.cross(np.subtract(bcpPropDict[key1]['xyz'],originAtomXYZ),np.subtract(bcpPropDict[key2]['xyz']),originAtomXYZ)[0]
                if cross < 0:
                    isClockwise = True
                else:
                    isClockwise = False
                if isClockwise:
                    tempDict={'Is Clockwise':isClockwise,'cross':cross}
                    bondAt1 = key1.replace(originAtomLabel+'-','')
                    bondAt1 = bondAt1.replace('-'+originAtomLabel,'')
                    bondAt2 = key2.replace(originAtomLabel+'-','')
                    bondAt2 = bondAt2.replace('-'+originAtomLabel,'')
                    crossDict.update({bondAt1 + '-To-' + bondAt2:tempDict})
    orderDict= {}                
    for cw in crossDict:
        string = cw.split('-')
        start = string[0]
        end = string[2]
        orderDict.update({cw:{'Start':start,'End':end}})
#         keyDict={} #only one loop with key
#         xyz = bcpPropDict[key]['xyz']-originAtomXYZ
#         keyDict.update({'AngleToY': _angle_btw_bcp(xyz,np.array([0.0,1.0,0.0]))})
#         keyDict.update({'AngleToZ': _angle_btw_bcp(xyz,np.array([0.0,0.0,1.0]))})
#         angleDict.update({key:keyDict})
    keysList = list(orderDict.keys())
    if len(keysList) == 3:
        if orderDict[keysList[0]]['End'] != orderDict[keysList[1]]['Start'] and orderDict[keysList[0]]['End'] == orderDict[keysList[2]]['Start']:
            reordered_dict = {k : orderDict[k] for k in [keysList[0],keysList[2],keysList[1]]}
        else:
            reordered_dict = orderDict
        ordered_bcp_props = {}
        for order in reordered_dict:
            for bcp in bcpPropDict:
                print(bcp)
                if reordered_dict[order]['Start'] in bcp:
                    ordered_bcp_props.update({bcp:bcpPropDict[bcp]})
                    continue
    else:
        ordered_bcp_props = bcpPropDict
    return ordered_bcp_props

def _set_origin(xyzArray,originAtom):
    """shifts all points in xyz geometry by originAtom coordinates, returns xyz geometry"""
    org = xyzArray[originAtom-1,]
    return xyzArray - org

def _get_lengths(xyzArray):
    """given xyz geometry returns np.array of magnitudes of squared distances from the origin"""
    lengths = np.array([])
    for atom in xyzArray:
        length=0
        for coord in atom:
            length += coord**2
        lengths = np.append(lengths,length)
    return lengths

def _zero_y_for_negx(t_xyz,negXAtom):
    """perform first rotation in setting -x axis to zero the y value"""
    if t_xyz[0,negXAtom-1] == 0 and t_xyz[1,negXAtom-1] == 0: #on xy
        print('hi')
        G= np.array([[0.,0.,1.],[0.,1.,0.],[-1.,0.,0.]])
    elif t_xyz[1,negXAtom-1] == 0 and t_xyz[2,negXAtom-1] == 0: #already on x axis atom 2 y and z=0
        print('here')
        G = np.array([[1.,0.,0.],[0.,1.,0.],[0.,0.,1.]])
    elif t_xyz[0,negXAtom-1] == 0 and t_xyz[2,negXAtom-1] == 0:
        print('ho there')
        G = np.array([[0.,1.,0.],[-1.,0.,0.],[0.,0.,1.]])
    else:
        print('no way')
        d = t_xyz[0,negXAtom-1]/math.sqrt(t_xyz[0,negXAtom-1]**2 + t_xyz[1,negXAtom-1]**2)
        s = t_xyz[1,negXAtom-1]/math.sqrt(t_xyz[0,negXAtom-1]**2 + t_xyz[1,negXAtom-1]**2)
        G = np.array([[d,s,0.],[-s,d,0.],[0.,0.,1.]])
    return np.matmul(G, t_xyz)   

def _zero_z_for_negx(t_xyz,negXAtom):
    """perform second rotation in setting -x axis to zero the z value"""
    #perform after y
    d = t_xyz[0,negXAtom-1]/math.sqrt(t_xyz[0,negXAtom-1]**2 + t_xyz[2,negXAtom-1]**2)
    s = t_xyz[2,negXAtom-1]/math.sqrt(t_xyz[0,negXAtom-1]**2 + t_xyz[2,negXAtom-1]**2)
    G = np.array([[d,0,s],[0,1,0],[-s,0,d]])
    return np.matmul(G,t_xyz)

def _set_xaxis(xyzArray,negXAtom):
    """Given xyz geometry with atom at origin, return xyz geometry with negXAtom on -x."""
    t_xyz = xyzArray.T
    
    #define initial xyz vector lengths. Should be unchanged after rotation
    tol = 0.0001 #tolerance for change
    initial_lengths= _get_lengths(xyzArray)

    #Set Givens matrix for first rotation to zero [2,2]
#     if t_xyz[0,1] == 0 and t_xyz[1,1] == 0: #on xy
#         G= np.array([[0,0,1],[0,1,0],[-1,0,0]])
#     elif t_xyz[1,1] == 0 and t_xyz[2,1] == 0: #already on x axis atom 2 y and z=0
#         G = np.array([[1,0,0],[0,0,1],[0,0,1]])
#     elif t_xyz[0,1] == 0 and t_xyz[2,1] == 0:
#         G = np.array([[0,1,0],[-1,0,0],[0,0,1]])
#     else:
#         d = t_xyz[0,1]/math.sqrt(t_xyz[0,1]**2 + t_xyz[1,1]**2)
#         s = t_xyz[1,1]/math.sqrt(t_xyz[0,1]**2 + t_xyz[1,1]**2)
#         G = np.array([[d,s,0],[-s,d,0],[0,0,1]])
    t_rot1 = _zero_y_for_negx(t_xyz,negXAtom)
    rot1_lengths = _get_lengths(t_rot1.T)
    if np.any((rot1_lengths - initial_lengths)>tol):
        print("Geometry perturbation exceeded")
        raise AssertionError('Geometry change after rotation exceeded tolerance')
    t_rot2 = _zero_z_for_negx(t_rot1,negXAtom)
    rot2_lengths = _get_lengths(t_rot2.T)
    
    if np.any((rot2_lengths - initial_lengths)>tol):
        print("Geometry perturbation exceeded")
        raise AssertionError('Geometry change after rotation exceeded tolerance')
    if t_rot2[0,negXAtom-1] > 0: #if negxatom is along +x, rotate 180
        G = np.array([[-1,0,0],[0,1,0],[0,0,-1]])
        t_rot_final = np.matmul(G,t_rot2)
        return t_rot_final.T
    else:
        return t_rot2.T
    
def _align_dicts(testDict,refDict,statDict):
    """
    Arguments:
        testDict and refDict:For two dictionaries that are ordered in same rotational sense
        statDict: scaling to use
    Returns:
        np.array of xyz point to use
    """
    testKeysList = list(testDict.keys())
    refDictKeysList = list(refDict.keys())
    #BCP distances for first match
    dif00 = _get_popelier_dif(testDict[testKeysList[0]],refDict[refDictKeysList[0]],statDict)
    dif11 = _get_popelier_dif(testDict[testKeysList[1]],refDict[refDictKeysList[1]],statDict)
    dif22 = _get_popelier_dif(testDict[testKeysList[2]],refDict[refDictKeysList[2]],statDict)
    #BCP distances for second match
    dif10 = _get_popelier_dif(testDict[testKeysList[1]],refDict[refDictKeysList[0]],statDict)
    dif21 = _get_popelier_dif(testDict[testKeysList[2]],refDict[refDictKeysList[1]],statDict)
    dif02 = _get_popelier_dif(testDict[testKeysList[0]],refDict[refDictKeysList[2]],statDict)
    #BCP distances for third match
    dif20 = _get_popelier_dif(testDict[testKeysList[2]],refDict[refDictKeysList[0]],statDict)
    dif01 = _get_popelier_dif(testDict[testKeysList[0]],refDict[refDictKeysList[1]],statDict)
    dif12 = _get_popelier_dif(testDict[testKeysList[1]],refDict[refDictKeysList[2]],statDict)
    #Total score list and index of closest space
    matchScores = [dif00+dif11+dif22,dif10+dif21+dif02,dif20+dif01+dif12]
    minInd = matchScores.index(min(matchScores))
    
    #identify the refDict BCP that is on the +y axis
    if refDict[refDictKeysList[0]]['xyz'][1] > 0 and abs(refDict[refDictKeysList[0]]['xyz'][2]) < 0.01:
        refPosY=0
    elif refDict[refDictKeysList[1]]['xyz'][1] > 0 and abs(refDict[refDictKeysList[1]]['xyz'][2]) < 0.01:
        refPosY=1
    elif refDict[refDictKeysList[2]]['xyz'][1] > 0 and abs(refDict[refDictKeysList[2]]['xyz'][2]) < 0.01: 
        refPosY=2
    #for the best match, set the posYPoint to the one that mapped to refDict +y point    
    if minInd==0:
        if refPosY==0:
            posYPoint = testDict[testKeysList[0]]['xyz']
        elif refPosY == 1:
            posYPoint = testDict[testKeysList[1]]['xyz']
        elif refPosY == 2:
            posYPoint = testDict[testKeysList[2]]['xyz']
    elif minInd==1:
        if refPosY==0:
            posYPoint = testDict[testKeysList[1]]['xyz']
        elif refPosY == 1:
            posYPoint = testDict[testKeysList[2]]['xyz']
        elif refPosY == 2:
            posYPoint = testDict[testKeysList[0]]['xyz']
    elif minInd==2:
        if refPosY==0:
            posYPoint = testDict[testKeysList[2]]['xyz']
        elif refPosY == 1:
            posYPoint = testDict[testKeysList[0]]['xyz']
        elif refPosY == 2:
            posYPoint = testDict[testKeysList[1]]['xyz']
    return posYPoint    

def _get_posy_point(sumFileNoExt,atomDict,attachedAtom,negXAtomLabel,default_stats=True):
    """returns point to put on +y axis matching definition in rotate_substituent."""
    ccProps = qt.get_cc_props(sumFileNoExt,attachedAtom)    
    if len(ccProps) > 0:
        vscc = qt.identify_vscc(ccProps,atomDict)
    else:
        vscc = {} 
    if len(vscc) == 1:
        #reorient setting vscc to +y
        vkeys = list(vscc.keys())
        posYPoint = vscc[vkeys[0]]['xyz']
    elif len(vscc) == 2 and "N" not in attachedAtom:
        vkeys = list(vscc.keys())
        
        posYPoint = (vscc[vkeys[0]]['xyz'] + vscc[vkeys[1]]['xyz'])/2
        # print(posYPoint)
        #reorient to average of vscc points for +y
    else:
        #refDict = _get_reference_map()
        sumFile = open(sumFileNoExt+".sum",'r')
        data = sumFile.readlines()
        sumFile.close()
        #bcpsToMatch is bcp dictionary, ordered for clockwise rot
        #data,originAtomXYZ,negXAtomLabel,originAtomLabel
        bcpsToMatch = _find_bcp_match(data,atomDict[attachedAtom]['xyz'],negXAtomLabel,attachedAtom)
        numBonds=len(bcpsToMatch)+1
        #on the assumption that if an atom has two bonds (_find_bap_match returns None), 
        # and does not have a lone pair, it is linear, so we do not do another rotation 
        # and posYPoint is None
        if len(bcpsToMatch) > 0:
            posYPoint = []
        else:
            atType = ''.join([i for i in attachedAtom if not i.isdigit()])
            matchDict = _get_bcp_reference(atType,numBonds)
            if default_stats:
                statDict = _DEFAULT_STAT_DICT
            posYPoint = _align_dicts(bcpsToMatch,matchDict,statDict)
        #reorient to another point
        #posy point will be the point that would lie along the y-axis in reference in maximal match case
        # print('not done yet')
    return posYPoint    

def _set_yaxis(xyzArray,posYArray):
    """rotate a geom positioned on -x axis so posYArray will lie on +y"""
    theta = math.atan2(posYArray[2],posYArray[1])
    # print(theta)
    c = math.cos(theta)
    s = math.sin(theta)
    G = np.array([[1,0,0],[0,c,s],[0,-s,c]])
    rot1 = np.matmul(G,xyzArray.T)
    rot1vec = np.matmul(G,posYArray)
    
    if rot1vec[1] < 0:
        G2 = np.array([[1,0,0],[0,-1,0],[0,0,-1]])
        final_geom = np.matmul(G2,rot1)
    else:
        final_geom = rot1
    return final_geom.T    


def _get_popelier_dif(bcpDictA,bcpDictB,statDict):
    """compute distance in BCP space between dictA and B scaling using statDict"""
    distancesq = 0.
    for prop in bcpDictA:
        if prop != 'xyz':
            scaledA = (bcpDictA[prop][0]-statDict[prop]['mean'])/statDict[prop]['sd']
            scaledB = (bcpDictB[prop][0]-statDict[prop]['mean'])/statDict[prop]['sd']
            distancesq += (scaledA - scaledB)**2
    return math.sqrt(distancesq)    


#commented out - have included reference bcps at data at top of file rather tahn calculating each time
# def _get_ref_bcps(sumfilenoext,atPairList,originAtom,originAtomXYZ=np.array([0.,0.,0.])):
#     """given reference sumfile, extract bcp properties for needed bcps"""
#     sumFile = open(sumfilenoext+".sum","r") #open file, read lines, close file
#     data = sumFile.readlines()
#     sumFile.close()
#     bcpDict = {}
#     for bcp in atPairList:
#         #block = qt.get_bcp_block(data,bcp)
#         bcpDict.update({'{at1}-{at2}'.format(at1=bcp[0],at2=bcp[1]):qt.get_bcp_properties(data,bcp)})
#     clockbcpDict = _find_clockwise_rot(bcpDict,originAtom,originAtomXYZ)    
#     return clockbcpDict  


# def _get_reference_map():
#     """creates dictionary with reference species to map coordinate system to"""
#     pathToRefs=sys.path[0].replace('subproptools','')+'/'+'referenceMaps'+'/'
#     print(pathToRefs)
#     refDict={}
#     carbonDict = {}
#     carbonDict.update({'sp3':_get_ref_bcps(pathToRefs + 'SubH_CFHCH3_reor_wb97xd_aug-cc-pvtz',[['C1','H3'],['C1','C4'],['C1','F8']],originAtom='C1')})
#     carbonDict.update({'sp2':_get_ref_bcps(pathToRefs +'SubH_CHCH2_wb97xd_aug-cc-pvtz_reor',[['C1','H3'],['C1','H2']],originAtom='C1')})
#     #carbonDict.update({'sp2':_get_ref_bcps('SubH_CHCH2-ReorJuly2-B3LYP-def2-TZVPPD-Field',[['C1','H3'],['C1','C4']],'C1')})
#     refDict.update({'C':carbonDict})
#     boronDict = {}
#     boronDict.update({'sp3':_get_ref_bcps(pathToRefs +'SubH_BHCH3BH2_wb97xd_aug-cc-pvtz_reor',[['B1','H3'],['B1','H4'],['B1','C5']],originAtom='B1')})
#     boronDict.update({'sp2':_get_ref_bcps(pathToRefs +'SubH_BHF_wb97xd_aug-cc-pvtz_reor',[['B1','H3'],['B1','F4']],originAtom='B1')})
#     refDict.update({'B':boronDict})
#     nitrogenDict = {}
#     nitrogenDict.update({'sp3':_get_ref_bcps(pathToRefs +'SubH_NHOCH3_wb97xd_aug-cc-pvtz_reor',[['N1','H3'],['N1','O4'],['N1','C5']],originAtom='N1')})
#     refDict.update({'N':nitrogenDict})
#     phosphorousDict = {}
#     phosphorousDict.update({'sp3':_get_ref_bcps(pathToRefs +'SubH_POCH3H_wb97xd_aug-cc-pvtz_reor',[['P1','H3'],['P1','O4'],['P1','C5']],originAtom='P1')})
#     refDict.update({'P':phosphorousDict})
#     siliconDict = {}
#     siliconDict.update({'sp3':_get_ref_bcps(pathToRefs +'SubH_SiFHCH3_wb97xd_augccpvtz_reor',[['Si1','H3'],['Si1','C4'],['Si1','F8']],originAtom='Si1')})
#     siliconDict.update({'sp2':_get_ref_bcps(pathToRefs +'SubH_Si2H3_wb97xd_augccpvtz_reor',[['Si1','H3'],['Si1','NNA7']],originAtom='Si1')})
#     refDict.update({'Si':siliconDict})
#     #update dictionary with further references
#     return refDict


def rotate_substituent(sumFileNoExt,originAtom,negXAtom,posYAtom=0):
    """
    Rotates a substituent to the defined coordinate system.

    Coordinate system defined as: 
        originAtom at (0,0,0)
        negXAtom at (-x,0,0)
        Atom on +y defined as:
            * average of lone pairs for 2 lone pairs on originAtom
            * Position of lone pair for 1 lone pair on originAtom
            * For no lone pairs: map BCPs onto reference for the atom type
            * Minimum distance in BCP space defined the atom to put on +y

    Args:
        sumFileNoExt (string): name of a sum file, without the .sum extension
        originAtom (int): the integer number of the atom to place at the origin
        negXAtom (int): the integer number of the atom to place along the -x axis
        posYAtom (int): override for above defined +y point, set to posYAtom instead

    Returns:
        pandas data frame of output geometry (columns Atom, x, y, z)

    Examples:
        >>> rotate_substituent('SubCH3_CFH2_anti2146_reor',1,1)
        Atom    x   y   z
        C1      0.  0.  0.
        H2 -{float} 0.  0.
        .(remaining geometry)
        .
        .
    """
    #read sum file
    sumFile = open(sumFileNoExt+".sum","r")
    data = sumFile.readlines()
    sumFile.close()
    atomDict = qt.get_atomic_props(data) #(needed for VSCC identification)

    molecule_xyz = qt.get_xyz(sumFileNoExt)
    #Labels format A1 etc
    negXAtomLabel = molecule_xyz['Atoms'][negXAtom-1]
    attachedAtom = molecule_xyz['Atoms'][originAtom-1]
    #perform reorientation
    molecule_orig = _set_origin(molecule_xyz['xyz'],originAtom)
    molecule_xaxis = _set_xaxis(molecule_orig,negXAtom)
    if posYAtom:
        posYPoint = molecule_xaxis[posYAtom-1] 
    else:    
        posYPoint = _get_posy_point(sumFileNoExt,atomDict,attachedAtom,negXAtomLabel)
    if len(posYPoint) > 0:    
        final_orientation = _set_yaxis(molecule_xaxis,posYPoint)
    else:
        final_orientation = molecule_xaxis    
    #Generate output
    outFrame = pd.DataFrame(final_orientation*0.529177,columns = ['x','y','z'])
    outFrame['Atom'] = molecule_xyz['Atoms']
    outFrame = outFrame[['Atom','x','y','z']]
    # with open(sumFileNoExt+'_reor.txt','w') as txt_file:
    #     of_string = outFrame.to_string(header=False,index=False)
    #     txt_file.write(of_string)
    return outFrame

def output_to_gjf(old_file_name,reor_geom,esm='wB97XD',basis_set="aug-cc-pvtz",add_label='',
                 n_procs=4,mem='3200MB',charge=0,multiplicity=1,wfx=True):
    """Given a rotated molecule, writes new geometry to single point Gaussian calculation
    
    Args:
        old_file_name - the file name of the sum file(no extension) before reorientation
        reor_geom - the data frame output of rotate_substituent
        esm - whatever electronic structure method is to be used in single point(HF/MP2/functional)
        basis_set - basis set to be used
        add_label - any extra label for file name, default empty
        n_procs=4 - numper of processors for remote clusters, set to 0 if not desired
        mem='3200MB' -amount of memory for remote clusters, set to 0 if not desired
        charge - charge of molecule
        multiplicity - multiplicity of molecule
        wfx - whether or not to write wfx, default True

    Returns:
        no return, but creates new gjf file old_file_name_reor_add_label.gjf
        File looks like:
        %chk=new_file_name.chk
        %nprocs=n_procs
        %mem=mem
        #p esm/basis_set output=wfx nosymmetry

        single point on old_file reoriented by subproptools

        charge multiplicity
        (xyz geom)

        new_file_name.wfx
        (blank lines)
    """
    new_file_name =old_file_name+'_reor'+add_label
    #delete file if already exists
    if os.path.exists(new_file_name+'.gjf'):
        # print('deleting')
        os.remove(new_file_name+'.gjf')
    with open(new_file_name+'.gjf', 'a') as f:
        f.write("%chk={chk}\n".format(chk=new_file_name+'.chk'))
        if n_procs:
            f.write('%nprocs={nprocs}\n'.format(nprocs=n_procs))
        if mem:
            f.write('%mem={mem}\n'.format(mem=mem))
        if wfx:
            f.write('#p {esm}/{bas} output=wfx nosymmetry\n'.format(esm=esm,bas=basis_set))
        else:
            f.write('#p {esm}/{bas} nosymmetry\n'.format(esm=esm,bas=basis_set))
        f.write('\n')
        f.write('single point on {file} reoriented by subproptools\n'.format(file=old_file_name))
        f.write('\n')
        f.write('{q} {mul}\n'.format(q=charge,mul=multiplicity))
        dfAsString = reor_geom.to_string(header=False, index=False)
        f.write(dfAsString)
        f.write('\n\n')
        if wfx:
            f.write(new_file_name+'.wfx\n\n\n')
        else:
            f.write('\n\n\n')    
    return

def rotate_sheet(csv_file,esm,basis,n_procs=4,mem='3200MB',wfx=True,extra_label=''):
    """
    Given csv file and Gaussian calculation options, reorient files in csv_file and output gjf

    Args:
        csv_file - csv file containining these columns:
        Substituent, originAtom, negXAtom, posYAtom, charge, multiplicity,  label1, label2,...
            Substituent: string label for substituent
            originAtom: numerical index(starting form 1) of atom to use as origin
            negXAtom: numerical index(starting form 1) of atom to place on -x
            posYAtom: usually 0, but if override desired, numerical index(starting form 1) of atom to place on +y
            charge: charge of the molecule
            multiplicity: multiplicity of the molecule
            label1, label2... label depicting situation for molecule(substrate, method)
                e.g. "SubH", "B3LYP/cc-pvDZ" etc
        esm: electronic structure method (HF/MP2/DFT functional/etc)
        basis: string for basis to be used
        n_procs: number of processors for use on Cedar, default to 4
        mem: memory to use on Cedar, default to 3200MB
        wfx: if we wish to write wfx, default True
        extra_label: an additional label for the reoriented file if needed, default none
    Returns:
        no return value, but outputs to gjf files in working directory(or directory in path of filenames)    

    """
    csvFrame = pd.read_csv(csv_file)
    ncolumns = csvFrame.shape[1]
    nrow = csvFrame.shape[0]
    for sub in range(0,nrow):
        charge=csvFrame['charge'][sub]
        multiplicity=csvFrame['multiplicity'][sub]
        originAtom = csvFrame['originAtom'][sub]
        negXAtom = csvFrame['negXAtom'][sub]
        posYAtom = csvFrame['posYAtom'][sub]
        for label in range(6,ncolumns):
            rot_file = csvFrame.iloc[sub,label]
            rot_geom = rotate_substituent(rot_file,originAtom=originAtom,negXAtom=negXAtom,posYAtom=posYAtom)
            output_to_gjf(rot_file,rot_geom,esm=esm,basis_set=basis,add_label=extra_label,
                         n_procs=n_procs,mem=mem,charge=charge,multiplicity=multiplicity,wfx=wfx)
    return
