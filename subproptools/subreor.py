import os # file system stuff
import sys
sys.path.append(sys.path[0].replace('subproptools','')+'/'+'referenceMaps')
import pandas as pd #data frames
import math #sqrt
import numpy as np #arrays
from subproptools import qtaimExtract as qt #sum file manipulation and property extraction
import regex as re

def get_bcp_reference(originAtom,numBonds,refDict):
    """Takes reference dictionary containing reference maps and chooses the right map to use."""
    if originAtom =='C' and numBonds==4: #sp3 carbon
        retDict = refDict['C']['sp3']
    elif originAtom =='C' and numBonds==3: #sp2 carbon
        retDict = refDict['C']['sp2']
        #note linear shouldn't get to this point
    elif originAtom =='B' and numBonds==3: #planar boron
        print('planar boron')
        retDict = refDict['B']['sp2']
    elif originAtom =='B' and numBonds==4: #sp3 boron
        print('sp3 boron')
        retDict = refDict['B']['sp3']
    elif originAtom == 'N' and numBonds==4:
        print('ammonium')
        retDict = refDict['N']['sp3']
    elif originAtom =='Al' and numBonds==3: #planar boron
        print('planar aluminum')
    elif originAtom =='Al' and numBonds==4: #sp3 boron
        print('sp3 aluminum') 
    elif originAtom =='Si' and numBonds==4: #sp3 carbon
        print('sp3 silicon')
        retDict = refDict['Si']['sp3']
    elif originAtom =='Si' and numBonds==3: #sp2 carbon
        print('sp2 silicon')
        retDict = refDict['Si']['sp2']
    elif originAtom == 'P' and numBonds==4:
        print('phosphonium')
        retDict = refDict['P']['sp3']
    return retDict 

def default_stat_dict():
    """defines mean and sd for scaling properties if no other available
    
    Scaling based on mean & sd of all bcps of files used in transferability paper(not just Substrate-Substituent)
    See:notebooks/bcpStatExtraction.ipynb
    Returns:
      a dictionary for each property, with names matching the names returned in a get_bcp_properties dictionary
    """
    stat_dict = {'rho': {'mean':0.290686,'sd':0.077290},
                 'lambda1':{'mean':-0.725552,'sd':0.299756},
                 'lambda2':{'mean':-0.678830,'sd':0.291123},
                 'lambda3':{'mean':0.583261,'sd':0.449474},
                 'DelSqRho':{'mean':-0.821120,'sd':0.570553},
                 'Ellipticity':{'mean':0.077722,'sd':0.137890},
                 'V':{'mean':-0.475963,'sd':0.332327},
                 'G':{'mean':0.135342,'sd':0.184923},
                 'H':{'mean':-0.340622,'sd':0.176627},
                 'DI(R,G)':{'mean':1.081894,'sd':0.369324},}
    return stat_dict

def find_bcp_match(data,originAtomXYZ,negXAtomLabel, originAtomLabel):
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
    originBCPs = qt.find_connected(data,negXAtomLabel,originAtomLabel)
    bcpPropDict = {}
    #get the bcp properties
    for bcp in originBCPs:
        bcpBlock = qt.get_bcp_block(data,atPair=bcp)
        bcpPropDict.update({bcp[0]+'-'+bcp[1]: qt.get_bcp_properties(bcpBlock)})
    clockwiseKeys = find_clockwise_rot(bcpPropDict,originAtomXYZ,originAtomLabel)
    #at this point have bcpDictionary ordered from 1st to last with clockwise bcp
    return clockwiseKeys #this is a dictionary of bcps

def angle_btw_bcp(xyzA,xyzB):
    """find angle between two bcp orientation vectors defined by x,y,z np.array, after flattening to yz plane."""
    xyzA[0] = 0.0
    xyzB[0]=0.0
    angle= math.acos((xyzA[0]*xyzB[0]+xyzA[1]*xyzB[1]+xyzA[2]*xyzB[2])/(math.sqrt(xyzA[0]**2+xyzA[1]**2+xyzA[2]**2)*math.sqrt(xyzB[0]**2+xyzB[1]**2+xyzB[2]**2)))
    return angle

def find_clockwise_rot(bcpPropDict,originAtomXYZ,originAtomLabel):
    #given dictionary of bcp properties, find which ones are in a clockwise rotation
    #return list of dictionary keys ordered for clockwise rotation
    crossDict = {}
    for key1 in bcpPropDict:
        for key2 in bcpPropDict:
            if key1 != key2:
                xyz1 = bcpPropDict[key1]['xyz']
                xyz1[0] = 0.#project to yz plane
                xyz2=bcpPropDict[key2]['xyz']
                xyz2[0]=0.
                cross=np.cross(bcpPropDict[key1]['xyz'],bcpPropDict[key2]['xyz'])[0]
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
#         keyDict.update({'AngleToY': angle_btw_bcp(xyz,np.array([0.0,1.0,0.0]))})
#         keyDict.update({'AngleToZ': angle_btw_bcp(xyz,np.array([0.0,0.0,1.0]))})
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

def set_origin(xyzArray,originAtom):
    """shifts all points in xyz geometry by originAtom coordinates, returns xyz geometry"""
    org = xyzArray[originAtom-1,]
    return xyzArray - org

def get_lengths(xyzArray):
    """given xyz geometry returns np.array of magnitudes of squared distances from the origin"""
    lengths = np.array([])
    for atom in xyzArray:
        length=0
        for coord in atom:
            length += coord**2
        lengths = np.append(lengths,length)
    return lengths

def zero_y_for_negx(t_xyz,negXAtom):
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

def zero_z_for_negx(t_xyz,negXAtom):
    """perform second rotation in setting -x axis to zero the z value"""
    #perform after y
    d = t_xyz[0,negXAtom-1]/math.sqrt(t_xyz[0,negXAtom-1]**2 + t_xyz[2,negXAtom-1]**2)
    s = t_xyz[2,negXAtom-1]/math.sqrt(t_xyz[0,negXAtom-1]**2 + t_xyz[2,negXAtom-1]**2)
    G = np.array([[d,0,s],[0,1,0],[-s,0,d]])
    return np.matmul(G,t_xyz)

def set_xaxis(xyzArray,negXAtom):
    """Given xyz geometry with atom at origin, return xyz geometry with negXAtom on -x."""
    t_xyz = xyzArray.T
    
    #define initial xyz vector lengths. Should be unchanged after rotation
    tol = 0.0001 #tolerance for change
    initial_lengths= get_lengths(xyzArray)

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
    t_rot1 = zero_y_for_negx(t_xyz,negXAtom)
    rot1_lengths = get_lengths(t_rot1.T)
    
    if np.any((rot1_lengths - initial_lengths)>tol):
        print("Geometry perturbation exceeded")
    t_rot2 = zero_z_for_negx(t_rot1,negXAtom)
    rot2_lengths = get_lengths(t_rot2.T)
    
    if np.any((rot2_lengths - initial_lengths)>tol):
        print("Geometry perturbation exceeded")
    if t_rot2[0,negXAtom-1] > 0: #if negxatom is along +x, rotate 180
        G = np.array([[-1,0,0],[0,1,0],[0,0,-1]])
        t_rot_final = np.matmul(G,t_rot2)
        return t_rot_final.T
    else:
        return t_rot2.T
    
def align_dicts(testDict,refDict,statDict):
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
    dif00 = get_popelier_dif(testDict[testKeysList[0]],refDict[refDictKeysList[0]],statDict)
    dif11 = get_popelier_dif(testDict[testKeysList[1]],refDict[refDictKeysList[1]],statDict)
    dif22 = get_popelier_dif(testDict[testKeysList[2]],refDict[refDictKeysList[2]],statDict)
    #BCP distances for second match
    dif10 = get_popelier_dif(testDict[testKeysList[1]],refDict[refDictKeysList[0]],statDict)
    dif21 = get_popelier_dif(testDict[testKeysList[2]],refDict[refDictKeysList[1]],statDict)
    dif02 = get_popelier_dif(testDict[testKeysList[0]],refDict[refDictKeysList[2]],statDict)
    #BCP distances for third match
    dif20 = get_popelier_dif(testDict[testKeysList[2]],refDict[refDictKeysList[0]],statDict)
    dif01 = get_popelier_dif(testDict[testKeysList[0]],refDict[refDictKeysList[1]],statDict)
    dif12 = get_popelier_dif(testDict[testKeysList[1]],refDict[refDictKeysList[2]],statDict)
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

def get_posy_point(sumFileNoExt,atomDict,attachedAtom,negXAtomLabel,default_stats=True):
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
        print(posYPoint)
        #reorient to average of vscc points for +y
    else:
        refDict = get_reference_map()
        sumFile = open(sumFileNoExt+".sum",'r')
        data = sumFile.readlines()
        sumFile.close()
        #bcpsToMatch is bcp dictionary, ordered for clockwise rot
        #data,originAtomXYZ,negXAtomLabel,originAtomLabel
        bcpsToMatch = find_bcp_match(data,atomDict[attachedAtom]['xyz'],negXAtomLabel,attachedAtom)
        numBonds=len(bcpsToMatch)+1
        atType = ''.join([i for i in attachedAtom if not i.isdigit()])
        matchDict = get_bcp_reference(atType,numBonds,refDict)
        if default_stats:
            statDict = default_stat_dict()
        posYPoint = align_dicts(bcpsToMatch,matchDict,statDict)
        #reorient to another point
        #posy point will be the point that would lie along the y-axis in reference in maximal match case
        print('not done yet')
    return posYPoint    

def set_yaxis(xyzArray,posYArray):
    """rotate a geom positioned on -x axis so posYArray will lie on +y"""
    theta = math.atan2(posYArray[2],posYArray[1])
    print(theta)
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


def get_popelier_dif(bcpDictA,bcpDictB,statDict):
    """compute distance in BCP space between dictA and B scaling using statDict"""
    distancesq = 0.
    for prop in bcpDictA:
        if prop != 'xyz':
            scaledA = (bcpDictA[prop][0]-statDict[prop]['mean'])/statDict[prop]['sd']
            scaledB = (bcpDictB[prop][0]-statDict[prop]['mean'])/statDict[prop]['sd']
            distancesq += (scaledA - scaledB)**2
    return math.sqrt(distancesq)    


def get_ref_bcps(sumfilenoext,atPairList,originAtom,originAtomXYZ=np.array([0.,0.,0.])):
    """given reference sumfile, extract bcp properties for needed bcps"""
    sumFile = open(sumfilenoext+".sum","r") #open file, read lines, close file
    data = sumFile.readlines()
    sumFile.close()
    bcpDict = {}
    for bcp in atPairList:
        block = qt.get_bcp_block(data,bcp)
        bcpDict.update({'{at1}-{at2}'.format(at1=bcp[0],at2=bcp[1]):qt.get_bcp_properties(block)})
    clockbcpDict = find_clockwise_rot(bcpDict,originAtomXYZ,originAtom)    
    return clockbcpDict  


def get_reference_map():
    """creates dictionary with reference species to map coordinate system to"""
    pathToRefs=sys.path[0].replace('subproptools','')+'/'+'referenceMaps'+'/'
    print(pathToRefs)
    refDict={}
    carbonDict = {}
    carbonDict.update({'sp3':get_ref_bcps(pathToRefs + 'SubH_CFHCH3_reor_wb97xd_aug-cc-pvtz',[['C1','H3'],['C1','C4'],['C1','F8']],originAtom='C1')})
    carbonDict.update({'sp2':get_ref_bcps(pathToRefs +'SubH_CHCH2_wb97xd_aug-cc-pvtz_reor',[['C1','H3'],['C1','H2']],originAtom='C1')})
    #carbonDict.update({'sp2':get_ref_bcps('SubH_CHCH2-ReorJuly2-B3LYP-def2-TZVPPD-Field',[['C1','H3'],['C1','C4']],'C1')})
    refDict.update({'C':carbonDict})
    boronDict = {}
    boronDict.update({'sp3':get_ref_bcps(pathToRefs +'SubH_BHCH3BH2_wb97xd_aug-cc-pvtz_reor',[['B1','H3'],['B1','H4'],['B1','C5']],originAtom='B1')})
    boronDict.update({'sp2':get_ref_bcps(pathToRefs +'SubH_BHF_wb97xd_aug-cc-pvtz_reor',[['B1','H3'],['B1','F4']],originAtom='B1')})
    refDict.update({'B':boronDict})
    nitrogenDict = {}
    nitrogenDict.update({'sp3':get_ref_bcps(pathToRefs +'SubH_NHOCH3_wb97xd_aug-cc-pvtz_reor',[['N1','H3'],['N1','O4'],['N1','C5']],originAtom='N1')})
    refDict.update({'N':nitrogenDict})
    phosphorousDict = {}
    phosphorousDict.update({'sp3':get_ref_bcps(pathToRefs +'SubH_POCH3H_wb97xd_aug-cc-pvtz_reor',[['P1','H3'],['P1','O4'],['P1','C5']],originAtom='P1')})
    refDict.update({'P':phosphorousDict})
    siliconDict = {}
    siliconDict.update({'sp3':get_ref_bcps(pathToRefs +'SubH_SiFHCH3_wb97xd_augccpvtz_reor',[['Si1','H3'],['Si1','C4'],['Si1','F8']],originAtom='Si1')})
    siliconDict.update({'sp2':get_ref_bcps(pathToRefs +'SubH_Si2H3_wb97xd_augccpvtz_reor',[['Si1','H3'],['Si1','NNA7']],originAtom='Si1')})
    refDict.update({'Si':siliconDict})
    #update dictionary with further references
    return refDict


def rotate_substituent(sumFileNoExt,originAtom,negXAtom,posYAtom=0):
    """
    Rotates a substituent to the defined coordinate system.
    Coordinate system defined as: 
    originAtom at (0,0,0)
    negXAtom at (-x,0,0)
    Atom on +y defined as:
        average of lone pairs for 2 lone pairs on originAtom
        Position of lone pair for 1 lone pair on originAtom
        For no lone pairs: map BCPs onto reference for the atom type
        Minimum distance in BCP space defined the atom to put on +y

    Args:
        sumFileNoExt (string): name of a sum file, without the .sum extension
        originAtom (int): the integer number of the atom to place at the origin
        negXAtom (int): the integer number of the atom to place along the -x axis

    Returns:
        pandas data frame of output geometry (columns Atom, x, y, z)
        Outputs geometry to .txt file of name sumFileNoExt_reor.txt
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
    molecule_orig = set_origin(molecule_xyz['xyz'],originAtom)
    molecule_xaxis = set_xaxis(molecule_orig,negXAtom)
    if posYAtom:
        posYPoint = molecule_xaxis[posYAtom-1] 
    else:    
        posYPoint = get_posy_point(sumFileNoExt,atomDict,attachedAtom,negXAtomLabel)
    final_orientation = set_yaxis(molecule_xaxis,posYPoint)
    #Generate output
    outFrame = pd.DataFrame(final_orientation*0.529177,columns = ['x','y','z'])
    outFrame['Atom'] = molecule_xyz['Atoms']
    outFrame = outFrame[['Atom','x','y','z']]
    with open(sumFileNoExt+'_reor.txt','w') as txt_file:
        of_string = outFrame.to_string(header=False,index=False)
        txt_file.write(of_string)
    return outFrame

def outputToGjf(old_file_name,reor_geom,esm='wB97XD',basis_set="aug-cc-pvtz",add_label='',n_procs=4,mem='3200MB',charge=0,multiplicity=1,wfx=True):
    new_file_name =old_file_name+'_reor'+add_label
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
        f.write('\n')
        f.write(new_file_name+'.wfx\n\n\n')
    return
