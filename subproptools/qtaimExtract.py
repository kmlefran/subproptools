import os # file system stuff
import pandas as pd #data frames
import math #sqrt
import numpy as np #arrays
def search_str(linesObj,word,searchStart=0,ignore=0):
    #lines obj - list, each containing a string. e.g. the output of file.readlines() method
    #word - string that we look to match
    #searchStart - integer, if we don't need to start searching at the beginning, choose \
    #line number to start at
    #ignore - integer on how many instances of the string to ignore 
    #(e.g. if you need the second instance, ignore=1)
    #Returns index of line that word occurs, or -1 if not found
    wordLine=-1 #will return -1 if string not found
    #start for
    for ln_num, line in enumerate(linesObj): #iterate over lines
        #print(line)
        
        if line.find(word) >-1 and linesObj.index(line) >= searchStart:
        # start outer if - True if line found past the searchStart    
            if ignore == 0:
            #Start inner if - if this is the instance of the string you want, get index    
                wordLine = ln_num
                break
            else:
            #if not instance of string you want, deduct from the count    
                ignore = ignore - 1
                continue
            #end inner if
       #end outer if     
    #end for            
    return wordLine

def get_di_table(data,tableHeader='2*D2(A,B)          DI(A,B)          %Deloc(A,B)       %Deloc(B,A)',endString='Areas of Interatomic Surfaces:'):
    tableStart = search_str(data,tableHeader,searchStart=0)+2 #data starts 2 lines after header
    tableEnd = search_str(data,endString,tableStart)-2 #-2 since there is a 2 line gap
    tableData=data[tableStart:tableEnd]
    headerNames = ['Atom A','Atom B']
    for colName in tableHeader.split():
        headerNames.append(colName) 
    table = pd.DataFrame(columns=headerNames)
    #Add each row of table to Data Frame
    for i,line in enumerate(tableData):
        table.loc[i] = line.split()[0:6]
    return table        

def get_table(data,tableHeader,ignored=0,endString='Total'):
    #data - list each containing a string. i.e. the output of sumfile.readlines() method
        #note this function is designed for .sum files and won't work on others
    #tableHeader - the header of the table, NOT including "Atom" or "Atom A"
    #ignored=0 - to be passed to search_str, find the (ignored)+1 instance of the string
    #endString - a string two past the last row where data in the table is
    #For most .sum file tables, this is Total, but for some (quadrupoles)
    #, there is no Total, so it will usually take the name of the next table
    #returns a DataFrame with Table matching the header
    #Find start and end lines of table
    #print(data[0])
    tableStart = search_str(data,tableHeader,searchStart=0,ignore=ignored)+2 #data starts 2 lines after header
    tableEnd = search_str(data,endString,tableStart)-1 #only deduct 1, because next line will get up to but not including that line
    #So the -1 here and the next line combine to get data up to two rows before endString
    tableData=data[tableStart:tableEnd] #get data
    #Generate column names of DataFrame. #Include atom here, as sumfile header is Atom A and the whitespace there will mess with how code is writting
    headerNames = ['Atom']
    for colName in tableHeader.split():
        headerNames.append(colName)
    #Initialize empty data frame    
    table = pd.DataFrame(columns=headerNames)
    #Add each row of table to Data Frame
    print(headerNames)
    
    for i,line in enumerate(tableData):
        #print(line)
        table.loc[i] = line.split()
        
    return table    

def get_bcp_block(data,atPair=["C1","H2"]):
    #find all the lines with BCP data for a pair of atoms
    #data - list each containing a string. i.e. the output of sumfile.readlines() method
        #note this function is designed for .sum files and won't work on others
    #atPair - list of length 2, containing atom element and numeric label, see default for ex    
    #returns list of strings with BCP data    
    #find start line of BCP data, search both permutations of atom labels    
    #Start -1 from the search line, because the coords are the line before
    bcpStart = search_str(data,word='Type = (3,-1) BCP {at1} {at2}'.format(at1=atPair[0],at2=atPair[1]))-1
    if bcpStart == -2:
        print('didnt find first')
        bcpStart = search_str(data,word='Type = (3,-1) BCP {at1} {at2}'.format(at1=atPair[1],at2=atPair[0]))-1
    bcpEnd = bcpStart+34 #search_str(data,word='GBL_IV',searchStart=bcpStart) #find end of bcp data at GBL_IV
    bcpBlock = data[bcpStart:bcpEnd]
    #print(bcpBlock)
    return bcpBlock # return the lines of the BCP data

def get_sub_di(data,subAtomLabels=[]):
    diTable = get_di_table(data)
    diTable = diTable.drop(['2*D2(A,B)','%Deloc(A,B)','%Deloc(B,A)'],axis=1)
    diClass = []
    #logic error - need to instead just iterate over tables and see if equals any of subAtomLabels for A/B
    #after meeting
    diTable['DI(A,B)'] = diTable['DI(A,B)'].astype(float)
    for ind in diTable.index:
        if any(x == diTable['Atom A'][ind] for x in subAtomLabels) and any(x == diTable['Atom B'][ind] for x in subAtomLabels):
            diClass.append('Substituent-Substituent')
        elif any(x == diTable['Atom A'][ind] for x in subAtomLabels) or any(x == diTable['Atom B'][ind] for x in subAtomLabels):
            diClass.append('Substituent-Substrate')
        else:
            diClass.append('Substrate-Substrate')
    diTable['Interaction'] = diClass
    return sum(diTable.loc[diTable['Interaction'] == 'Substituent-Substrate', 'DI(A,B)'])

def get_bcp_properties(bcpBlock):
    #bcpBlock: takes output of get_bcp_block
    #initialize empyt dict
    bcpDict={}
    for line in bcpBlock: #iterate over lines
        splitLine = line.split() #split the line into individual words based on whitespace
        #note splitLine[i] - i is from manual inspection of data
        if 'Coords' in splitLine: #if line contains the coordinates, get coordinates
            bcpDict.update({'xyz': np.array([float(splitLine[4]),float(splitLine[5]),float(splitLine[6])])})
            #bcpDict.update({'y': splitLine[5]})
            #bcpDict.update({'z': splitLine[6]})
        elif 'Rho' in splitLine: #get density at BCP
            bcpDict.update({'rho': [float(splitLine[2])]})
        elif 'HessRho_EigVals' in splitLine: #get lambda_i
            bcpDict.update({'lambda1': [float(splitLine[2])]})
            bcpDict.update({'lambda2': [float(splitLine[3])]})
            bcpDict.update({'lambda3': [float(splitLine[4])]})
        elif 'DelSqRho' in splitLine: #get DelSqRho
            bcpDict.update({'DelSqRho': [float(splitLine[2])]})
        elif 'Ellipticity' in splitLine: #get ellipticity
            bcpDict.update({'Ellipticity': [float(splitLine[3])]})
        elif 'V' in splitLine: #get V
            bcpDict.update({'V': [float(splitLine[2])]})
        elif 'G' in splitLine: # get G
            bcpDict.update({'G': [float(splitLine[2])]})
        #elif line.find('Rho') 
    bcpDict.update({'H': [bcpDict['V'][0] + bcpDict['G'][0]]})  #get total energy density  
    return bcpDict

def get_atomic_props(data):
    #data - - list each containing a string. i.e. the output of sumfile.readlines() method
        #note this function is designed for .sum files and won't work on others
    #returns a dictionary of dictionaries. The keys for the outer dictionary are atom labels
    #the keys for the inner dictionary are properties
    #e.g. can get data as allatomDict['C1']['Mu_X'], which gets x dipole for atom C1
    
    #find tables where desired data is located    
    xyzFrame = get_table(data,'Charge                X                  Y                  Z',endString='Some Atomic Properties:')
    eFrame = get_table(data,'q(A)              L(A)              K(A)          K_Scaled(A)      |Mu_Intra(A)|')
    muXFrame = get_table(data,'Mu_Intra_X(A)     Mu_Bond_X(A)        Mu_X(A)')
    muYFrame = get_table(data,'Mu_Intra_Y(A)     Mu_Bond_Y(A)        Mu_Y(A)')
    muZFrame = get_table(data,'Mu_Intra_Z(A)     Mu_Bond_Z(A)        Mu_Z(A)')
    muMagFrame = get_table(data,'|Mu_Intra(A)|     |Mu_Bond(A)|        |Mu(A)|')
    quadFrame = get_table(data,'Q_XX(A)           Q_XY(A)           Q_XZ(A)           Q_YY(A)           Q_YZ(A)           Q_ZZ(A)',endString='Eigenvalues and Eigenvectors of Atomic Traceless Quadrupole Moment Tensors')
    #need these radial distoritions for calculating atomic quadrupole contributions
    radFrame = get_table(data,'R-2(A)            R-1(A)            R0(A)             R+1(A)            R+2(A)',endString='Atomic Radial Distortion Moments')
    #currently volume is getting the 0.0004 isosurface, and i don't know how to change that
    #like in theory this should work, but it's not. I'll fix later
    volFrame = get_table(data,'Area(A)           Vol(A)          N(Vol(A))      N(Vol(A))/Vol(A)     %N(Vol(A))',ignored=1)
    
    #atomicFrame = pd.concat([eFrame,muXFrame,muYFrame,muZFrame,quadFrame,volFrame,radFrame],axis=1)
    allatomDict = {} #initialize empty dictionary to store dictionaries of the atoms
    for atom in eFrame['Atom']:
        atomDict = {} #create a dictionary for each atom
        #the following lines go through the tables, and find the row and column with desired property
        # the row matches Atom, the column matches property we want, and we use .iloc[0] to get it as a string, which we then change to float. 
        #If no iloc[0] it is a data frame object and float will eventually return an error on it
        #xyz is a np array - element 0 is x, element 1 is y, element 2 is z
        atomDict.update({'xyz': np.array([float(xyzFrame[xyzFrame['Atom'] == atom]['X'].iloc[0]),float(xyzFrame[xyzFrame['Atom'] == atom]['Y'].iloc[0]),float(xyzFrame[xyzFrame['Atom'] == atom]['Z'].iloc[0])])})
        atomDict.update({'q': [float(eFrame[eFrame['Atom'] == atom]['q(A)'].iloc[0])]})
        atomDict.update({'K': [float(eFrame[eFrame['Atom'] == atom]['K(A)'].iloc[0])]})
        atomDict.update({'K_Scaled': [float(eFrame[eFrame['Atom'] == atom]['K_Scaled(A)'].iloc[0])]})
        atomDict.update({'Mu_Intra_X': [float(muXFrame[muXFrame['Atom']==atom]['Mu_Intra_X(A)'].iloc[0])]})
        atomDict.update({'Mu_Intra_Y': [float(muYFrame[muYFrame['Atom']==atom]['Mu_Intra_Y(A)'].iloc[0])]})
        atomDict.update({'Mu_Intra_Z': [float(muZFrame[muZFrame['Atom']==atom]['Mu_Intra_Z(A)'].iloc[0])]})
        atomDict.update({'Mu_Bond_X': [float(muXFrame[muXFrame['Atom']==atom]['Mu_Bond_X(A)'].iloc[0])]})
        atomDict.update({'Mu_Bond_Y': [float(muYFrame[muYFrame['Atom']==atom]['Mu_Bond_Y(A)'].iloc[0])]})
        atomDict.update({'Mu_Bond_Z': [float(muZFrame[muZFrame['Atom']==atom]['Mu_Bond_Z(A)'].iloc[0])]})
        atomDict.update({'Mu_X': [float(muXFrame[muXFrame['Atom']==atom]['Mu_X(A)'].iloc[0])]})
        atomDict.update({'Mu_Y': [float(muYFrame[muYFrame['Atom']==atom]['Mu_Y(A)'].iloc[0])]})
        atomDict.update({'Mu_Z': [float(muZFrame[muZFrame['Atom']==atom]['Mu_Z(A)'].iloc[0])]})
        atomDict.update({'|Mu_Intra|': [float(muMagFrame[muMagFrame['Atom']==atom]['|Mu_Intra(A)|'].iloc[0])]})
        atomDict.update({'|Mu_Bond|': [float(muMagFrame[muMagFrame['Atom']==atom]['|Mu_Bond(A)|'].iloc[0])]})
        atomDict.update({'|Mu|': [float(muMagFrame[muMagFrame['Atom']==atom]['|Mu(A)|'].iloc[0])]})
        atomDict.update({'Q_XX': [float(quadFrame[quadFrame['Atom']==atom]['Q_XX(A)'].iloc[0])]})
        atomDict.update({'Q_XY': [float(quadFrame[quadFrame['Atom']==atom]['Q_XY(A)'].iloc[0])]})
        atomDict.update({'Q_XZ': [float(quadFrame[quadFrame['Atom']==atom]['Q_XZ(A)'].iloc[0])]})
        atomDict.update({'Q_YY': [float(quadFrame[quadFrame['Atom']==atom]['Q_YY(A)'].iloc[0])]})
        atomDict.update({'Q_YZ': [float(quadFrame[quadFrame['Atom']==atom]['Q_YZ(A)'].iloc[0])]})
        atomDict.update({'Q_ZZ': [float(quadFrame[quadFrame['Atom']==atom]['Q_ZZ(A)'].iloc[0])]})
        atomDict.update({'R+2': [float(radFrame[radFrame['Atom']==atom]['R+2(A)'].iloc[0])]})
        atomDict.update({'R+1': [float(radFrame[radFrame['Atom']==atom]['R+1(A)'].iloc[0])]})
        atomDict.update({'Vol': [float(volFrame[volFrame['Atom']==atom]['Vol(A)'].iloc[0])]})
        allatomDict.update({atom: atomDict}) #add the atomic dictionary to the dictionary of dictionaries    
    return allatomDict

def get_dist(coordListA,coordListB):
    #takes 2 points, coordList are 3 element list/array of xyz coords.
    #returns distance between them
    dist = math.sqrt((coordListA[0]-coordListB[0])**2 + (coordListA[1]-coordListB[1])**2 + (coordListA[2]-coordListB[2])**2)
    return dist

def find_closest_nuclei(xyz,atomDict):
    #probably need to update in symmetry case in case of equidistant nuclei
    # xyz - takes a point as a 3 length list of xyzcoords
    # atomDict - object from get_atomic_props
    #returns the two atoms that are closest to the xyz point given
    # original design purpose was to find nuclei closest to a charge concentration
    distList=[] #initialize empty list for distances
    atList = [] #initialize empty list for atom labels
    for atom in atomDict: #get list of distances from point xyz to all atoms, and list of atom labels
        atList.append(atom)
        distList.append(get_dist(xyz,atomDict[atom]['xyz']))
    npDistList= np.array(distList) #convert distance list to array for np.where functionality
    minDist = np.partition(npDistList,1)[0:2] #find the two lowest distances
    #return atList with indices of npDistList that are equal to the two lowest values
    return [atList[np.where(npDistList==minDist[0])[0][0]],atList[np.where(npDistList==minDist[1])[0][0]]]

def identify_vscc(multiccDict,atomDict, thresh=0.7):
    #multiccDict - dictionary of charge concentrations - output of get_cc_props
    #atomDict - dictionary of atomic properties - output of get_atomic_props
    #thresh - numeric threshold for distance between shells 
    #(eg inner shell charge concentration is more than 0.7 au closer to nuclei than VSCC)
    
    #takes a list of all charge concentrations
    #filters it to return dictionary of only points identified as VSCC based on below criteria
    #will be a vscc if its the farthest out, if its not on a bond line, and if it is not close to nuclei
    
    vsccDict={} #initialize empty dictionary for cc
    nucDistList = [] #empty list to store distances
    potentialccList = [] #empty list to store keys for potential VSCC after eliminating some criteria
    for cc in multiccDict: #for each cc
        if multiccDict[cc]['distFromNuc'] > 0.1: # this will filter out the innermost cc that are on the nuclei
            ccXYZ = multiccDict[cc]['xyz'] #get the xyz coords of the cc
            nucPair = find_closest_nuclei(ccXYZ,atomDict) #find the two closest atoms to the cc
            #if it is a bonded cc, it will be the nuclei to which it is bonded
            #check if the cc is on the bond between those two nuclei
            isBonded = is_on_line(atomDict[nucPair[0]]['xyz'],atomDict[nucPair[1]]['xyz'],ccXYZ)
            if not isBonded: #if it is not on the line, it is potentially a VSCC, store it, and its distance
                nucDistList.append(multiccDict[cc]['distFromNuc'])
                potentialccList.append(cc)
    #at this point, non-bonded, non-core CC remain, so just need to check which CC are the outermost
    #do this by comparing to maximum distance from nuclei
    if len(nucDistList) > 0:
        outerShellDist = max(nucDistList) #note this only checks nonbonded noncore
    for cc in potentialccList:
        if abs(multiccDict[cc]['distFromNuc']-outerShellDist) < thresh: 
            #for given cc, if close (within thresh) to outermost, it is a VSCC, store it as such
            #So, going up to 3rd row for now. May work for 4th but I haven't checked
            #C/O/F easy - any not core not bonded should be VSCC
            #P/S/Si harder
            # in example data inner~0.25-0.3, outer ~ 1.2-1.3
            #arbitrary default threshold seems like 0.7 will work
            vsccDict.update({cc: multiccDict[cc]})
    return vsccDict    

def is_on_line(lineStartList,lineEndList,pointToCheckList,epsilon=0.1):
    #return boolean if pointToCheck is on bond between lineStart and lineEnd
    #takes 3 lists of coordinates of 3 length xyz
    #line connecting atoms: (x1 + t(x2-x1),y1 + t(y2-y1),z1 + t(z2-z1))
    #find t - written out on a piece of paper somewhere based on dot product of vectors
    #to reconstruct - create equations of lines. (pointToCheck-pointWitht) and (lineStart-pointWitht)
    # dot product of those vectors should be 0(closest point is. 90 degrees)
    #rearrange to solve for t
    t = ((pointToCheckList[0]-lineStartList[0])*(lineEndList[0]-lineStartList[0]) + 
        (pointToCheckList[1]-lineStartList[1])*(lineEndList[1]-lineStartList[1]) + 
        (pointToCheckList[2]-lineStartList[2])*(lineEndList[2]-lineStartList[2]))/(((lineEndList[2]-lineStartList[2])**2) + ((lineEndList[1]-lineStartList[1])**2)+((lineEndList[0]-lineStartList[0])**2))
    #Find the point defined by t
    perpPoint = np.array([lineStartList[0]+t*(lineEndList[0]-lineStartList[0]),
                          lineStartList[1]+t*(lineEndList[1]-lineStartList[1]),
                          lineStartList[2]+t*(lineEndList[2]-lineStartList[2])])
    #get distance between test point and closest line point
    distToLine = get_dist(pointToCheckList,perpPoint) 
    #define thresholds to check if t lies between nuclei or beyond.
    #note - can I just check if t < 0 or t > 1 for that?
    #leaving as is for now
    lowX = min(lineStartList[0],lineEndList[0])
    highX = max(lineStartList[0],lineEndList[0])
    lowY = min(lineStartList[1],lineEndList[1])
    highY = max(lineStartList[1],lineEndList[1])
    lowZ = min(lineStartList[2],lineEndList[2])
    highZ = max(lineStartList[2],lineEndList[2])
    #if t is between the end points of the line, and the pointToCheck is within epsilon of that point, it is a bonded cc
    #return boolean based on evaluation
    if perpPoint[0] >= lowX and perpPoint[0] <= highX and perpPoint[1] >= lowY and perpPoint[1] <= highY and perpPoint[2] >= lowZ and perpPoint[2] <= highZ and distToLine < epsilon:
        return True
    else:
        return False
    
def get_cc_props(filename,atomLabel,type='vscc'):
    #takes a sumfilename with no extension and an atom label that we want cc for (eg F4)
    #returns dictionary of all cc for atomLabelwith all the (3,+3) cc and properties
    #create path to subdirectory
    currdir = os.getcwd()
    #example dir name SubH_CCF-ReorPosY-B3LYP-def2-TZVPPD-Field_atomicfiles
    subdirname = filename + "_atomicfiles"
    pathToSubdir = currdir+ "/"+subdirname+"/"
    lowerAtomLabel = atomLabel.lower() #the laprhocp files are eg f4.agpviz - note lower case
    #open the file, read data, and close file
    atFile = open(pathToSubdir+lowerAtomLabel+".agpviz",'r')
    atData = atFile.readlines()
    atFile.close()
    #initialize empty dict to store all cc
    allccDict = {}
    alldictcounter=1 #counter since the key will be numeric
    
    for lnNum,line in enumerate(atData):
        if "Type = (3,+3)" in line: #if we're looking at CC
            oneccDict = {} # create empty dict for storing properties for one cc
            xyzSplit = atData[lnNum+1].split() #line after label is xyz, split it and store xyz
            oneccDict.update({'xyz':np.array([float(xyzSplit[2]),float(xyzSplit[3]),float(xyzSplit[4])])})
            #next line is distance from nuc, split and store
            distFromNucSplit = atData[lnNum+2].split()
            oneccDict.update({'distFromNuc': float(distFromNucSplit[2])})
            #next is rho, then delsqrho - split and store those
            rhoSplit = atData[lnNum+3].split()
            oneccDict.update({'rho': float(rhoSplit[2])})
            delsqrhoSplit = atData[lnNum+4].split()
            oneccDict.update({'delsqrho': float(delsqrhoSplit[2])})
            #add the cc to the dictionary, update counter
            allccDict.update({alldictcounter: oneccDict})
            alldictcounter += 1
    #find all (3,+3) cps and get their properties
    return allccDict

def get_atomic_quad_contrib(atomDict):  
    #atomDict - one of the dictionary object returned from get_atomic_props
    secondMomentDict = {}
    #calculate the second moment based on the moments from dictionary
    #return dictionary
    secondMomentDict.update({'Qxx': atomDict['q'][0]*(atomDict['xyz'][0]**2) + 
                            (atomDict['Q_XX'][0] + atomDict['R+2'][0])/3 + atomDict['xyz'][0]*atomDict['Mu_Intra_X'][0] + 
                            atomDict['xyz'][0]*atomDict['Mu_Intra_X'][0]})
    secondMomentDict.update({'Qxy': atomDict['q'][0]*atomDict['xyz'][0]*atomDict['xyz'][1] + 
                            atomDict['Q_XY'][0]/3+ atomDict['xyz'][1]*atomDict['Mu_Intra_X'][0] + 
                            atomDict['xyz'][0]*atomDict['Mu_Intra_Y'][0]})
    secondMomentDict.update({'Qxz': atomDict['q'][0]*atomDict['xyz'][0]*atomDict['xyz'][2] + 
                            atomDict['Q_XZ'][0]/3+ atomDict['xyz'][2]*atomDict['Mu_Intra_X'][0] + 
                            atomDict['xyz'][0]*atomDict['Mu_Intra_Z'][0]})
    secondMomentDict.update({'Qyy': atomDict['q'][0]*(atomDict['xyz'][1]**2) + 
                            (atomDict['Q_YY'][0] + atomDict['R+2'][0])/3+ atomDict['xyz'][1]*atomDict['Mu_Intra_Y'][0] + 
                            atomDict['xyz'][1]*atomDict['Mu_Intra_Y'][0]})
    secondMomentDict.update({'Qyz': atomDict['q'][0]*atomDict['xyz'][1]*atomDict['xyz'][2] + 
                            atomDict['Q_YZ'][0]/3+ atomDict['xyz'][2]*atomDict['Mu_Intra_Y'][0] + 
                            atomDict['xyz'][1]*atomDict['Mu_Intra_Z'][0]})
    secondMomentDict.update({'Qzz': atomDict['q'][0]*(atomDict['xyz'][2]**2) + 
                            (atomDict['Q_ZZ'][0] + atomDict['R+2'][0])/3 + atomDict['xyz'][2]*atomDict['Mu_Intra_Z'][0] + 
                            atomDict['xyz'][2]*atomDict['Mu_Intra_Z'][0]})
    secondMomentDict.update({'trace': secondMomentDict['Qxx'] + secondMomentDict['Qyy']+secondMomentDict['Qzz']})
    atomicQuadrupoleDict = {}
    #get the atomic contribution to quadrupole
    atomicQuadrupoleDict.update({'Q_xx': [0.5*(3*secondMomentDict['Qxx']-secondMomentDict['trace'])]})
    atomicQuadrupoleDict.update({'Q_xy': [0.5*(3*secondMomentDict['Qxy'])]})
    atomicQuadrupoleDict.update({'Q_xz': [0.5*(3*secondMomentDict['Qxz'])]})
    atomicQuadrupoleDict.update({'Q_yy': [0.5*(3*secondMomentDict['Qyy']-secondMomentDict['trace'])]})
    atomicQuadrupoleDict.update({'Q_yz': [0.5*(3*secondMomentDict['Qyz'])]})
    atomicQuadrupoleDict.update({'Q_zz': [0.5*(3*secondMomentDict['Qzz']-secondMomentDict['trace'])]})
    return atomicQuadrupoleDict  #return dictionary

def trial_atomic_contrib_rotation(atomDict, rotMat):
    #trial atomic contrib with rotation. Rotate second moment matrix after construction from original coordinates.
    #did not seem to work
    secondMomentMat = np.array([[(atomDict['Q_XX'][0] + atomDict['R+2'][0])/3,atomDict['Q_XY'][0]/3, atomDict['Q_XZ'][0]/3],
                                [atomDict['Q_XY'][0]/3,(atomDict['Q_YY'][0] + atomDict['R+2'][0])/3,atomDict['Q_YZ'][0]/3],
                                [atomDict['Q_XZ'][0]/3,  atomDict['Q_YZ'][0]/3,  (atomDict['Q_ZZ'][0] + atomDict['R+2'][0])/3]])
    secondRotMat = np.matmul(rotMat,secondMomentMat)
    secondMomentDict={}
    secondMomentDict.update({'Qxx': atomDict['q'][0]*(atomDict['xyz'][0]**2) + 
                            secondRotMat[0,0] + atomDict['xyz'][0]*atomDict['Mu_Intra_X'][0] + 
                            atomDict['xyz'][0]*atomDict['Mu_Intra_X'][0]})
    secondMomentDict.update({'Qxy': atomDict['q'][0]*atomDict['xyz'][0]*atomDict['xyz'][1] + 
                           secondRotMat[0,1]+ atomDict['xyz'][1]*atomDict['Mu_Intra_X'][0] + 
                            atomDict['xyz'][0]*atomDict['Mu_Intra_Y'][0]})
    secondMomentDict.update({'Qxz': atomDict['q'][0]*atomDict['xyz'][0]*atomDict['xyz'][2] + 
                            secondRotMat[0,2]+ atomDict['xyz'][2]*atomDict['Mu_Intra_X'][0] + 
                            atomDict['xyz'][0]*atomDict['Mu_Intra_Z'][0]})
    secondMomentDict.update({'Qyy': atomDict['q'][0]*(atomDict['xyz'][1]**2) + 
                            secondRotMat[1,1]+ atomDict['xyz'][1]*atomDict['Mu_Intra_Y'][0] + 
                            atomDict['xyz'][1]*atomDict['Mu_Intra_Y'][0]})
    secondMomentDict.update({'Qyz': atomDict['q'][0]*atomDict['xyz'][1]*atomDict['xyz'][2] + 
                            secondRotMat[1,2]+ atomDict['xyz'][2]*atomDict['Mu_Intra_Y'][0] + 
                            atomDict['xyz'][1]*atomDict['Mu_Intra_Z'][0]})
    secondMomentDict.update({'Qzz': atomDict['q'][0]*(atomDict['xyz'][2]**2) + 
                            secondRotMat[2,2]+ atomDict['xyz'][2]*atomDict['Mu_Intra_Z'][0] + 
                            atomDict['xyz'][2]*atomDict['Mu_Intra_Z'][0]})
    secondMomentDict.update({'trace': secondMomentDict['Qxx'] + secondMomentDict['Qyy']+secondMomentDict['Qzz']})
    atomicQuadrupoleDict = {}
    #get the atomic contribution to quadrupole
    atomicQuadrupoleDict.update({'Q_xx': [0.5*(3*secondMomentDict['Qxx']-secondMomentDict['trace'])]})
    atomicQuadrupoleDict.update({'Q_xy': [0.5*(3*secondMomentDict['Qxy'])]})
    atomicQuadrupoleDict.update({'Q_xz': [0.5*(3*secondMomentDict['Qxz'])]})
    atomicQuadrupoleDict.update({'Q_yy': [0.5*(3*secondMomentDict['Qyy']-secondMomentDict['trace'])]})
    atomicQuadrupoleDict.update({'Q_yz': [0.5*(3*secondMomentDict['Qyz'])]})
    atomicQuadrupoleDict.update({'Q_zz': [0.5*(3*secondMomentDict['Qzz']-secondMomentDict['trace'])]})
    return atomicQuadrupoleDict  #return dictionary

def get_sub_props(atomDict,subAtoms,atomList):
    # atom dict - object from get_atomic_props, after updating with get_atomic_quad_contrib - ALL atom 
    #subAtoms - indices of atoms in group. Starting at 1
    #    atomList - list of atom labels for all atoms in dictionary 
    # returns dictionary of substituent properties
    groupDict={} #create empty dict
    first=True #flag for first iteration through loop
    for atom in subAtoms:
        if first: #if first time, create all elements in dicionary
            first=False #don't come here again
            #could do this with a for loop, but only want select keys, not all
            #add the atomic property to the dictionary - at this point just atomic, will update later
            groupDict = {'q': [atomDict[atomList[atom-1]]['q'][0]], 
                         'K': [atomDict[atomList[atom-1]]['K'][0]],
                         'K_Scaled': [atomDict[atomList[atom-1]]['K_Scaled'][0]],
                         'Mu_Intra_X': [atomDict[atomList[atom-1]]['Mu_Intra_X'][0]],
                         'Mu_Intra_Y': [atomDict[atomList[atom-1]]['Mu_Intra_Y'][0]],
                         'Mu_Intra_Z': [atomDict[atomList[atom-1]]['Mu_Intra_Z'][0]],
                         'Mu_Bond_X': [atomDict[atomList[atom-1]]['Mu_Bond_X'][0]],
                         'Mu_Bond_Y': [atomDict[atomList[atom-1]]['Mu_Bond_Y'][0]],
                         'Mu_Bond_Z': [atomDict[atomList[atom-1]]['Mu_Bond_Z'][0]],
                         'Mu_X': [atomDict[atomList[atom-1]]['Mu_X'][0]],
                         'Mu_Y': [atomDict[atomList[atom-1]]['Mu_Y'][0]],
                         'Mu_Z': [atomDict[atomList[atom-1]]['Mu_Z'][0]],
                         'Q_xx': [atomDict[atomList[atom-1]]['quadContrib']['Q_xx'][0]],
                         'Q_xy': [atomDict[atomList[atom-1]]['quadContrib']['Q_xy'][0]],
                         'Q_xz': [atomDict[atomList[atom-1]]['quadContrib']['Q_xz'][0]],
                         'Q_yy': [atomDict[atomList[atom-1]]['quadContrib']['Q_yy'][0]],
                         'Q_yz': [atomDict[atomList[atom-1]]['quadContrib']['Q_yz'][0]],
                         'Q_zz': [atomDict[atomList[atom-1]]['quadContrib']['Q_zz'][0]],
                         'Vol': [atomDict[atomList[atom-1]]['Vol'][0]]}
        else: #for the rest, add the atomic property to the group dictionary element
            for prop in groupDict:
                if 'Q' not in prop:
                    groupDict[prop][0] += atomDict[atomList[atom-1]][prop][0]
                else:
                    groupDict[prop][0] += atomDict[atomList[atom-1]]['quadContrib'][prop][0]
    groupDict.update({'|Mu_Intra|' : [math.sqrt(groupDict['Mu_Intra_X'][0]**2 + groupDict['Mu_Intra_Y'][0]**2 + groupDict['Mu_Intra_Z'][0]**2)]}) 
    groupDict.update({'|Mu_Bond|' : [math.sqrt(groupDict['Mu_Bond_X'][0]**2 + groupDict['Mu_Bond_Y'][0]**2 + groupDict['Mu_Bond_Z'][0]**2)]}) 
    groupDict.update({'|Mu_Bond|' : [math.sqrt(groupDict['Mu_X'][0]**2 + groupDict['Mu_Y'][0]**2 + groupDict['Mu_Z'][0]**2)]})     
    return groupDict

def extract_sub_props(sumFileNoExt,subAtoms,groupProps=True,bcpId = [[1,2]],lapRhoCpAtoms=[]):
    #returns a dictionary of all group properties - bcp, group, atomic, and vscc
    #sumFileNoExt - a sumfile name, without the .sum
    #subAtoms - indices of atoms in the molecule comprising group - starts at 1
    #groupProps - boolean - do you want to compute group properties
#    bcpIdx - list of 2 length lists. Pass empty list if no bcp properties wanted. 2 length lists are indices of atoms that you want BCPs for
    #lapRhoCpAtoms = list of atom indices that you want to find laprhoCPs for. Defaults to empty
    #vol surface is 0.001 au isodensity
    sumFile = open(sumFileNoExt+".sum","r") #open file, read lines, close file
    data = sumFile.readlines()
    sumFile.close()
    #find atom labels in atomList
    atomTable = get_table(data,'q(A)              L(A)              K(A)          K_Scaled(A)      |Mu_Intra(A)|')
    atomList = list(atomTable['Atom'])
    subatomLabels = []
    for i,atom in enumerate(atomList):
        if any(x == i+1 for x in subAtoms):
            subatomLabels.append(atomList[i])
    #if bcpIdx not empty, get bcp properties for each listed bcp
#right now every time I run, it passes the bcpIdx that it starts with to the next loop. Weird af
    bcpIdx = []
    #so here, code was doing something weird. bcpId was stored in memory and not reset to [[1,2]] every function call
        #every iteration updated to to -1 for each index, so I'm starting at [1,2] in bcpId and copying the values(
        #not pointers to bcpIdx)
        #don't know why it's necessary but it works now.
        #anyways that's why these 5 lines are here
    atomicProps = get_atomic_props(data) # get atomic dictionary
    for bcp in bcpId: 
        atList=[]
        for at in bcp:
            atList.append(at)
        bcpIdx.append(atList)    

    if len(bcpIdx) > 0:
        bcpProperties = {} #initialize empty dictionary
        for bcpPair in bcpIdx:
            bcpPair[0] -= 1 # lists start from 0, indices would be passed from molecule
            bcpPair[1] -= 1 # which would start at 1, so adjust index
            bcpBlock=[] #reset bcpBlock so won't be copied into future loop if block unable to be found
            bcpBlock = get_bcp_block(data, [atomList[i] for i in bcpPair])
            #add the bcp properties to the dictionary under the atom labels
            prop = get_bcp_properties(bcpBlock)
            prop.update({'DI(R,G)': [get_sub_di(data,subatomLabels)]})
            for key in atomicProps:
                if key == atomList[bcpPair[0]] or key == atomList[bcpPair[1]]:
                    if '1' in key:
                        keynum=1
                    elif '2' in key:
                        keynum=2
                    prop.update({'r(BCP-{at})'.format(at=keynum):[get_dist(atomicProps[key]['xyz'],prop['xyz'])]})
            bcpProperties.update({'{at1}-{at2}'.format(at1=atomList[bcpPair[0]],at2=atomList[bcpPair[1]]):prop })        
    else: #
        bcpProperties = "BCP Properties not requested"
    if groupProps: #if you want group properties
        for atom in atomicProps: #update atomic dictionary with quadrupole contribution
            atomicProps[atom].update({'quadContrib': get_atomic_quad_contrib(atomicProps[atom])})
        subDict = get_sub_props(atomicProps,subAtoms,atomList) # get substituent properties
    else:
        subDict = "Group Properties not requested"

    if len(lapRhoCpAtoms) > 0: #if we want laprhocps for at least one atom
        vsccProps={}
        for atom in lapRhoCpAtoms:#for each atom requested, get lapRhoCps
            allCC = get_cc_props(sumFileNoExt,atomList[atom-1])
            vsccProps.update({atomList[atom-1]: identify_vscc(allCC,atomicProps)})
    else:
        vsccProps={"VSCC Props not requested"}
    #create output dictionary to return all requested properties, if a property not requested return a string stating that    
    outDict = {'Group':subDict,'Atomic':atomicProps,'BCP':bcpProperties,'VSCC':vsccProps}    
    return outDict

def find_connected(data,negXAtomLabel,originAtomLabel):
    #find all atoms to which bonded
    bcpLines=[]
    for line in data:
        if 'Type = (3,-1)' in line and negXAtomLabel not in line and originAtomLabel in line:
            bcpLines.append(line)
    bcpList = []        
    for bcp in bcpLines:
        splitbcp = bcp.split()
        bcpList.append([splitbcp[4],splitbcp[5]])
    return bcpList


def sub_prop_frame(csvFile):
    csvFrame = pd.read_csv(csvFile)
   
    #Take csv file for format and extract properties:
    #Substituent, subAtoms, label1,label2...
    #Substituent: string of substituent
    #subAtoms: string of space separated substituent atoms eg '1 3 4'
    #label1 - contains .sum file with no extension
    #labeli could be e.g. a model chemistry, or substrate
    #that is the column "SubH" would have sum files for the substituents attached to H
    #SubC6H5 would have sum files for substituents attached to C6H5 etc
    all_label_dict = {} #initialize output dictionary. Key structure: label:{Group,BCP}
    
    ncolumns = csvFrame.shape[1]
    nrow = csvFrame.shape[0]
    subAtoms = []
    for sub in range(0,nrow): #get subAtoms from csv in list format
        subAtomsString = csvFrame.loc[sub]['subAtoms'].split()
        subAtomsInt = [eval(i) for i in subAtomsString]
        subAtoms.append(subAtomsInt)
    for col in range(2,ncolumns): #2 is the start of the label columns, so for all labels:
        ind_label_dict = {} # create empty dictionary for this label
        count=0 #For first iteration for each label will create the data frame
        for row_num,sub in enumerate(csvFrame['Substituent']): #iterate over substituents
            #extract sumfile name from table
            sumFile = csvFrame[csvFrame['Substituent']==sub][csvFrame.columns[col]].iloc[0]
            #get properties
            extracted_props = extract_sub_props(sumFile,subAtoms[row_num])
            # dont want to store array in data frame
            excludeKeys = ['xyz']
            #add a substituent label to data frame
            extracted_props['Group'].update({'Substituent': sub})
            if count ==0:
                count = 1 #don't come here again for this label (don't need to make data frame again)
                #create data frame
                groupFrame = pd.DataFrame.from_dict(extracted_props['Group'],orient='columns')
                #for each bcp properties gotten for - currently only do this for one
                #should work generally soon
                for bnum,bcp in enumerate(extracted_props['BCP']): #currently only use for one bcp query
                    tempbcpDict = {k: extracted_props['BCP'][bcp][k] for k in set(list(extracted_props['BCP'][bcp].keys()))-set(excludeKeys)}
                    tempbcpDict.update({'x': [extracted_props['BCP'][bcp]['xyz'][0]]})
                    tempbcpDict.update({'y': [extracted_props['BCP'][bcp]['xyz'][1]]})
                    tempbcpDict.update({'z': [extracted_props['BCP'][bcp]['xyz'][2]]})
                    tempbcpDict.update({'Substituent' : sub})
                    tempbcpDict.update({'BCP': bcp})
                    if bnum ==0: #create data frame for first bcp in list
                        bcpFrame = pd.DataFrame.from_dict(tempbcpDict,orient='columns')
                    else: #else add to data frame
                        bcpFrame = pd.concat([bcpFrame,pd.DataFrame(tempbcpDict)],ignore_index=True)            
            else: #add to data frame after first iteration
                groupFrame = pd.concat([groupFrame,pd.DataFrame(extracted_props['Group'])],ignore_index=True)
                for bcp in extracted_props['BCP']:
                    tempbcpDict = {k: extracted_props['BCP'][bcp][k] for k in set(list(extracted_props['BCP'][bcp].keys()))-set(excludeKeys)}
                    tempbcpDict.update({'x': [extracted_props['BCP'][bcp]['xyz'][0]]})
                    tempbcpDict.update({'y': [extracted_props['BCP'][bcp]['xyz'][1]]})
                    tempbcpDict.update({'z': [extracted_props['BCP'][bcp]['xyz'][2]]})
                    tempbcpDict.update({'Substituent' : sub})
                    tempbcpDict.update({'BCP': bcp})
                    bcpFrame = pd.concat([bcpFrame,pd.DataFrame(tempbcpDict)],ignore_index=True)            
#create output dictioanry of data frames {label1:{Group,BCP},label2:{Group,BCP},....}
        ind_label_dict.update({'Group': groupFrame}) 
        ind_label_dict.update({'BCP': bcpFrame})
        all_label_dict.update({csvFrame.columns[col]: ind_label_dict })
    return all_label_dict         

def get_xyz(sumfile):
    sumFile = open(sumfile+".sum","r") #open file, read lines, close file
    data = sumFile.readlines()
    sumFile.close()
    xyzTable = get_table(data,'Charge                X                  Y                  Z',endString='Some Atomic Properties:')
    xyzTable['X'] = pd.to_numeric(xyzTable['X'],downcast='float')
    xyzTable['Y'] = pd.to_numeric(xyzTable['Y'],downcast='float')
    xyzTable['Z'] = pd.to_numeric(xyzTable['Z'],downcast='float')
    return {'xyz':xyzTable[['X','Y','Z']].to_numpy(),'Atoms':xyzTable['Atom']}



