"""reference_maps
Dictionaries containing maps used in substituent reorientation
"""
import numpy as np
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