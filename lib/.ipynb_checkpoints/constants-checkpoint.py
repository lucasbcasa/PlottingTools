kB = 8.617333262145E-5 # [eV / K]
t = 2.7 # [eV]
ratio = 2.7 / kB # [K]

tablesPATH = './tables/'
indicesFile = tablesPATH + 'PRINT/55_sites_list'
positionsFile = tablesPATH + 'PRINT/55_sites_positions'
vecsFile = tablesPATH + 'PRINT/55_moire_vectors'
kspaceFile = tablesPATH + 'kp_list'

# High Symmetry Points
HSP = {'indices':[0,243-1,454-1,574-1], 'labels':[r'$K$',r'$\Gamma$',r'$M$',r'$K$'], 'fontsize':16}

import string
panel_labels = ['('+letter+')' for letter in string.ascii_lowercase]
