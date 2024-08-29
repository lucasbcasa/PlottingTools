import numpy as np
import matplotlib.pyplot as plt
from matplotlib.pyplot import figure

# Get the moire unit vectors in component indexing
def getMoireVecsComp(filename):
    with open(filename) as file:
        # The vectors are written as " 55 -27", so we neglect the first entry of the split
        mVecs = np.array([[int(num) for num in line.split(' ')[1:]] for line in file])
    return mVecs

# Get the site list for the first cell
def getSiteList(filename):
    with open(filename) as file:
        # The sites are written as "0 0 0 0"
        siteList = np.array([[int(num) for num in line.split(' ')] for line in file])
    return siteList

# Span all unit cells and add the equivalent sites from the siteList
def spanCells(mVecs, siteList, supercellSize):
    (M, N) = supercellSize
    # Get all inequivalent lattice vectors
    mVecsAll = [ [m * mVecs[0] + n * mVecs[1], m * mVecs[2] + n * mVecs[3]] for m in range(M) for n in range(N)]
    siteListSupercell = [] #np.zeros((len(sizeList), len(mVecsAll)))
    # Loop over all sites and vectors
    for site in siteList:
        for vec in mVecsAll:
            # Get the layer and corresponding vector
            layer = site[0]
            lvec = np.array([0,0,*vec[layer]])
            equivSite = site + lvec
            # Append to the list
            siteListSupercell += [tuple(equivSite)]
    return siteListSupercell

# Creates dicts that translate between linear and component index notation
def getToLinearDict(componentFilename):
    with open(componentFilename) as file:
        toLinear = dict([ [tuple([int(num) for num in line.split(' ')[1:]]), i] for i, line in enumerate(file)])
    return toLinear

# Converts a whole list from component to linear index notation
def toLinearList(siteListComponent, toLinearDict):
    siteListLinear = []
    for site in siteListComponent:
        x = toLinearDict[site]
        siteListLinear += [x]
    siteListLinear = np.array(siteListLinear)
    return siteListLinear

# Builds supercell list of sites in linear indexing from the list of the first cell; to be used by PFM method
def getSiteListSupercell(mVecsCompFile, siteListFile, supercellSize, toLinearDict):
    # Get the moire unit vectors in component indexing
    mVecs = getMoireVecsComp(mVecsCompFile)
    # Get the list of sites in the first cell
    siteList = getSiteList(siteListFile)
    # Span all unit cells
    siteListSupercell = spanCells(mVecs, siteList, supercellSize)
    # Translate from component to linear indexing
    #print(len(siteListSupercell))
    siteListSupercellLinear = toLinearList(siteListSupercell, toLinearDict)
    # Add to final list
    #saveList()
    return siteListSupercellLinear



def isIn(a, b):
    truthArray =[]
    for el in b:
        if el in a:
            truthArray += [True]
        else:
            truthArray += [False]
    return np.array(truthArray)

def plotSelectedSites(selectedSites, allSites, indices, positions, layers, sublattices, sets, plotsfolder, **kwargs):

    fig, axes = plt.subplots(1, 4,
                           gridspec_kw = {'wspace':0, 'hspace':0.1}, 
                           subplot_kw={'aspect':'equal'})
        
    x = positions[:,0]
    y = positions[:,1]
    truthArray = isIn(selectedSites, allSites)
    c = np.where(truthArray, ['red'], ['blue'])

    top = np.nonzero(layers-1)
    bottom = np.nonzero(layers)
    
    # List of sites in sublattices A/B
    a = np.nonzero(sublattices-1)
    b = np.nonzero(sublattices)

    # List of sites separated by sublattice and layer
    top_a = np.intersect1d(top, a)
    top_b = np.intersect1d(top, b)

    bottom_a = np.intersect1d(bottom, a)
    bottom_b = np.intersect1d(bottom, b)

    x_top, x_bottom = x[top], x[bottom]
    y_top, y_bottom = y[top], y[bottom]
    c_top, c_bottom = c[top], c[bottom]
    
    x_top_a, x_bottom_a = x[top_a], x[bottom_a]
    y_top_a, y_bottom_a = y[top_a], y[bottom_a]
    c_top_a, c_bottom_a = c[top_a], c[bottom_a]
    
    x_top_b, x_bottom_b = x[top_b], x[bottom_b]
    y_top_b, y_bottom_b = y[top_b], y[bottom_b]
    c_top_b, c_bottom_b = c[top_b], c[bottom_b]
    
    transp = kwargs.get('transp', False)
    if transp:
        a_top, a_bottom = abs(c_top), abs(c_bottom)
    else:
        a_top, a_bottom = 1, 1
            
#     n_top=colors.TwoSlopeNorm(vmin=np.amin(c_top), vcenter=0., vmax=np.amax(c_top))
#     n_bottom=colors.TwoSlopeNorm(vmin=np.amin(c_bottom), vcenter=0., vmax=np.amax(c_bottom))
#     divnorm = colors.TwoSlopeNorm(vmin=np.amin(c), vcenter=0., vmax=np.amax(c))
    
    # Mirroring scale
#     n_top=colors.TwoSlopeNorm(vmin=-np.amax(abs(c_top)), vcenter=0., vmax=np.amax(abs(c_top)))
#     n_bottom=colors.TwoSlopeNorm(vmin=-np.amax(abs(c_bottom)), vcenter=0., vmax=np.amax(abs(c_bottom)))
#     divnorm = colors.TwoSlopeNorm(vmin=-np.amax(abs(c)), vcenter=0., vmax=np.amax(abs(c)))


    dotSize = kwargs.get('dotSize', 4)
    #cmap = kwargs.get('cmap', 'cool')
    cmap = kwargs.get('cmap', 'seismic')

    ax1 = plt.subplot(1,4,1)
    ax1.set_title('Top A')
    #leftPanel = ax1.scatter(x_top,y_top,c=c_top, s=dotSize, cmap=cmap)
    #leftPanel = ax1.scatter(x_top,y_top,c=c_top, s=dotSize * a_top/np.amax(a_top), cmap=cmap, norm=divnorm)
    #leftPanel = ax1.scatter(x_top_a,y_top_a,c=c_top_a, s=dotSize, cmap=cmap, norm=divnorm)
    leftPanel = ax1.scatter(x_top_a,y_top_a,c=c_top_a, s=dotSize)

    #plt.colorbar(leftPanel, ax=ax1)

    ax2 = plt.subplot(1,4,2)
    ax2.set_title('Bottom A')
    #rightPanel = ax2.scatter(x_bottom_a,y_bottom_a,c=c_bottom_a, s=dotSize, cmap=cmap, norm=divnorm)
    rightPanel = ax2.scatter(x_bottom_a,y_bottom_a,c=c_bottom_a, s=dotSize)
    #plt.colorbar(rightPanel, ax=ax2)
    
    ####
    ax3 = plt.subplot(1,4,3)
    ax3.set_title('Top B')
    #leftPanel = ax3.scatter(x_top,y_top,c=c_top, s=dotSize, cmap=cmap)
    #leftPanel_b = ax3.scatter(x_top_b,y_top_b,c=c_top_b, s=dotSize, cmap=cmap, norm=divnorm)
    leftPanel_b = ax3.scatter(x_top_b,y_top_b,c=c_top_b, s=dotSize)

   # plt.colorbar(leftPanel_b, ax=ax3)

    ax4 = plt.subplot(1,4,4)
    ax4.set_title('Bottom B')
    #rightPanel_b = ax4.scatter(x_bottom_b,y_bottom_b,c=c_bottom_b, s=dotSize, cmap=cmap, norm=divnorm)
    rightPanel_b = ax4.scatter(x_bottom_b,y_bottom_b,c=c_bottom_b, s=dotSize)
    #plt.colorbar(rightPanel_b, ax=ax4)
    ####
    
    
#     if vacPos is not None:
#         vacIndex = getVacancyIndex(vacPos, indices)
#         vacX, vacY = positions[vacIndex]
#         vacInfo = (vacPos[0], vacPos[1])
#         # Mark the vacancy with a black dot
#         axes_dict = {
#                     '(0, 0)': ax1,
#                     '(1, 0)': ax2,
#                     '(0, 1)': ax3, 
#                     '(1, 1)':ax4
#                     }
#         vacAx = axes_dict[str(vacInfo)]
#         vacAx.scatter(vacX, vacY,c='k',s=dotSize*5, marker='x') 
    
    plt.tight_layout()
    
    fig.set_size_inches((16 * kwargs.get('supercellSize', [1])[0], 12))
    
    #baseline_string = ''
    #if kwargs.get('baseline') is not None: baseline_string='; vs Baseline'
    #plt.suptitle(kwargs.get('suptitle', "Vacancy @ {} - LDOS{}".format(vacPos, baseline_string)))
    #plt.savefig(plotsfolder + kwargs.get('figname', 'LDOS-v@{}{}'.format(vacPos, baseline_string)), dpi=750)