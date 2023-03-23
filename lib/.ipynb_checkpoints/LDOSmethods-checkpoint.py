import numpy as np
import matplotlib.pyplot as plt
from matplotlib.pyplot import figure
from matplotlib import colors
import scipy.optimize as opt

# This function normalizes the energy-integrated LDOS to 1
def normalizeLDOS(LDOS, kpoints=None):
    if kpoints: return LDOS/kpoints
    return LDOS/len(LDOS) # If no information on kpoints is given, normalize by total number of states

# Get the LDOS from file; optionally shift and normalize it
def getLDOS(filename, normalize=False, baseline=None):
    with open(filename) as file:
        LDOS = np.array([float(line) for line in file])
    if normalize: LDOS = normalizeLDOS(LDOS)
    if baseline is not None: LDOS -= baseline
    return LDOS

def getLatticeIndices(indicesFile):
    with open(indicesFile) as file:
        indices = np.array([[int(val) for val in line.split()] for line in file])
    return indices

def getLatticePositions(positionsFile):
    with open(positionsFile) as file:
        positions = np.array([[float(val) for val in line.split()] for line in file])
    return positions

def getLatticeParts(indices):
    layers, sublattices = (indices[:,0], indices[:, 1])
    return layers, sublattices

def getLatticeVecs(vecsFile):
    with open(vecsFile) as file:
        mVecs = np.array([[float(val) for val in line.split()] for line in file])
    return mVecs

# Fitting functions
def oneOverR(x, a0, a1, a2):
    return a0 * ( np.cos(a1*x + a2) ) / x**2

def gaussian(x, a, sigma):
    return a * np.exp(-0.5 * x**2 / sigma**2)

def makeFits(functions, x, y, cutoff, ax):

    for label, function in functions.items():
        params = opt.curve_fit(function, x, y, p0=(1,0,0))

        x_array = np.linspace(*cutoff,100)
        y_array = function(x_array, *params[0])

        ax.plot(x_array, y_array, label=label, lw=5)
    ax.legend()

# Return the indices where the condition is satisfied
def getValidPoints(x, xmin, xmax):
    return np.argwhere( (x>=xmin) & (x<=xmax) )

# Apply the cuttofs
def cut(x,f,cutoff):
    validPoints = getValidPoints(x, *cutoff).ravel()
    return (x[validPoints], f[validPoints])

def getVacancyPosition():
    return (0,0,0,0)

def getVacancyIndex(vacPos, indices):
    if vacPos is None: return None
    vacIndex = np.where((indices == vacPos).all(axis=1))[0][0]
    return vacIndex

def minDist(positions, moireVecs, vacX, vacY):
    
    mindist = np.array([np.sqrt((x-vacX)**2 + (y-vacY)**2) for x,y in zip(positions[:,0],positions[:,1])])
    
    # Loop over neighboring cells
    for m in (-1,0,1):
        for n in (-1,0,1):
            # Get a multiple of the moire unit vectors
            vx, vy = m*moireVecs[0] + n*moireVecs[1]
            # Get the distance to the corresponding site at neighbor cell
            dist = np.array([np.sqrt((x-vacX + vx)**2 + (y-vacY + vy)**2) for x,y in zip(positions[:,0],positions[:,1])])
            # Find the sites for which this multiple of the moire vector minimizes the distance
            locs = np.where(dist < mindist)
            # Save that
            mindist[locs] = dist[locs]
    return mindist
            

def getDist(positions, indices, vacPos, moireVecs=None):
    
    if vacPos is None:
        vacPos = (0,0,0,0)
        vacIndex = getVacancyIndex(vacPos, indices)
    vacIndex = getVacancyIndex(vacPos, indices)
    vacX, vacY = positions[vacIndex]
    
    if moireVecs is not None:
        dist = minDist(positions, moireVecs, vacX, vacY)
    else:
        dist = np.array([np.sqrt((x-vacX)**2 + (y-vacY)**2) for x,y in zip(positions[:,0],positions[:,1])])
    
    return dist
    
def getAngle(positions, indices, vacPos, sym=1):
    
    if vacPos is None:
        vacPos = (0,0,0,0)
        vacIndex = getVacancyIndex(vacPos, indices)

    
    vacIndex = getVacancyIndex(vacPos, indices)
    vacX, vacY = positions[vacIndex]
    angles = np.array([np.angle((x-vacX) + 1j*(y-vacY)) for x,y in zip(positions[:,0],positions[:,1])])
    return np.mod(angles, 2*np.pi/sym)

def getSets(layers, sublattices):
    # List of sites in the top/bottom layer
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

    # Useful dictionaries
    top = {'a': top_a, 'b': top_b}
    bottom = {'a': bottom_a, 'b': bottom_b}
    sets = {'top': top, 'bottom': bottom}
    
    return sets

def plotLDOSmap(LDOS, indices, positions, layers, sublattices, sets, vacPos, plotsfolder, **kwargs):

    fig, axes = plt.subplots(1, 4,
                           gridspec_kw = {'wspace':0, 'hspace':0.1}, 
                           subplot_kw={'aspect':'equal'})
        
    x = positions[:,0]
    y = positions[:,1]
    c = LDOS

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
    n_top=colors.TwoSlopeNorm(vmin=-np.amax(abs(c_top)), vcenter=0., vmax=np.amax(abs(c_top)))
    n_bottom=colors.TwoSlopeNorm(vmin=-np.amax(abs(c_bottom)), vcenter=0., vmax=np.amax(abs(c_bottom)))
    divnorm = colors.TwoSlopeNorm(vmin=-np.amax(abs(c)), vcenter=0., vmax=np.amax(abs(c)))


    dotSize = kwargs.get('dotSize', 4)
    #cmap = kwargs.get('cmap', 'cool')
    cmap = kwargs.get('cmap', 'seismic')

    ax1 = plt.subplot(1,4,1)
    ax1.set_title('Top A')
    #leftPanel = ax1.scatter(x_top,y_top,c=c_top, s=dotSize, cmap=cmap)
    #leftPanel = ax1.scatter(x_top,y_top,c=c_top, s=dotSize * a_top/np.amax(a_top), cmap=cmap, norm=divnorm)
    leftPanel = ax1.scatter(x_top_a,y_top_a,c=c_top_a, s=dotSize, cmap=cmap, norm=divnorm)

    #plt.colorbar(leftPanel, ax=ax1)

    ax2 = plt.subplot(1,4,2)
    ax2.set_title('Bottom A')
    rightPanel = ax2.scatter(x_bottom_a,y_bottom_a,c=c_bottom_a, s=dotSize, cmap=cmap, norm=divnorm)
    #plt.colorbar(rightPanel, ax=ax2)
    
    ####
    ax3 = plt.subplot(1,4,3)
    ax3.set_title('Top B')
    #leftPanel = ax3.scatter(x_top,y_top,c=c_top, s=dotSize, cmap=cmap)
    leftPanel_b = ax3.scatter(x_top_b,y_top_b,c=c_top_b, s=dotSize, cmap=cmap, norm=divnorm)

   # plt.colorbar(leftPanel_b, ax=ax3)

    ax4 = plt.subplot(1,4,4)
    ax4.set_title('Bottom B')
    rightPanel_b = ax4.scatter(x_bottom_b,y_bottom_b,c=c_bottom_b, s=dotSize, cmap=cmap, norm=divnorm)
    #plt.colorbar(rightPanel_b, ax=ax4)
    ####
    
    
    if vacPos is not None:
        vacIndex = getVacancyIndex(vacPos, indices)
        vacX, vacY = positions[vacIndex]
        vacInfo = (vacPos[0], vacPos[1])
        # Mark the vacancy with a black dot
        axes_dict = {
                    '(0, 0)': ax1,
                    '(1, 0)': ax2,
                    '(0, 1)': ax3, 
                    '(1, 1)':ax4
                    }
        vacAx = axes_dict[str(vacInfo)]
        vacAx.scatter(vacX, vacY,c='k',s=dotSize*5, marker='x') 
    
    plt.tight_layout()
    
    fig.set_size_inches((16 * kwargs.get('supercellSize', [1])[0], 12))
    
    baseline_string = ''
    if kwargs.get('baseline') is not None: baseline_string='; vs Baseline'
    plt.suptitle(kwargs.get('suptitle', "Vacancy @ {} - LDOS{}".format(vacPos, baseline_string)))
    plt.savefig(plotsfolder + kwargs.get('figname', 'LDOS-v@{}{}'.format(vacPos, baseline_string)), dpi=750)
    
# Plot LDOS vs distance, for both layers (columns) and sublattices (rows)
def  plotLDOSradial(LDOS, dist, angles, sets, vacPos, plotsfolder, functions=None, **kwargs):
    #plt.figure(figsize=(12, 6), dpi=160);
    
    cutoff = kwargs.get('cutoff', (1, 25))
    cmap = kwargs.get('cmap', 'cool')
    dotSize = kwargs.get('dotSize', 4)
    lineCutAngle = kwargs.get('lineCutAngle', None)
    
    fig, axes = plt.subplots(nrows=2,
                             ncols=2,
                             sharex=True, 
                             sharey=True, 
                             squeeze=False,
                             figsize=(12, 6),
                             dpi=160)
        
    for i, layer in enumerate( ('top', 'bottom') ):
        for j, sublattice in enumerate( ('a', 'b') ):
            
            ax = axes[j, i]

            if sublattice=='a': # Top subplots: sublattice A; make title
                ax.set_title(layer)
            if sublattice=='b': # Top subplots: sublattice A; make xlabel
                ax.set_xlabel("Distance from vacancy")
            if layer=='top': # Left subplots: top layer ; make ylabel
                ax.set_ylabel("LDOS - subl. {}".format(sublattice))
            
            sites = sets[layer][sublattice]
            
            x_all = dist[sites]
            y_all = LDOS[sites]
            c_all = angles[sites]
            
            if lineCutAngle is not None:
                epsilon = kwargs.get('epsilon',0.01)
                angleCutoff = (lineCutAngle-epsilon, lineCutAngle+epsilon)
                _, x_all = cut(c_all, x_all, angleCutoff)
                c_all, y_all = cut(c_all, y_all, angleCutoff)


            x_cut, y_cut = cut(x_all, y_all, cutoff)
            x_cut, c_cut = cut(x_all, c_all, cutoff)
            
            ax.scatter(x_cut, y_cut, c=c_cut, cmap=cmap, s=dotSize)
            ax.tick_params(labelleft=False)
            
            if functions: makeFits(functions, x_cut, y_cut, cutoff, ax)
    baseline_string = ''
    if kwargs.get('baseline') is not None: baseline_string='; vs Baseline'
    plt.suptitle(kwargs.get('suptitle', "Vacancy @ {} - LDOS{}".format(vacPos, baseline_string)))
    plt.savefig(plotsfolder + kwargs.get('figname', 'RadialLDOS-v@{}{}'.format(vacPos, baseline_string)))

def minDist_AA(x, y, moireVecs, AA_x, AA_y):
    
    mindist = np.sqrt((x-AA_x)**2 + (y-AA_y)**2)
    
    # Loop over neighboring cells
    for m in (-1,0,1):
        for n in (-1,0,1):
            # Get a multiple of the moire unit vectors
            vx, vy = m*moireVecs[0] + n*moireVecs[1]
            # Get the distance to the corresponding site at neighbor cell
            dist = np.sqrt((x-AA_x + vx)**2 + (y-AA_y + vy)**2)
            # Find the sites for which this multiple of the moire vector minimizes the distance
            if dist<mindist:
                mindist = dist
    return mindist

    
def getLDOS_AA(x_all, y_all, LDOS_all, uVecs, supercell, radius):
    # uVecs = unit cell vectors
    (u1, u2) = uVecs
    # supercell = size of supercell
    (m,n) = supercell
    # mVec = moire vectors
    mVec = (m * u1, n * u2)
    
    AA_pos = np.zeros((m,n,2))
    LDOS_AA = np.zeros((m,n))
    
    for i in range(m):
        for j in range(n):
            AA_pos[i,j] = (u1 * i + u2 * j)
            for (x, y, LDOS) in zip(x_all, y_all, LDOS_all):
                dist_to_AA = minDist_AA(x, y, mVec, *AA_pos[i,j])
                #if ( (x-AA_pos[i,j][0])**2 + (y-AA_pos[i,j][1])**2 <= radius**2): # If site near AA region, count its LDOS
                if ( dist_to_AA <= radius): # If site near AA region, count its LDOS
                    LDOS_AA[i,j] += LDOS
            
    return (AA_pos, LDOS_AA)
    
    
def plotLDOSradialAA(LDOS, positions, layers, sublattices, indices, uVecs, supercell, radius, dist, angles, sets, vacPos, plotsfolder, functions=None, **kwargs):
    #plt.figure(figsize=(12, 6), dpi=160);
    
    cutoff = kwargs.get('cutoff', (1, 25))
    cmap = kwargs.get('cmap', 'cool')
    dotSize = kwargs.get('dotSize', 4)
    lineCutAngle = kwargs.get('lineCutAngle', None)
    
    fig, axes = plt.subplots(nrows=2,
                             ncols=2,
                             sharex=True, 
                             sharey=True, 
                             squeeze=False,
                             figsize=(12, 6),
                             dpi=160)
        
    x = positions[:,0]
    y = positions[:,1]
    c = LDOS

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
    
    xs = (x_top_a, x_bottom_a, x_top_b, x_bottom_b)
    ys = (y_top_a, y_bottom_a, y_top_b, y_bottom_b)
    cs = (c_top_a, c_bottom_a, c_top_b, c_bottom_b)
        
        
    for i, layer in enumerate( ('top', 'bottom') ):
        for j, sublattice in enumerate( ('a', 'b') ):
            
            ax = axes[j, i]

            if sublattice=='a': # Top subplots: sublattice A; make title
                ax.set_title(layer)
            if sublattice=='b': # Top subplots: sublattice A; make xlabel
                ax.set_xlabel("Distance from vacancy")
            if layer=='top': # Left subplots: top layer ; make ylabel
                ax.set_ylabel("LDOS - subl. {}".format(sublattice))
            
            sites = sets[layer][sublattice]
            
            x_all = dist[sites]
            y_all = LDOS[sites]
            c_all = angles[sites]
            
            if lineCutAngle is not None:
                epsilon = kwargs.get('epsilon',0.01)
                angleCutoff = (lineCutAngle-epsilon, lineCutAngle+epsilon)
                _, x_all = cut(c_all, x_all, angleCutoff)
                c_all, y_all = cut(c_all, y_all, angleCutoff)


            x_cut, y_cut = cut(x_all, y_all, cutoff)
            x_cut, c_cut = cut(x_all, c_all, cutoff)
            
            x = xs[2*j + i]
            y = ys[2*j + i]
            c = cs[2*j + i]
            
            (AA_pos, AA_LDOS) = getLDOS_AA(x, y, c, uVecs, supercell, radius)
            
            vacIndex = getVacancyIndex(vacPos, indices)
            #print(vacIndex, positions[vacIndex])
            vacX, vacY = positions[vacIndex]
            print(AA_pos - (vacX, vacY))
            AA_rel_pos = AA_pos - (vacX, vacY)
            
            #AA_dist = np.sqrt(AA_rel_pos[:,:,0]**2 + AA_rel_pos[:,:,1]**2)
            
            (u1, u2) = uVecs
            # supercell = size of supercell
            (m,n) = supercell
            # mVec = moire vectors
            mVec = (m * u1, n * u2)
            
            AA_dist = [[minDist_AA(vacX, vacY, mVec, AA_pos[k,l,0], AA_pos[k,l,1]) for l in range(len(AA_pos[0]))] for k in range(len(AA_pos))]
            #AA_dist = [np.sqrt(a**2 + b**2) for (a,b) in AA_rel_pos]
            
            ax.scatter(AA_dist, AA_LDOS, marker='^', s=dotSize)
            ax.hlines(0, 0, np.amax(AA_dist), color='r')
            ax.tick_params(labelleft=False)
            
            if functions: makeFits(functions, x_cut, y_cut, cutoff, ax)
    baseline_string = ''
    if kwargs.get('baseline') is not None: baseline_string='; vs Baseline'
    plt.suptitle(kwargs.get('suptitle', "Vacancy @ {} - LDOS{}".format(vacPos, baseline_string)))
    plt.savefig(plotsfolder + kwargs.get('figname', 'AA LDOS-v@{}{}'.format(vacPos, baseline_string)))
    
def plotLDOSradialAA_singlePlot(LDOS, positions, layers, sublattices, indices, uVecs, supercell, radius, dist, angles, sets, vacPos, plotsfolder, functions=None, **kwargs):
    #plt.figure(figsize=(12, 6), dpi=160);
    
    cutoff = kwargs.get('cutoff', (1, 25))
    cmap = kwargs.get('cmap', 'cool')
    dotSize = kwargs.get('dotSize', 4)
    lineCutAngle = kwargs.get('lineCutAngle', None)
    
    fig, axes = plt.subplots(nrows=1,
                             ncols=1,
                             sharex=True, 
                             sharey=True, 
                             squeeze=False,
                             figsize=(6, 4),
                             dpi=160)
        
    x = positions[:,0]
    y = positions[:,1]
    c = LDOS

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
    
    xs = (x_top_a, x_bottom_a, x_top_b, x_bottom_b)
    ys = (y_top_a, y_bottom_a, y_top_b, y_bottom_b)
    cs = (c_top_a, c_bottom_a, c_top_b, c_bottom_b)
        
    arts = []
    minimum = 0
    for i, layer in enumerate( ('top', 'bottom') ):
        for j, sublattice in enumerate( ('a', 'b') ):
            
            #ax = axes[j, i]
            ax = axes[0,0]

            ax.set_title(kwargs.get('title',r"$\Delta$LDOS in AA regions"))
            ax.set_xlabel("Distance from vacancy")
            ax.set_ylabel(r"$\Delta$LDOS")
            
            sites = sets[layer][sublattice]
            
            x_all = dist[sites]
            y_all = LDOS[sites]
            c_all = angles[sites]
            
            if lineCutAngle is not None:
                epsilon = kwargs.get('epsilon',0.01)
                angleCutoff = (lineCutAngle-epsilon, lineCutAngle+epsilon)
                _, x_all = cut(c_all, x_all, angleCutoff)
                c_all, y_all = cut(c_all, y_all, angleCutoff)


            x_cut, y_cut = cut(x_all, y_all, cutoff)
            x_cut, c_cut = cut(x_all, c_all, cutoff)
            
            x = xs[2*j + i]
            y = ys[2*j + i]
            c = cs[2*j + i]
            
            (AA_pos, AA_LDOS) = getLDOS_AA(x, y, c, uVecs, supercell, radius)
            
            vacIndex = getVacancyIndex(vacPos, indices)
            vacX, vacY = positions[vacIndex]
            AA_rel_pos = AA_pos - (vacX, vacY)
            
            #AA_dist = np.sqrt(AA_rel_pos[:,:,0]**2 + AA_rel_pos[:,:,1]**2)
            
            (u1, u2) = uVecs
            # supercell = size of supercell
            (m,n) = supercell
            # mVec = moire vectors
            mVec = (m * u1, n * u2)
            
            AA_dist = [[minDist_AA(vacX, vacY, mVec, AA_pos[k,l,0], AA_pos[k,l,1]) for l in range(len(AA_pos[0]))] for k in range(len(AA_pos))]
            #AA_dist = [np.sqrt(a**2 + b**2) for (a,b) in AA_rel_pos]
            
            pts = ax.scatter(AA_dist, AA_LDOS, marker=kwargs.get('marker', 'D'), s=dotSize)
            #arts += pts
            pts.set_label(layer + ' ' + sublattice)
            
            ax.set_xlim(left = 0)
            xlim = ax.get_xlim()
            
            ax.hlines(0, *xlim, color='k', lw=4, ls='dashdot')
            #ax.tick_params(labelleft=False)
            
            minimum = min (np.amin(AA_LDOS), minimum)
            
            
            if functions: makeFits(functions, x_cut, y_cut, cutoff, ax)
    #ax.legend(arts,['Top A','Top B','Bottom A','Bottom B'])
    
    ax.set_yticks([minimum,0])
    ax.legend(loc='lower left')
    baseline_string = ''
    if kwargs.get('baseline') is not None: baseline_string='; vs Baseline'
    #plt.suptitle(kwargs.get('suptitle', "Vacancy @ {} - LDOS{}".format(vacPos, baseline_string)))
    plt.savefig(plotsfolder + kwargs.get('figname', 'AA LDOS-v@{}{}'.format(vacPos, baseline_string)))