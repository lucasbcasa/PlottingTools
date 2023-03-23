def plotLDOSmap(LDOS, indices, positions, layers, sublattices, vacPos, plotsfolder, **kwargs):

    fig, axes = plt.subplots(2, 2,
                           gridspec_kw = {'wspace':0, 'hspace':0}, 
                           subplot_kw={'aspect':'equal'})
        
    x = positions[:,0]
    y = positions[:,1]
    c = LDOS

    top = np.nonzero(layers-1)
    bottom = np.nonzero(layers)
    
#     top_a = np.nonzero((layers - 1) * (sublattices-1))
#     top_b = np.nonzero((layers - 1) * (sublattices))
    
#     bottom_a = np.nonzero((layers) * (sublattices-1))
#     bottom_b = np.nonzero((layers) * (sublattices))
    
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
            
    n_top=colors.TwoSlopeNorm(vmin=np.amin(c_top), vcenter=0., vmax=np.amax(c_top))
    n_bottom=colors.TwoSlopeNorm(vmin=np.amin(c_bottom), vcenter=0., vmax=np.amax(c_bottom))
    divnorm = colors.TwoSlopeNorm(vmin=np.amin(c), vcenter=0., vmax=np.amax(c))


    dotSize = kwargs.get('dotSize', 4)
    #cmap = kwargs.get('cmap', 'cool')
    cmap = kwargs.get('cmap', 'seismic')

    ax1 = plt.subplot(2,2,1)
    ax1.set_title('Top A')
    #leftPanel = ax1.scatter(x_top,y_top,c=c_top, s=dotSize, cmap=cmap)
    #leftPanel = ax1.scatter(x_top,y_top,c=c_top, s=dotSize * a_top/np.amax(a_top), cmap=cmap, norm=divnorm)
    leftPanel = ax1.scatter(x_top_a,y_top_a,c=c_top_a, s=dotSize, cmap=cmap, norm=divnorm)

    if vacPos is not None:
        vacIndex = getVacancyIndex(vacPos, indices)
        vacX, vacY = positions[vacIndex]
        #print(np.sort(abs(x_top-vacX)**2 + abs(y_top-vacY)**2))
        ax1.scatter(vacX, vacY,c='k',s=dotSize) # Mark the vacancy with a black dot

    #plt.colorbar(leftPanel, ax=ax1)

    ax2 = plt.subplot(2,2,2)
    ax2.set_title('Bottom A')
    rightPanel = ax2.scatter(x_bottom_a,y_bottom_a,c=c_bottom_a, s=dotSize, cmap=cmap, norm=divnorm)
    #plt.colorbar(rightPanel, ax=ax2)
    
    ####
    ax3 = plt.subplot(2,2,3)
    ax3.set_title('Top B')
    #leftPanel = ax3.scatter(x_top,y_top,c=c_top, s=dotSize, cmap=cmap)
    leftPanel_b = ax3.scatter(x_top_b,y_top_b,c=c_top_b, s=dotSize, cmap=cmap, norm=divnorm)

    if vacPos is not None:
        vacIndex = getVacancyIndex(vacPos, indices)
        vacX, vacY = positions[vacIndex]
        #print(np.sort(abs(x_top-vacX)**2 + abs(y_top-vacY)**2))
        ax3.scatter(vacX, vacY,c='k',s=dotSize) # Mark the vacancy with a black dot

   # plt.colorbar(leftPanel_b, ax=ax3)

    ax4 = plt.subplot(2,2,4)
    ax4.set_title('Bottom B')
    rightPanel_b = ax4.scatter(x_bottom_b,y_bottom_b,c=c_bottom_b, s=dotSize, cmap=cmap, norm=divnorm)
    #plt.colorbar(rightPanel_b, ax=ax4)
    ####
    
    baseline_string = ''
    if kwargs.get('baseline') is not None: baseline_string='; vs Baseline'
    plt.suptitle(kwargs.get('suptitle', "Vacancy @ {} - LDOS{}".format(vacPos, baseline_string)))
    plt.savefig(plotsfolder + kwargs.get('figname', 'LDOS-v@{}{}'.format(vacPos, baseline_string)))
    
    fig.set_size_inches((12, 12))
    