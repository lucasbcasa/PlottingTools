import numpy as np
import seaborn as sb
import matplotlib.pyplot as plt
import constants

def getDOS(filename):
    with open(filename) as f:
        DOS = [line.split() for line in f] # Split the file into lines (strings), lines into lists
        print(type(DOS))
        print(len(DOS))
        DOS = [float(item) for line in DOS for item in line] # Flatten the list of lists
    return np.asarray(DOS)

def filterDOS(DOS, xmin, xmax):
    return DOS[(DOS>=xmin) & (DOS<=xmax)] # Filter the DOS for values in a given range
    
def printBandNumber(DOS, kpoints):
    print(len(DOS)/kpoints)
    
def eVconvert(DOS):
    return DOS * constants.t

def subtractBaseline(DOS, baseline): # not working
    nbins = 100
    hmin = np.minimum(DOS, baseline)
    hmax = np.maximum(DOS, baseline)
    
    #histogram is build with fixed min and max values
    hist1, _ = np.histogram(DOS,range=(hmin,hmax), bins=nbins)
    hist2, _ = np.histogram(baseline,range=(hmin,hmax), bins=nbins)

    #makes sense to have only positive values 
    diff = np.clip(hist1 - hist2, 0, None)
    return diff

def prepareDOSList(files, eV=False, filterRange=None, baseline=False):
    DOSList = []
    for filename in files:
        DOSList.append(getDOS(filename))
    for i, DOS in enumerate(DOSList):
        #if baseline: DOS = DOS
        if filterRange: DOS = filterDOS(DOS, *filterRange)
        if eV: DOS = eVconvert(DOS)
        #DOSList[i] = np.histogram(DOS)
        DOSList[i] = DOS
    return DOSList
    
def plotDOS(DOSList, **kwargs):
    plt.figure(figsize=(16, 8), dpi=80)

    bw = kwargs.get('bws', 0.05)
    kde_kws = kwargs.get('kde_kws', {'bw':bw})
    nbin = kwargs.get('nbins', 100)
    xlim = kwargs.get('xlim')
    ylim = kwargs.get('ylim')
    yticks = kwargs.get('yticks')
    xticks = kwargs.get('xticks')
    titles = kwargs.get('titles')
    color = kwargs.get('color', 'b')
    rug = kwargs.get('rug', False)


    for i, DOS in enumerate(DOSList):

        plt.subplot(1,len(DOSList), i+1)

        sb.distplot(DOS, 
                    bins= nbin, 
                    kde_kws=kde_kws, 
                    color=color, 
                    rug=rug)
        if xlim: plt.xlim(*xlim)
        if ylim: plt.ylim(*ylim)
        if xticks: plt.xticks(xticks)
        if yticks: plt.yticks(yticks)

        if titles: plt.title(titles[i], fontsize=kwargs.get('fontsize', 10))
        plt.xlabel("E [t]")
        if i==0: plt.ylabel('DOS')
        else: plt.ylabel('')