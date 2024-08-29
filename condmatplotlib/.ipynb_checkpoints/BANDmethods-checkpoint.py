import numpy as np
import matplotlib.pyplot as plt
from matplotlib import animation, rc, patches
import constants
import math
import re

def getBAND(filename):
    with open(filename) as f:
        BAND = [line.split() for line in f] # Split the file into lines (strings), lines into lists
        BAND = [[float(item) for item in line] for line in BAND] # Flatten the list of lists
        BAND = np.transpose(np.asarray(BAND))
    return BAND

def prepareBANDList(files):
    BANDList = []
    for filename in files:
        BANDList.append(getBAND(filename))
    for i, BAND in enumerate(BANDList):
        BANDList[i] = BAND
    return BANDList

# def getKspace(file):
#     kspace = None
#     return kspace

def plotBAND(BANDList, **kwargs):
    plt.figure(figsize=kwargs.get('figsize',(10,10)), dpi=kwargs.get('dpi', 200), facecolor='w')
    
    plt.suptitle(kwargs.get('suptitle', r'Suptitle'))

    fontdict = kwargs.get('fontdict', dict(fontsize=kwargs.get('fontsize', 20), titlefontisze=15, labelfontdict=15))
    fontsize=fontdict['fontsize']
    titlefontsize=fontdict['titlefontsize']
    labelfontsize=fontdict['labelfontsize']
    
    ylim = kwargs.get('ylim')
    titles = kwargs.get('titles')
    
    hsp = kwargs.get('hsp', constants.HSP)
    indices, labels = hsp['indices'],  hsp['labels']
    
    ms=kwargs.get('markerSize',1)
    ylabel = kwargs.get('ylabel', 'Energy [t]')
    moireEnergies=kwargs.get('moireEnergies', None)
    yticks = kwargs.get('yticks',[0.004, 0.006])
    ytick_labels = kwargs.get('ytick_labels', [0.004, 0.006])
    panel_labels = kwargs.get('panel_labels', False)
    if kwargs.get('draw_lines', False): ls = 'solid'
    else: ls = 'none'
    nrows=kwargs.get('nrows', 1)
    ncols = -(len(BANDList) // -nrows) # Inverse floor division = ceiling division
    
    def isLeftPanel(index):
        return index%ncols==0

    def isBottomPanel(index):
        return index//ncols==nrows-1

    for i, BAND in enumerate(BANDList):

        ax = plt.subplot(nrows,ncols, i+1)
        kspace = kwargs.get('kspace', np.linspace(1, len(BAND), len(BAND)))
        xlim = kwargs.get('xlim', (kspace[0],kspace[-1]))
        ax.plot(kspace, BAND, linestyle=ls, marker='o', c='k', ms=ms)
        #ax.set_xticks([])
        if xlim: ax.set_xlim(*xlim)
        if ylim: ax.set_ylim(*ylim)
        if titles: ax.set_title(titles[i], fontsize=titlefontsize)

        for j in indices:
            ax.vlines(kspace[j], *ax.get_ylim(), color='k')
        if isBottomPanel(i):
            ax.set_xticks(kspace[indices])
            # ax.set_xticklabels(labels,fontsize=hsp.get('fontsize',8))
            ax.set_xticklabels(labels,fontsize=labelfontsize)
        else: 
            ax.set_xticks(())
        if not isLeftPanel(i):
            ax.set_yticks(())
        if isLeftPanel(i):
            ax.set_yticks(yticks, fontsize=fontsize)
            ax.set_yticklabels(ytick_labels,fontsize=labelfontsize)
            ax.set_ylabel(ylabel, fontsize=labelfontsize)
        if moireEnergies: 
            for h in moireEnergies:
                ax.hlines(h, *ax.get_xlim(), color='r')
                
        if panel_labels: 
            # ax.text(xlim[0], 0.95* ylim[-1], constants.panel_labels[i], fontsize=fontsize)
            ax.text(0,1, constants.panel_labels[i], fontsize=fontsize, transform = ax.transAxes, verticalalignment='top')
    plt.tight_layout()

# Takes string, delimiters and returns the first occurence
def getValue(string, delimiters, **kwargs):
    leftDelimiter = re.compile(delimiters[0])
    rightDelimiter = re.compile(delimiters[1])
    string = leftDelimiter.split(string,1)
    if len(string) == 1: return kwargs.get('default')
    string = rightDelimiter.split(string[1],1)[0]
    return float(string)

def prepareBANDDict(files, **kwargs):
    delimiters = kwargs.get('delimiters', (r'_en', r'_'))
    default = kwargs.get('default', math.inf)
        
    # dict(potential:file)
    potentialDict = {getValue(file, delimiters, default=default): file for file in files}
    potentials = list(potentialDict.keys())
    potentials.sort(reverse=kwargs.get('reverse', False))

    potentialDict = {value: potentialDict[value] for value in potentials}
    
    # dict(potential:band)
    BANDDict = {potential: getBAND(filename) for (potential, filename) in potentialDict.items()}
    
    return BANDDict

def processBANDDict(BANDDict):
    parameters = np.asarray(list(BANDDict.keys()))
    data = np.asarray(list(BANDDict.values()))
    kspace = len(data[0,0]) * list(np.linspace(1, len(data[0]), len(data[0])))

    data = np.transpose(data,axes=(0,2,1))
    
    return data, kspace, parameters

def getAnimation(data, x, **kwargs):
    fig, ax = plt.subplots(figsize=(8,6), dpi=200)
    if kwargs.get('title'): plt.title(kwargs.get('title'))
    legend = kwargs.get('legend')
    ax.set_xlim(kwargs.get('xlim', (x[0],x[-1])))
    ax.set_ylim(kwargs.get('ylim', (0,0.01)))
    ax.set_xticks([])
    
    vec = np.linspace(1, 574, 574)
    hsp = kwargs.get('hsp')
    if hsp:
        indices = hsp['indices']
        labels = hsp['labels']
        for i in indices:
            ax.vlines(vec[i], *ax.get_ylim(), color='k')
        ax.set_xticks(vec[indices])
        ax.set_xticklabels(labels,fontsize=hsp.get('fontsize',8))
    

#     rect =  patches.Rectangle((0.8,0.8), 0.1, 0.1, label='1')
#     ax.add_patch(rect)
    
    line, = ax.plot([], [], linestyle='none', marker='o', c='k', ms=2)
    
    if legend: text = ax.text(0, ax.get_ylim(), legend.format(0), c='r', backgroundcolor='w')

    def init():
        line.set_data([], []);
        return (line,)

    def getAnimator(data, x):
        return lambda i: animate(data, x, i)
    
    def animate(data, x, i):
        y = data[i]
        line.set_data(x, y)
        if legend: text.set_text(legend.format(parameters[i]))
        return (line,)
    
    animator = getAnimator(data, x)

    anim = animation.FuncAnimation(fig, animator, init_func=init,
                                   frames=len(data), interval=kwargs.get('interval',200), 
                                   blit=True);
    return anim