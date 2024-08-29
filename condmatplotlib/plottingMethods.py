import numpy as np
import matplotlib.pyplot as plt

############
# Plotting #
############

# Creates a single plot in the subplot mosaic.
# If you have differnt types of data you want to have in your project
# (like band structures, LDOS maps, etc)
# You should define different functions like this and then call them 
# as needed from the figureGenerator script
def makeTemplatePanel(ax, data, **kwargs):
    fontdict = kwargs.get('fontdict', dict(fontsize=kwargs.get('fontsize', 20), titlefontsize=15, labelfontsize=15, tickfontsize=10))
    fontsize=fontdict['fontsize']
    titlefontsize=fontdict['titlefontsize']
    labelfontsize=fontdict['labelfontsize']
    ticksfontsize=fontdict['tickfontsize']

    ax.tick_params(labelsize=ticksfontsize)

    plotPanel(ax, data, **kwargs)

    # Here can go a lot of different customization of the plot layout
    xlim = kwargs.get('xlim', ax.get_xlim())
    ax.set_xlim(*xlim)
    xlabel = kwargs.get('xlabel', r'X-Label')
    ax.set_xlabel(xlabel, fontsize=labelfontsize, labelpad=0)
    xticks = kwargs.get('xticks', ax.get_xticks())
    ax.set_xticks(xticks)
    xtick_labels = kwargs.get('xtick_labels', ax.get_xticklabels())
    ax.set_xticks(xticks, labels=xtick_labels)
    # if kwargs.get('xtick_labels', False):
    #     xtick_labels = kwargs.get('xtick_labels', ax.get_xticklabels())

    ylim = kwargs.get('ylim', ax.get_ylim())
    ax.set_ylim(*ylim)
    ylabel = kwargs.get('ylabel', r'Y-Label')
    ax.set_ylabel(ylabel, fontsize=labelfontsize, labelpad=0)
    yticks = kwargs.get('yticks', ax.get_yticks())
    ax.set_yticks(yticks)
    ax.set_xticks(xticks, labels=xtick_labels)
    
    if kwargs.get('xtick_labels', False):
        ytick_labels = kwargs.get('ytick_labels', ax.get_yticklabels())
        ax.set_yticklabels(ytick_labels)
    ax.yaxis.set_label_coords(-0.015, 0.6)
    
    title = kwargs.get('title')
    if title: ax.set_title(title, fontsize=titlefontsize)

# The actual data plotting comes here
def plotPanel(ax, data, **kwargs):
    if kwargs.get('draw_lines', False): ls = 'solid'
    else: ls = 'none'
    ax.plot(data['x'], data['y'], linestyle=ls, marker='o', c='k', ms=kwargs.get('markerSize',1))