import matplotlib
import tkinter as tk

###########################
# File managing functions #
###########################

# This methods are responsible for saving figures
def saveFigToFile(fig, fmt='png', dpi=200):
    root = tk.Tk()
    root.withdraw() # Hide the main window.
    root.call('wm', 'attributes', '.', '-topmost', True) # Raise the root to the top of all windows.
    file = tk.filedialog.asksaveasfilename() # Create a file instead
    if file=='': raise RuntimeError('No file selected!')
    else: fig.savefig(file, facecolor='w', bbox_inches = 'tight', format=fmt, dpi=dpi)
    return file

def getFigDir(file):
    path = file.rsplit(sep='/', maxsplit=1)[0] + '/'
    figPath = path.replace('/data/', '/figures/')
    return figPath

def getFigPath(figDir):
    figPath = figDir.rsplit(sep='/', maxsplit=1)[0]
    return figPath

def saveFigToPath(fig, figPath, fmt='png', dpi=200):
    fig.savefig(figPath + '.' + fmt, facecolor='w', bbox_inches = 'tight', format=fmt, dpi=dpi)
    return figPath

###################
# Actual plotters #
###################

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import plottingMethods

# These methods start the heavy work of creating figures
# They create subplot mosaiscs with the structure you ask for
# They then call the appropriate plotting functions depending on 
# which data type you are plotting
def makePanels(mosaic, **kwargs):
    panel_font = kwargs.get('panel_font', {'size': 15})
    fig_kw = dict(figsize=kwargs.get('figsize', (12,6)), dpi = kwargs.get('dpi', 200), facecolor='white')
    fig, axes = plt.subplot_mosaic(mosaic, sharex=False, sharey=False, **fig_kw)
    for i,(key,ax) in enumerate(axes.items()):
        labelToTheRight=(kwargs.get('fMosaic')[key]['panel_kwargs']).get('labelToTheRight',False)
        if labelToTheRight: panel_label_ha = 'right'
        else: panel_label_ha = 'left'
        if kwargs.get('panel_label', True):
            ax.text(*kwargs.get('panel_label_pos', (1*labelToTheRight,0.995)), 
                ax.get_label(), 
                font=panel_font, 
                transform = ax.transAxes, 
                verticalalignment=kwargs.get('panel_label_va', 'top'), 
                horizontalalignment=panel_label_ha)
    return fig, axes

def makeFigure(mosaic, fMosaic, **kwargs):
    fig, axes = makePanels(mosaic, fMosaic=fMosaic, **kwargs)
    for panel, content in fMosaic.items():
        ax = axes[panel]
        panel_kwargs = content['panel_kwargs']
        print(panel, content['type'])
        if content['type'] == 'templateDATA':
            plottingMethods.makeTemplatePanel(ax, content['data'], **panel_kwargs)
        elif content['type'] == 'anothertypeofDATA':
            # Call different plotting function here for a different data type
            pass
        else: raise ValueError('Not expected data type.')
    # fig.tight_layout(w_pad=0)
    return fig, axes