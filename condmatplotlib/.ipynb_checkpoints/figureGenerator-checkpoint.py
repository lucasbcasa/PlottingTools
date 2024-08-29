import matplotlib
import tkinter as tk


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
    # figName = figDir.rsplit(sep='/', maxsplit=2)[-2]
    # figName = figName.split(sep='-', maxsplit=1)[1]
    # figPath = figDir + figName
    figPath = figDir.rsplit(sep='/', maxsplit=1)[0]
    return figPath

def saveFigToPath(fig, figPath, fmt='png', dpi=200):
    fig.savefig(figPath + '.' + fmt, facecolor='w', bbox_inches = 'tight', format=fmt, dpi=dpi)
    return figPath