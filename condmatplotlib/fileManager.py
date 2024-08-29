# File importing and saving
import tkinter as tk
import tkinter.filedialog as fd
import numpy as np

# This method opens a new prompt window where you can select files from which to import your data
def selectFiles(**kwargs): # Currently the 'multiple' functionality breaks the 'dialog' functionality
    # If multiple files are wanted, call the selectFiles method for single file until no file is selected
    if kwargs.get('multiple'): 
        files=[]
        while(True):
            file = selectFiles() # Returns list of str with possibly several file paths
            if file == '': break # If no file is selected exit loop
            files += file # If one or more files are selected, add them all to the list and repeat
    # If the selectFiles method is selected for without the multiple condition just return the file path
    else:
        root = tk.Tk()
        root.withdraw() # Hide the main window.
        root.call('wm', 'attributes', '.', '-topmost', True) # Raise the root to the top of all windows.
        files = fd.askopenfilenames(parent=root, title=kwargs.get('dialog', None))
    return files

# This method creates a list with each element representing a dataset
def prepareDATAList(files):
    DATAList = []
    for filename in files:
        DATAList.append(getDATA(filename))
    for i, DATA in enumerate(DATAList):
        DATAList[i] = DATA
    return DATAList

# This method processes each file to import the data.
# Modify this in order to adapt to your data format (.npy, hdf5, .txt, etc)
# Currently it assumes a simple .txt file, with x-data in the top row and y-data in the bottom row
def getDATA(filename):
    with open(filename) as f:
        DATAset = [line.split() for line in f] # Split the file into lines (strings), lines into lists
        DATAset = [[float(item) for item in line] for line in DATAset] # Flatten the list of lists
    DATA = dict(x=DATAset[0], y=DATAset[1])
    return DATA

