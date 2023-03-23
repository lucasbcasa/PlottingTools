from importlib import reload
# from types import SimpleNamespace
import numpy as np
import matplotlib.pyplot as plt
import sys
sys.path.append('./lib')

# Constants and Definitions
from constants import *
import constants

# Bands Tools
from plottingMethods import *
import plottingMethods as pm

# Figure Generation
from figureGenerator import *
import figureGenerator as fg

# File Management
from fileManager import *
import fileManager as fm

# Useful for debugging and testing new features
def reload_modules():
    reload(constants)
    reload(pm)
    reload(fg)
    reload(fm)