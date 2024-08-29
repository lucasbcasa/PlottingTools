from importlib import reload

from . import constants

def reload_modules():
    reload(constants)