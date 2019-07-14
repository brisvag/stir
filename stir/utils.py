# Copyright 2019-2019 the mtools authors. See copying.md for legal info.

"""
Collection of utilities
"""

import os
import argparse
import psutil


# used by the argument parser to check if input files are valid
def valid_str(param):
    _, ext = os.path.splitext(param)
    if ext.lower() not in ('.gro',) or not os.path.isfile(param):
        raise argparse.ArgumentTypeError(f'File {param} must be a valid gromacs structure file.')
    return param


# used by the argument parser to check if input files are valid
def valid_traj(param):
    _, ext = os.path.splitext(param)
    if ext.lower() not in ('.xtc',) or not os.path.isfile(param):
        raise argparse.ArgumentTypeError(f'File {param} must be a valid gromacs trajectory file.')
    return param


# used by the argument parser to check if input files are valid
def valid_top(param):
    _, ext = os.path.splitext(param)
    if ext.lower() not in ('.top', '.itp', '.tpr') or not os.path.isfile(param):
        raise argparse.ArgumentTypeError(f'File {param} must be a valid tpr or topology file.')
    return param


def clean_path(path_in):
    """
    cleans up paths and resolves ~ and symlinks
    """
    return os.path.realpath(os.path.expanduser(os.path.expandvars(path_in)))


def enough_ram(traj_list, skip):
    # check if there's enough memory to load the requested traj and warn the user if needed:
    freemem = psutil.virtual_memory().available
    traj_size = 0
    for traj in traj_list:
        traj_size += os.path.getsize(clean_path(traj))
    # TODO: add keepwater info when we can actually use it to load selectively
    # check if there's enough free memory: the number 5 is based on some testing
    if freemem < 5*(traj_size/skip):
        return False
    return True
