"""
Collection of utilities
"""

import os
import argparse


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
