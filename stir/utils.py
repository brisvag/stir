# Copyright 2019-2019 the stir authors. See copying.md for legal info.

"""
Collection of utilities
"""

from pathlib import Path
import argparse
import psutil


# functions used by the argument parser to check if input files are valid
def valid_str(param):
    p = clean_path(param)
    if p.suffix not in ('.gro', '.pdb',) or not p.is_file():
        raise argparse.ArgumentTypeError(f'File {p} is not a valid structure file.')
    return p


def valid_traj(param):
    p = clean_path(param)
    if p.suffix not in ('.xtc',) or not p.is_file():
        raise argparse.ArgumentTypeError(f'{p} is not a valid gromacs trajectory file.')
    return p


def valid_top(param):
    p = clean_path(param)
    if p.suffix not in ('.top', '.itp', '.tpr') or not p.is_file():
        raise argparse.ArgumentTypeError(f'{p} is not a valid tpr or topology file.')
    return p


def clean_path(path):
    """
    cleans up paths and resolves ~ and symlinks
    """
    return Path(path).expanduser().resolve()


def enough_ram(traj_list, skip):
    """
    check if there's enough memory to load the requested trjectories
    """
    freemem = psutil.virtual_memory().available
    traj_size = 0
    for traj in traj_list:
        traj_size += clean_path(traj).stat().st_size
    # TODO: add keepwater info when we can actually use it to load selectively
    # check if there's enough free memory: the number 5 is based on some testing
    if freemem < 5*(traj_size/skip):
        return False
    return True


def stir_help():
    """
    STIR functions:

    - garnish: draw bonds between atoms based on topology
    - nice: apply a nice visualization preset
    - supercell: show neighboring periodic/symmetric cells
    - cheese: take a picture or a video of the system

    Simply run `help stir` or just `stir` to print this help again.

    Other useful PyMOL commands:
    - intra_fit: align a trajectory to a specific selection and frame
    """
    print(stir_help.__doc__)
