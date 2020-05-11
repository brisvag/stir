# Copyright 2019-2019 the stir authors. See copying.md for legal info.

"""
Collection of utilities
"""

from pathlib import Path


def clean_path(path):
    """
    cleans up paths and resolves ~ and symlinks
    """
    return Path(path).expanduser().resolve()


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
