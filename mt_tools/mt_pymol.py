#!/usr/bin/env python

"""
Python wrapper for pymol to easily visualize martini trajectories
"""

import os
import argparse
import pymol
from pymol import cmd
import __main__


def valid_str(param):
    _, ext = os.path.splitext(param)
    if ext.lower() not in ('.gro',) or not os.path.isfile(param):
        raise argparse.ArgumentTypeError(f'File {param} must be a valid gromacs structure file.')
    return param


def valid_traj(param):
    _, ext = os.path.splitext(param)
    if ext.lower() not in ('.xtc',) or not os.path.isfile(param):
        raise argparse.ArgumentTypeError(f'File {param} must be a valid gromacs trajectory file.')
    return param


def valid_top(param):
    _, ext = os.path.splitext(param)
    if ext.lower() not in ('.top', '.itp', '.tpr') or not os.path.isfile(param):
        raise argparse.ArgumentTypeError(f'File {param} must be a valid tpr or topology file.')
    return param


parser = argparse.ArgumentParser(prog='mt_pymol')

parser.add_argument(dest='S', type=valid_str,
                    help='gro or similar file containing a martini structure')
parser.add_argument(dest='P', type=valid_top, default=None, nargs='?',
                    help='top or tpr file with the topology of the system')
parser.add_argument(dest='T', type=valid_traj, default=None, nargs='?',
                    help='corresponding trajectory file')
parser.add_argument('-s', '--skip', type=int, default=None)
parser.add_argument('-g', '--gmx', type=str, dest='G', default=None,
                    help='path to the gromacs executable')

args = parser.parse_args()

__main__.pymol_argv = ['pymol']
pymol.finish_launching()

this_script_dir = os.path.dirname(os.path.realpath(__file__))

mt_dir = os.path.realpath(os.path.join(this_script_dir, os.pardir))

cmd.run(os.path.join(mt_dir, 'pycg_bonds', 'pycg_bonds.py'))
cmd.run(os.path.join(mt_dir, 'mt_tools', 'mt_sele.py'))

cmd.load(args.S)
cmd.sync()
cmd.load_traj(args.T)
cmd.sync()
cmd.do(f'cg_bonds {args.P}')
cmd.sync()
cmd.do('mt_sele')
cmd.sync()
