#!/usr/bin/env python3

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

parser.add_argument(dest='struct', type=valid_str,
                    help='gro or similar file containing a martini structure')
parser.add_argument(dest='topol', type=valid_top, default=None, nargs='?',
                    help='top or tpr file with the topology of the system')
parser.add_argument(dest='traj', type=valid_traj, default=None, nargs='?',
                    help='corresponding trajectory file')
parser.add_argument('-s', '--skip', dest='skip', type=int, default=1,
                    help='when loading a trajectory, load frames with this rate')   # TODO
parser.add_argument('-g', '--gmx', dest='gmx', type=str, default=None,
                    help='path to the gromacs executable')  # TODO
parser.add_argument('--keepwater', dest='keepwater', action='store_true',
                    help='do not delete waters from the system. Decreases performance')
# TODO: add more options (load_traj start/end...)
# TODO: passing arguments to pymol

args = parser.parse_args()

# TODO: check size of trajectory and compare with ram size. If not a good ratio (?), warn the user
#       and ask for confirmation. Also, let them know they can use -s to skip frames

__main__.pymol_argv = ['pymol']
pymol.finish_launching()

this_script_dir = os.path.dirname(os.path.realpath(__file__))
mt_dir = os.path.realpath(os.path.join(this_script_dir, os.pardir))

cmd.run(os.path.join(mt_dir, 'config_files', 'pymolrc.py'))
cmd.run(os.path.join(mt_dir, 'pycg_bonds', 'pycg_bonds.py'))
cmd.run(os.path.join(mt_dir, 'mt_tools', 'mt_sele.py'))
cmd.run(os.path.join(mt_dir, 'mt_tools', 'mt_supercell.py'))
cmd.run(os.path.join(mt_dir, 'mt_tools', 'mt_movie.py'))

cmd.load(args.struct)
if not args.keepwater:
    cmd.remove('resname W or resname WN')
cmd.sync()

cmd.sync()

if args.traj:
    traj_args = [args.traj]
    if not args.keepwater:
        traj_args.append('selection="not resname W and not resname WN"')
    cmd.run(os.path.join(mt_dir, 'config_files', 'trajectory.py'))
    cmd.sync()
    cmd.load_traj(*traj_args)
cmd.sync()

cg_bond_args = []
if args.topol:
    cg_bond_args.append(args.topol)
if args.gmx:
    cg_bond_args.append(f'gmx={args.gmx}')
cg_bond_args = ', '.join(cg_bond_args)

cmd.do(f'cg_bonds {cg_bond_args}')
cmd.sync()

cmd.do('mt_sele')
cmd.sync()
if not args.keepwater:
    cmd.delete('solv')

mt_help = '''
Martini Tools:
    
help cg_bonds
help mt_sele
help mt_supercell
help mt_movie
'''

print(mt_help)
