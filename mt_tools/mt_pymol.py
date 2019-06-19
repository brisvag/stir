#!/usr/bin/env python3

"""
Python wrapper for pymol to easily visualize martini trajectories
"""

import os
import argparse
import pymol
from pymol import cmd
import __main__
import psutil


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
parser.add_argument(dest='traj', type=valid_traj, default=None, nargs='*',
                    help='corresponding trajectory file. If multiple files are given, '
                         'they are concatenated')
parser.add_argument('-s', '--skip', dest='skip', type=int, default=1,
                    help='when loading a trajectory, load frames with this rate')
parser.add_argument('-g', '--gmx', dest='gmx', type=str, default=None,
                    help='path to the gromacs executable')
parser.add_argument('--keepwater', dest='keepwater', action='store_true',
                    help='do not delete waters from the system. Decreases performance')
# TODO: add more options (load_traj start/end...)
# TODO: passing arguments to pymol

args = parser.parse_args()

if args.traj:
    freemem = psutil.virtual_memory().available
    traj_size = 0
    for traj in args.traj:
        traj_size += os.path.getsize(traj)
    water_ratio = 1
    if not args.keepwater:
        # TODO: VERY arbitrary number. When pycg_bonds parsing is a module, use that!
        water_ratio = 1/3
    # check if there's enough free memory: 5 is based on some testing
    if freemem < 5*(traj_size/args.skip):
        ok = False
        inp = input('WARNING: You may not have enough free memory to open this big trajectory.\n'
                    'Consider using the trajectory options (-s, ...).\n'
                    'Otherwise, continue at your own risk ;) [y/N] ')
        while not ok:
            if inp.lower() in ['yes', 'y']:
                ok = True
            elif inp.lower() in ['no', 'n']:
                parser.print_help()
                exit(0)
            else:
                print(f'"{inp}" is not a valid choice.')

__main__.pymol_argv = ['pymol']
pymol.finish_launching()

this_script_dir = os.path.dirname(os.path.realpath(__file__))
mt_dir = os.path.realpath(os.path.join(this_script_dir, os.pardir))

cmd.run(os.path.join(mt_dir, 'config_files', 'pymolrc.py'))
cmd.run(os.path.join(mt_dir, 'pycg_bonds', 'pycg_bonds.py'))
cmd.run(os.path.join(mt_dir, 'mt_tools', 'mt_sele.py'))
cmd.run(os.path.join(mt_dir, 'mt_tools', 'mt_supercell.py'))
cmd.run(os.path.join(mt_dir, 'mt_tools', 'mt_movie.py'))
cmd.sync()

cmd.load(args.struct)
cmd.sync()
sys_obj = cmd.get_object_list()[0]

if args.traj:
    cmd.run(os.path.join(mt_dir, 'config_files', 'trajectory.py'))
    for traj in args.traj:
        cmd.sync()
        cmd.load_traj(traj, sys_obj, interval=args.skip)
    cmd.sync()

# TODO: "selection" in load_traj seems not to work as planned. Can we get it to work?
#       Other option: call trjconv to get rid of the waters before loading
if not args.keepwater:
    cmd.remove('resname W or resname WN')
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

cmd.sync()
print(mt_help)
