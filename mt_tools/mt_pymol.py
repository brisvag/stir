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
from skewer import skewer

# local imports
from mt_tools import config
from mt_tools import mt_movie, mt_nice, mt_supercell
from mt_tools.utils import valid_str, valid_top, valid_traj, clean_path


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


# check if there's enough memory to load the requested traj and warn the user if needed
if args.traj:
    freemem = psutil.virtual_memory().available
    traj_size = 0
    for traj in args.traj:
        traj_size += os.path.getsize(clean_path(traj))
    water_ratio = 1
    if not args.keepwater:
        # TODO: VERY arbitrary number. When skewer's parsing is a module, use that!
        water_ratio = 1/2
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
                print(f'"{inp}" is not a valid choice. [y/N]')

# initialize pymol
__main__.pymol_argv = ['pymol']
pymol.finish_launching()

# run pymolrc and load all the mt_tools
config.pymolrc()
mt_nice.load()
mt_supercell.load()
mt_movie.load()

# load skewer
skewer.load()
cmd.sync()

# open the structure
cmd.load(clean_path(args.struct))
cmd.sync()
# get the loaded object's name, so we can load the traj into it as new states
sys_obj = cmd.get_object_list()[0]

# load trajectories
if args.traj:
    config.trajectory()
    for traj in args.traj:
        cmd.sync()
        cmd.load_traj(clean_path(traj), sys_obj, interval=args.skip)
    cmd.sync()

# TODO: "selection" in load_traj seems not to work as planned. Can we get it to work?
#       Other option: call trjconv to get rid of the waters before loading
# delete waters, unless they are needed
if not args.keepwater:
    cmd.remove('resname W or resname WN')
    cmd.sync()

# run skewer with as many arguments as we got
skewer_args = []
if args.topol:
    skewer_args.append(clean_path(args.topol))
if args.gmx:
    skewer_args.append(f'gmx={args.gmx}')
skewer_args = ', '.join(skewer_args)

cmd.do(f'skewer {skewer_args}')
cmd.sync()

# run mt_nice with the `clean` setting
cmd.do(f'mt_nice not *_elastics')
cmd.sync()

# print some help after everything is loaded
mt_help = '''
Martini Tools functions:

- skewer
- mt_nice, mt_sele, mt_color
- mt_supercell
- mt_movie
'''
cmd.sync()
print(mt_help)
