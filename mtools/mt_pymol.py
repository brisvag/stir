# Copyright 2019-2019 the mtools authors. See copying.md for legal info.

"""
Python wrapper for pymol to easily visualize martini trajectories
"""

import os
import sys
import argparse
import pymol
from pymol import cmd
import __main__
import psutil
from garnish import garnish

# local imports
from . import config
from . import mt_movie, mt_nice, mt_supercell
from .utils import valid_str, valid_top, valid_traj, clean_path


class MyParser(argparse.ArgumentParser):
    # print help if calling the program results in an error
    def error(self, message):
        sys.stderr.write(f'ERROR: {message}\n\n')
        self.print_help()
        sys.exit(2)


def main(*args, **kwargs):
    parser = MyParser(prog='mt_pymol', description='A python wrapper of martinitools for pymol.',
                      add_help=False)

    help_group = parser.add_argument_group('HELP')
    help_group.add_argument('-h', '--help', action='help',
                            help='show this help message and exit')

    req_group = parser.add_argument_group('required arguments')
    req_group.add_argument(dest='struct', type=valid_str,
                           help='gro or similar file containing a martini structure')

    pos_group = parser.add_argument_group('positional arguments')
    pos_group.add_argument(dest='topol', type=valid_top, default=None, nargs='?',
                           help='top or tpr file with the topology of the system')
    pos_group.add_argument(dest='traj', type=valid_traj, default=None, nargs='*',
                           help='corresponding trajectory file. If multiple files are given, '
                                'they are concatenated')

    opt_group = parser.add_argument_group('optional mt_pymol arguments')
    opt_group.add_argument('-g', '--gmx', dest='gmx', type=str, default=None,
                           help='path to the gromacs executable')
    opt_group.add_argument('--keepwater', dest='keepwater', action='store_true',
                           help='do not delete waters from the system. Decreases performance')

    traj_group = parser.add_argument_group('trajectory arguments')
    traj_group.add_argument('-s', '--skip', dest='skip', type=int, default=1,
                            help='when loading a trajectory, load frames with this rate')

    more_group = parser.add_argument_group('additional arguments')
    more_group.add_argument('-p', '--pymol', dest='pymol', default=[], nargs=argparse.REMAINDER,
                            help='remaining arguments will be passed to pymol. Accepts options '
                                 'and .pml scripts')
    # TODO: add more options:
    #       - load_traj start/end...
    #       - running subscripts automatically (mt_supercell...)

    args = parser.parse_args()

    # check if there's enough memory to load the requested traj and warn the user if needed
    if args.traj:
        freemem = psutil.virtual_memory().available
        traj_size = 0
        for traj in args.traj:
            traj_size += os.path.getsize(clean_path(traj))
        water_ratio = 1
        if not args.keepwater:
            # TODO: VERY arbitrary number. When garnish's parsing is a module, use that!
            #       EDIT: I will probably leave it like this. Unnecessary and complex to use parse_tpr
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

    pymol_args = []
    scripts = []
    for arg in args.pymol:
        if os.path.splitext(arg)[1] in ('.pml', '.py'):
            scripts.append(clean_path(arg))
        else:
            pymol_args.append(arg)

    # initialize pymol
    __main__.pymol_argv = ['pymol'] + pymol_args
    pymol.finish_launching()

    # run pymolrc and load all the mtools
    config.pymolrc()
    mt_nice.load()
    mt_supercell.load()
    mt_movie.load()

    # load garnish
    garnish.extend_garnish()
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

    # run garnish with as many arguments as we got
    garnish_args = []
    if args.topol:
        garnish_args.append(clean_path(args.topol))
    if args.gmx:
        garnish_args.append(f'gmx={args.gmx}')
    garnish_args = ', '.join(garnish_args)

    cmd.do(f'garnish {garnish_args}')
    cmd.sync()

    # run mt_nice with the `clean` setting
    cmd.do(f'mt_nice not *_elastics')
    cmd.sync()

    # finally run user-provided scripts
    for scr in scripts:
        cmd.run(scr)

    # print some help after everything is loaded
    mt_help = '''
    Martini Tools functions:

    - garnish
    - mt_nice, mt_sele, mt_color
    - mt_supercell
    - mt_movie
    '''
    cmd.sync()
    print(mt_help)
