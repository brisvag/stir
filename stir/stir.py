# Copyright 2019-2019 the stir authors. See copying.md for legal info.

"""
Python wrapper for pymol to easily visualize martini trajectories
"""

import os
import sys
import argparse
import pymol
from pymol import cmd
import __main__
from garnish import garnish

# local imports
from . import config
from . import render, view, supercell
from .utils import valid_str, valid_top, valid_traj, clean_path, enough_ram, stir_help


class MyParser(argparse.ArgumentParser):
    # print help if calling the program results in an error
    def error(self, message):
        sys.stderr.write(f'ERROR: {message}\n\n')
        self.print_help()
        sys.exit(2)


def main():
    parser = MyParser(prog='stir', description='A python wrapper for pymol and several tools and scripts.',
                      formatter_class=argparse.RawDescriptionHelpFormatter,
                      epilog='Examples:\n'
                             '\tstir system.gro topol.top md.xtc\n'
                             '\tstir system.gro --keep-water -r supercell 3,3,1\n'
                             '\tstir system.gro topol.tpr --pymol -qi myscript.pml',
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

    opt_group = parser.add_argument_group('optional stir arguments')
    opt_group.add_argument('--keep-water', dest='keepwater', action='store_true',
                           help='do not delete waters from the system. Decreases performance')
    opt_group.add_argument('-g', '--gmx', dest='gmx', type=str, default=None,
                           help='path to the gromacs executable')
    opt_group.add_argument('-r', '--run-tool', dest='runtool', type=str, default=[], nargs='*',
                           action='append',
                           help='a command to be run after loading. (e.g.: supercell 3,3,1).'
                                'Can be specified multiple times')

    traj_group = parser.add_argument_group('optional trajectory arguments')
    traj_group.add_argument('-s', '--skip', dest='skip', type=int, default=1,
                            help='load frames skipping this interval. Useful to reduce memory load')

    more_group = parser.add_argument_group('advanced arguments')
    more_group.add_argument('-p', '--pymol', dest='pymol', default=[], nargs=argparse.REMAINDER,
                            help='all following arguments will be passed directly to pymol. '
                                 'Accepts options and .pml scripts')
    # TODO: add more options:
    #       - load_traj start/end...

    args = parser.parse_args()

    if args.traj:
        if not enough_ram(args.traj, args.skip):
            ok = False
            while not ok:
                inp = input('WARNING: You may not have enough free memory to open '
                            'this big trajectory.\nConsider using the trajectory options to reduce '
                            'the memory load.\nOtherwise, proceed at your own risk ;) [y/N] ')
                if inp.lower() in ['yes', 'y']:
                    ok = True
                elif inp.lower() in ['no', 'n', '']:
                    parser.print_help()
                    exit(0)
                else:
                    print(f'ERROR: "{inp}" is not a valid choice')

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

    # run pymolrc and load all the stir
    config.pymolrc()
    view.load()
    supercell.load()
    render.load()

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

    # run nice with the `clean` setting
    cmd.do(f'nice selection="not *_elastics"')
    cmd.sync()

    # finally run user-requested tools
    for tool in args.runtool:
        command = ' '.join(tool)
        cmd.do(command)
        cmd.sync()
    # and user-provided scripts
    for scr in scripts:
        cmd.run(scr)
        cmd.sync()

    # print some help after everything is loaded
    stir_help()
    cmd.sync()
