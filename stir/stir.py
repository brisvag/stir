# Copyright 2019-2019 the stir authors. See copying.md for legal info.

"""
Python wrapper for pymol to easily visualize martini trajectories
"""

import sys
import argparse
import pymol
from pymol import cmd
import __main__
import garnish

# local imports
from . import __version__
from . import config
from . import render, view, supercell, edit
from .utils import clean_path, stir_help


class FilesAction(argparse.Action):
    def __call__(self, parser, namespace, values, option_strings):
        types = {
            'struct': ('.gro', '.pdb',),
            'topol': ('.top', '.itp', '.tpr',),
            'traj': ('.xtc',),
            'scene': ('.pse',),
        }
        sorted = {k: [] for k in types}
        for v in values:
            p = clean_path(v)
            if not p.is_file():
                parser.error(f'{p} does not exist.')
            for file, ext in types.items():
                if p.suffix in ext:
                    sorted[file].append(p)
                    break
            else:
                parser.error(f'{p.suffix} is not a valid file type.')

        setattr(namespace, self.dest, sorted)


class HelpfulParser(argparse.ArgumentParser):
    def error(self, message):
        # print usage and red error message
        sys.stderr.write(f'\033[91mERROR: {message}\033[0m\n\n')
        self.print_usage()
        sys.exit(2)


def main():
    """
    parses arguments user input and initializes stir based on user input
    provides help for usage
    """
    parser = HelpfulParser(
        prog='stir', formatter_class=argparse.RawDescriptionHelpFormatter, add_help=False,
        description='A python wrapper for pymol and several helpful tools and scripts\n'
                    'mainly focused on martini coarse-grained trajectories.\n\n'
                    'The accepted file formats are:\n'
                    '- structure: gro, pdb\n'
                    '- scene: pse\n'
                    '- topology: top, itp, tpr\n'
                    '- trajectory: xtc',
        epilog='Examples:\n'
               '\tstir system.gro topol.top md.xtc\n'
               '\tstir system.gro --keep-water -r supercell 3,3,1\n'
               '\tstir system.gro topol.tpr --pymol -qi myscript.pml',
        )

    files_group = parser.add_argument_group('file arguments')
    files_group.add_argument(dest='files', action=FilesAction, nargs='+',
                             help='a structure or scene file is required. Topology files allow for topology '
                                  'reconstruction. Trajectory files are also accepted. Order does not matter')

    opt_group = parser.add_argument_group('optional stir arguments')
    opt_group.add_argument('--keep-water', dest='keepwater', action='store_true',
                           help='do not delete waters from the system. Decreases performance')
    opt_group.add_argument('-g', '--gmx', dest='gmx', type=str, default=None,
                           help='path to the gromacs executable')
    opt_group.add_argument('-r', '--run-tool', dest='runtool', metavar='token', type=str, default=[], nargs='*',
                           action='append',
                           help='a command to be run after loading. (e.g.: supercell 3,3,1). '
                                'Can be specified multiple times')

    gar_group = parser.add_argument_group('optional garnish arguments')
    gar_group.add_argument('--no-fix', dest='nofix', action='store_false',  # store false!
                           help='disable the atom-id-based fix for the elastic network in garnish '
                                '(use if your system has messy, non-sequential numbering.')
    gar_group.add_argument('--no-prot', dest='noprot', action='store_false',    # store false!
                           help='do not guess protein backbone beads (use if normal guessing makes mistakes')
    gar_group.add_argument('--no-garnish', dest='nogarnish', action='store_true',
                           help='do not run garnish on the system (use with atomistic systems)')

    traj_group = parser.add_argument_group('optional trajectory arguments')
    traj_group.add_argument('-s', '--skip', dest='skip', metavar='n', type=int, default=1,
                            help='load trajectory frames skipping this interval. Useful to reduce memory load')
    traj_group.add_argument('-b', '--begin', dest='begin', metavar='frame', type=int, default=1,
                            help='first frame to load from trajectory. Only acts on first trajectory file')
    traj_group.add_argument('-e', '--end', dest='end', metavar='frame', type=int, default=-1,
                            help='last frame to load from trajectory. Only acts on last trajectory file')
    traj_group.add_argument('-m', '--max', dest='max', metavar='n', type=int, default=0,
                            help='maximum number of frames to load')

    more_group = parser.add_argument_group('pymol arguments')
    more_group.add_argument('-p', '--pymol', dest='pymol', default=[], nargs=argparse.REMAINDER,
                            help='all following arguments will be passed directly to pymol. '
                                 'Accepts options and .pml scripts')

    help_group = parser.add_argument_group('info')
    help_group.add_argument('-h', '--help', action='help',
                            help='show this help message and exit')
    help_group.add_argument('-V', '--version', action='version', version=f'%(prog)s {__version__}')

    args = parser.parse_args()

    struct = args.files['struct']
    scene = args.files['scene']
    topol = args.files['topol']
    traj = args.files['traj']

    # make sure we only have 1 structure OR scene file and at most one topol
    if not bool(struct) ^ bool(scene):  # not xor
        parser.error('you must provide either a structure or scene file')
    elif len(struct) > 1 or len(scene) > 1 or len(topol) > 1:
        parser.error('only one system can be opened at once... for now!')  # TODO?

    # sanitize trajectory args
    if args.skip < 1:
        args.skip = 1
    if args.begin < 1:
        args.begin = 1
    if args.end < 1 and args.end != -1:
        args.end = 1
    if args.max < 0:
        args.max = 0

    pymol_args = []
    scripts = []
    for arg in args.pymol:
        p = clean_path(arg)
        if p.suffix in ('.pml', '.py'):
            if not p.is_file():
                raise FileNotFoundError(f'{p} does not exist')
            scripts.append(str(p))
        else:
            pymol_args.append(str(p))

    # initialize pymol
    __main__.pymol_argv = ['pymol'] + pymol_args
    pymol.finish_launching()

    # run pymolrc and load all the stir tools
    config.pymolrc()
    view.load()
    supercell.load()
    render.load()
    edit.load()

    # load garnish
    garnish.extend_garnish()
    cmd.sync()

    # open the structure
    if scene:
        cmd.load(scene[0])
    elif struct:
        cmd.load(struct[0])
    cmd.sync()
    # get the loaded object's name, so we can load the traj into it as new states
    sys_obj = cmd.get_object_list()[0]

    # load trajectories, leaving out waters if not asked for
    if traj:
        skip = args.skip
        max_states = args.max
        selection = 'all'
        if not args.keepwater:
            selection = 'not resname W+WN'
        config.trajectory()
        for i, t in enumerate(traj):
            cmd.sync()
            start = 1
            if i == 0:
                start = args.begin
            stop = -1
            if i == len(traj) - 1:
                stop = args.end
            if args.max != 0:
                max_states = args.max - cmd.count_states()
                if max_states < 1:
                    break
            cmd.load_traj(t, sys_obj, interval=skip, start=start, stop=stop,
                          max=max_states, selection=selection)
        cmd.sync()

    # also, delete waters from first frame
    if not args.keepwater:
        cmd.remove('resname W+WN')
        cmd.sync()

    if not args.nogarnish:
        # sanitize topol
        if not topol:
            topol = [None]
        garnish.garnish(file=topol[0], gmx=args.gmx, fix_elastics=args.nofix,
                        guess_prot=args.noprot, show=False)
        cmd.sync()

        # load garnish data into pymol
        view.nicesele()
        cmd.sync()
        view.set_vdw()
        cmd.sync()
        view.set_chains()
        cmd.sync()

        # run nice with the default settings, or with balls if no topol was given
        view.nice()
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

    # add command for simple help function
    cmd.extend('stir', stir_help)

    # print some help after everything is loaded
    stir_help()
    cmd.sync()
