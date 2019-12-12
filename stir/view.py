# Copyright 2019-2019 the stir authors. See copying.md for legal info.

"""
Tools and functions for better visualization and selection
"""

from pymol import cmd, stored
import random
import re


def nice_settings():
    """
    sets up a series of settings and stores in pymol.stored to be used at later points by stir

    parts:
    - a set of nicely distinguishable colors for nicecolor
    - a dictionary of selection expressions for nicesele
    - a dictionary of representation information for nice
    """
    # colors are based on this SO answer:
    # https://graphicdesign.stackexchange.com/questions/3682/
    # where-can-i-find-a-large-palette-set-of-contrasting-colors-for-coloring-many-d
    colors = [
        (240, 163, 255),
        (0, 117, 220),
        (153, 63, 0),
        (76, 0, 92),
        (25, 25, 25),
        (0, 92, 49),
        (43, 206, 72),
        (255, 204, 153),
        (128, 128, 128),
        (148, 255, 181),
        (143, 124, 0),
        (157, 204, 0),
        (194, 0, 136),
        (0, 51, 128),
        (255, 164, 5),
        (255, 168, 187),
        (66, 102, 0),
        (255, 0, 16),
        (94, 241, 242),
        (0, 153, 143),
        (224, 255, 102),
        (116, 10, 255),
        (153, 0, 0),
        (255, 255, 128),
        (255, 255, 0),
        (255, 80, 5),
    ]

    nicecolors = {}
    for i, rgb in enumerate(colors):
        col_name = f'nicecolor_{i}'
        cmd.set_color(col_name, rgb)
        nicecolors[col_name] = rgb

    # store color ids in a list to be used by nicecolor
    stored.nicecolors = []
    for color in nicecolors:
        idx = cmd.get_color_index(color)
        stored.nicecolors.append(idx)

    # selection expressions
    stored.niceselectors = {
        'prot': 'polymer.protein',
        'BB': 'polymer.protein and name BB',
        'SC': 'polymer.protein and not name BB',
        'solv': 'resn W+WN+ION',
        'ions': 'resn ION',
        'lip': 'organic and not ions',
        'nucl': 'polymer.nucleic'
    }

    # settings to be used by nice. Values are lists with a function as first element and
    # its arguments following
    stored.nice_set = {
        'clean': {
            'prot': {
                'color': None,
                'style': None,
            },
            'BB': {
                'color': [nicecolor, 'chain'],
                'style': [cmd.show_as, 'sticks'],
            },
            'SC': {
                'color': None,
                'style': [cmd.hide, 'everything'],
            },
            'solv': {
                'color': None,
                'style': [cmd.hide, 'everything'],
            },
            'ions': {
                'color': None,
                'style': [cmd.hide, 'everything'],
            },
            'lip': {
                'color': [nicecolor, 'resn'],
                'style': [cmd.show_as, 'sticks'],
            },
            'nucl': {
                'color': [nicecolor, 'resi'],
                'style': [cmd.show_as, 'sticks'],
            },
        },
        'rainbow': {
            'prot': {
                'color': None,
                'style': None,
            },
            'BB': {
                'color': [nicecolor, 'chain'],
                'style': [cmd.show_as, 'sticks'],
            },
            'SC': {
                'color': [nicecolor, 'resn'],
                'style': [cmd.show_as, 'sticks'],
            },
            'solv': {
                'color': [cmd.color, 'blue'],
                'style': [cmd.show_as, 'nb_spheres'],
            },
            'ions': {
                'color': [nicecolor, 'name'],
                'style': [cmd.show_as, 'nb_spheres'],
            },
            'lip': {
                'color': [nicecolor, 'resi'],
                'style': [cmd.show_as, 'sticks'],
            },
            'nucl': {
                'color': [nicecolor, 'resn'],
                'style': [cmd.show_as, 'sticks'],
            },
        },
        'balls': {
            'prot': {
                'color': None,
                'style': None,
            },
            'BB': {
                'color': [cmd.color, 'purple'],
                'style': [cmd.show_as, 'spheres'],
            },
            'SC': {
                'color': [cmd.color, 'red'],
                'style': [cmd.show_as, 'spheres'],
            },
            'solv': {
                'color': [cmd.color, 'blue'],
                'style': [cmd.show_as, 'nb_spheres'],
            },
            'ions': {
                'color': [nicecolor, 'resn'],
                'style': [cmd.show_as, 'nb_spheres'],
            },
            'lip': {
                'color': [nicecolor, 'resn'],
                'style': [cmd.show_as, 'spheres'],
            },
            'nucl': {
                'color': [nicecolor, 'resi'],
                'style': [cmd.show_as, 'spheres'],
            },
        },
    }


def nicesele(delete=None):
    """
DESCRIPTION

    easily select relevant groups of atoms in a martini system

USAGE

    nicesele [, delete]

ARGUMENTS

    delete = if set, will delete the selections instead of creating them
    """
    if delete:
        if delete == 'delete':
            for sel in stored.niceselectors:
                cmd.delete(sel)
        else:
            print('Unknown option. Type `nicesele delete` to delete all the pre-made selections')
    else:
        for sel, logic in stored.niceselectors.items():
            cmd.select(sel, logic)
    cmd.sync()
    # disable last selection to avoid accidental modifications to it
    cmd.deselect()
    cmd.sync()


def nicecolor(method, selection='all'):
    """
DESCRIPTION

    color a selection with unique colors based on atom properties exposed by iterate

USAGE

    nicecolor method [, selection]

ARGUMENTS

    method = property that will be used as a method for coloring (e.g.: "resi")
    selection = selection to act upon
    """
    stored.tmp_dict = {}
    stored.r_choice = random.choice
    cmd.iterate(selection, f'stored.tmp_dict[{method}] = stored.r_choice(stored.nicecolors)')
    cmd.alter(selection, f'color = stored.tmp_dict[{method}]')
    cmd.sync()
    stored.tmp_dict = {}
    cmd.recolor()
    cmd.sync()


def set_vdw(selection='all'):
    """
    alters vdw radii based on bead type provided by garnish
    """
    def alter_vdw(elem, vdw):
        vdw_patterns = {
            re.compile('[QPNCX][\w\d]|W'): 2.35,
            re.compile('S([QPNCX][\w\d]|W)'): 2.00,
            re.compile('T([QPNCX][\w\d]|W)'): 1.65,
        }
        for pattern, radius in vdw_patterns.items():
            if pattern.match(elem):
                return radius
        # TODO: if anything goes wrong, just return the original vdw for now
        return vdw

    stored.alter_vdw = alter_vdw

    cmd.alter(selection, 'vdw=stored.alter_vdw(elem, vdw)')
    cmd.sync()


def nice(style='clean', selection='all'):
    """
DESCRIPTION

    apply a nice preset for the visualization of a martini molecular system
    NOTE: ignores object names matching `*_elastics`

USAGE

    nice [style [, selection]]

ARGUMENTS

    style = one of clean|rainbow|balls (default='clean')
    selection (default='all')
    """
    # sanitize input
    if style not in stored.nice_set.keys():
        print(f'Error: {style} is not a valid option. Choose from: {list(stored.nice_set.keys())}')
        return

    selection = f'{selection} and not *_elastics'

    nicesele()

    cmd.set('stick_radius', 0.7)
    cmd.sync()
    # set correct vdw radii
    set_vdw(selection)
    cmd.sync()

    settings = stored.nice_set[style]
    for sel_type, commands in settings.items():
        for command in commands.values():
            if command:
                # run function with its arguments. All functions must have `selection` as
                # valid argument for this to work!
                command[0](*command[1:], selection=f'{selection} and {sel_type}')
            else:
                continue


def load():
    """
    adds nice, nicesele and nicecolor commands to pymol
    """
    nice_settings()
    cmd.extend('nicesele', nicesele)
    cmd.extend('nicecolor', nicecolor)
    cmd.extend('nice', nice)

    # tab completion
    cmd.auto_arg[0]['nice'] = [lambda: cmd.Shortcut(['clean', 'rainbow', 'balls']), 'style', '']
    cmd.auto_arg[1]['nice'] = [cmd.selection_sc, 'selection', '']
