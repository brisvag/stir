from pymol import cmd, stored
import random


def nice_settings():
    """
    sets up a series of settings and stores in pymol.stored to be used at later points by mt_tools

    parts:
    - a set of nicely distinguishable colors for mt_color
    - a dictionary of selection expressions for mt_sele
    - a dictionary of representation information for mt_nice
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

    mt_colors = {}
    for i, rgb in enumerate(colors):
        col_name = f'mt_color_{i}'
        cmd.set_color(col_name, rgb)
        mt_colors[col_name] = rgb

    # store color ids in a list to be used by mt_color
    stored.mt_colors = []
    for color in mt_colors:
        idx = cmd.get_color_index(color)
        stored.mt_colors.append(idx)

    # selection expressions
    stored.mt_selectors = {
        'prot': 'polymer.protein',
        'BB': 'polymer.protein and name BB',
        'SC': 'polymer.protein and name SC*',
        'solv': 'resn W or resn WN or resn ION',
        'ions': 'resn ION',
        'lip': 'organic',
        'nucl': 'polymer.nucleic'
    }

    # settings to be used by mt_nice. Values are lists with a function as first element and
    # its arguments following
    stored.mt_nice_set = {
        'rainbow': {
            'prot': {
                'color_method': None,
                'style': None,
            },
            'BB': {
                'color_method': [mt_color, 'chain'],
                'style': [cmd.show_as, 'sticks'],
            },
            'SC': {
                'color_method': [mt_color, 'resn'],
                'style': [cmd.show_as, 'sticks'],
            },
            'solv': {
                'color_method': [cmd.color, 'blue'],
                'style': [cmd.hide, 'everything'],
            },
            'ions': {
                'color_method': [mt_color, 'resn'],
                'style': [cmd.show_as, 'spheres'],
            },
            'lips': {
                'color_method': [mt_color, 'resi'],
                'style': [cmd.show_as, 'sticks'],
            },
            'nucl': {
                'color_method': [mt_color, 'resn'],
                'style': [cmd.show_as, 'sticks'],
            },
        },
        'clean': {
            'prot': {
                'color_method': None,
                'style': None,
            },
            'BB': {
                'color_method': [mt_color, 'chain'],
                'style': [cmd.show_as, 'sticks'],
            },
            'SC': {
                'color_method': None,
                'style': [cmd.hide, 'everything'],
            },
            'solv': {
                'color_method': None,
                'style': [cmd.hide, 'everything'],
            },
            'ions': {
                'color_method': None,
                'style': [cmd.hide, 'everything'],
            },
            'lip': {
                'color_method': [mt_color, 'resn'],
                'style': [cmd.show_as, 'sticks'],
            },
            'nucl': {
                'color_method': [mt_color, 'resi'],
                'style': [cmd.show_as, 'sticks'],
            },
        },
    }


def mt_sele(delete=None):
    """
    easily select relevant groups of atoms

    Usage:

    to create selections: mt_sele
    to delete them: mt_sele delete

    See also: mt_desele
    """
    if delete:
        if delete == 'delete':
            for sel in stored.mt_selectors:
                cmd.delete(sel)
        else:
            print('Unknown option. Type `mt_sele delete` to delete all the mt selections')
    else:
        for sel, logic in stored.mt_selectors.items():
            cmd.select(sel, logic)


def mt_color(method, selection='all'):
    stored.tmp_dict = {}
    stored.r_choice = random.choice
    cmd.iterate(selection, f'stored.tmp_dict[{method}] = stored.r_choice(stored.mt_colors)')
    cmd.alter(selection, f'color = stored.tmp_dict[{method}]')
    cmd.sync()
    stored.tmp_dict = {}
    cmd.recolor()


def mt_nice(main_sele='all', style='clean'):
    mt_sele()

    cmd.set('stick_radius', 0.7)

    # TODO: a nicer way to do this?
    cmd.color('blue', f'{main_sele} and solv')

    settings = stored.mt_nice_set[style]
    for sele, commands in settings.items():
        for command in commands.values():
            if command:
                # run function with its arguments. All functions must have `selection` as
                # valid argument for this to work!
                command[0](*command[1:], selection=f'{main_sele} and {sele}')
            else:
                continue


def load():
    nice_settings()
    cmd.extend('mt_sele', mt_sele)
    cmd.extend('mt_color', mt_color)
    cmd.extend('mt_nice', mt_nice)
