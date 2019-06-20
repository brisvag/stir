from pymol import cmd, stored
import random


selectors = {
    'prot': 'polymer.protein',
    'BB': 'polymer.protein and name BB',
    'SC': 'polymer.protein and name SC*',
    'solv': 'resname W or resname WN or resname ION',
    'ions': 'resname ION',
    'lip': 'organic',
    'nucl': 'polymer.nucleic'
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
            for sel in selectors:
                cmd.delete(sel)
        else:
            print('Unknown option. Type `mt_sele delete` to delete all the mt selections')
    else:
        for sel, logic in selectors.items():
            cmd.select(sel, logic)


def nice_colors():
    """
    creates a set of nicely distinguishable colors and stores their pymol ids in stored.mt_colors
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

    stored.mt_colors = []
    for color in mt_colors:
        idx = cmd.get_color_index(color)
        stored.mt_colors.append(idx)


def mt_color(method, selection='all'):
    stored.tmp_dict = {}
    cmd.iterate(selection, f'stored.tmp_dict[{method}] = random.choice(stored.mt_colors)')
    cmd.alter(selection, f'color = stored.tmp_dict[{method}]')
    cmd.sync()
    stored.tmp_dict = {}
    cmd.recolor()


def mt_nice(selection='all'):
    # TODO: add multiple choice option, e.g: pretty|bilayer|vdw?
    mt_sele()

    cmd.set('stick_radius', 0.7)

    nice_settings = {
        'prot': {
            'color_method': None,
            'style': 'stick'
        },
        'BB': {
            'color_method': 'chain',
            'style': 'stick'
        },
        'SC': {
            'color_method': 'resi',
            'style': 'stick'
        },
        'solv': {
            'color_method': None,
            'style': 'sphere'
        },
        'ions': {
            'color_method': 'name',
            'style': 'sphere'
        },
        'lip': {
            'color_method': 'resi',
            'style': 'stick'
        },
        'nucl': {
            'color_method': 'chain',
            'style': 'stick'
        },
    }

    # TODO: a nicer way to do this?
    cmd.color('blue', f'{selection} and solv')

    for sel, data in nice_settings.items():
        if data['color_method']:
            mt_color(data['color_method'], f'{selection} and {sel}')
        if data['style']:
            cmd.show_as(data['style'], f'{selection} and {sel}')


def load():
    nice_colors()
    cmd.extend('mt_sele', mt_sele)
    cmd.extend('mt_color', mt_color)
    cmd.extend('mt_nice', mt_nice)
