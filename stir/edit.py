# Copyright 2019-2019 the stir authors. See copying.md for legal info.

"""
functions for editing of the system
"""

from pymol import cmd, stored

# local imports
from .view import nicesele


def split(rule='molecules', selection='all'):
    """
DESCRIPTION

    split the system into multiple objects according to a specific rule

USAGE

    split [, rule [, selection]]

ARGUMENTS

    rule = one of groups|molecules (default='groups')
           groups: splits protein, lipids, solvent, nucleic
           molecules: every molecule in a different object (also chains)
    selection (default='all')
    """
    # sanitize input
    if rule not in ('groups', 'molecules'):
        print(f'Error: {rule} is not a valid option. Choose from: groups, molecules')
        return

    nicesele()
    cmd.sync()

    if rule == 'groups':
        for obj in cmd.get_object_list(selection):
            groups = [sel for sel in ('prot', 'lip', 'solv', 'nucl') if cmd.count_atoms(f'{obj} and {sel}') > 0]
            for sel in groups:
                cmd.create(f'{obj}_{sel}', f'{obj} and {sel}')
                cmd.sync()

    elif rule == 'molecules':
        for obj in cmd.get_object_list(selection):
            stored.tmp_set = set()
            cmd.iterate(f'{selection} and {obj}', f'stored.tmp_set.add(int(segi))')
            cmd.sync()
            zeros = len(str(max(stored.tmp_set)))
            for mol_id in stored.tmp_set:
                cmd.create(f'{obj}_n{mol_id:0{zeros}}', f'segi {mol_id}')


def load():
    """
    adds split command to pymol
    """
    cmd.extend('split', split)

    # tab completion
    cmd.auto_arg[0]['split'] = [lambda: cmd.Shortcut(['groups', 'molecules']), 'rule', '']
    cmd.auto_arg[1]['split'] = [cmd.selection_sc, 'selection', '']
