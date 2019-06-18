from pymol import cmd


selectors = {
    'prot': 'polymer.protein',
    'BB': 'polymer.protein and name BB',
    'SC': 'polymer.protein and name SC*',
    'solv': 'resname W or resname WN or resname ION',
    'ions': 'resname ION',
    'lip': 'organic',
    'nucl': 'polymer.nucleic'
}


def mt_sele():
    """
    easily select relevant groups of atoms

    Usage: mt_sele

    See also: mt_desele
    """
    for sel, logic in selectors.items():
        cmd.select(sel, logic)
        if cmd.count_atoms(sel) == 0:
            cmd.delete(sel)


def mt_desele():
    """
    deletes all selections created by mt_sele

    Usage: mt_desele

    See also: mt_sele
    """
    for sel in selectors:
        cmd.delete(sel)


cmd.extend('mt_sele', mt_sele)
cmd.extend('mt_desele', mt_desele)
