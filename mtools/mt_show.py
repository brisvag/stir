from pymol import cmd

# local imports
from pycg_bonds.pycg_bonds import cg_bonds
from
from mtools import mt_sele


#script_dir = os.path.dirname(os.path.realpath(__file__))
#pycgb = os.path.join(script_dir, 'pycg_bonds/pycg_bonds.py')

def mt_show(topfile=None, selection='all'):
    print('stuff is happening!')
    cg_bonds(topfile, selection)
    mt_sele()


cmd.extend('mt_show', mt_show)
