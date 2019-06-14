from pymol import cmd

# local imports
import pycg_bonds.pycg_bonds
import mtools.mt_sele


#script_dir = os.path.dirname(os.path.realpath(__file__))
#pycgb = os.path.join(script_dir, 'pycg_bonds/pycg_bonds.py')

def mt_show(topfile=None, selection='all'):
    print('stuff is happening!')
    cmd.cg_bonds(topfile, selection)
    cmd.mt_sele()


cmd.extend('mt_show', mt_show)
