"""
initial configuration file for pymol
"""

from pymol import cmd
import psutil


def config():
    # multithreading
    cmd.set('max_threads', psutil.cpu_count())

    # color space
    cmd.space('pymol')

    # nicer visualization
    cmd.set('valence', 'on')
    cmd.set('cartoon_side_chain_helper', 'on')
    cmd.set('cartoon_highlight_color', 'grey50')
    cmd.set('sphere_mode', '9')
