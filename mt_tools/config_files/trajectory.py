"""
configuration file for trajectory visualization
"""

from pymol import cmd


def config():
    # improves performance drastically when loading a lot of states
    cmd.set('defer_builds_mode', 3)
