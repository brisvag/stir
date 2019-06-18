from pymol import cmd
import os


def mt_movie():
    """
    loads nice settings for rendering
    records a nice movie of the trajectory and saves it as a file
    """
    this_script_dir = os.path.dirname(os.path.realpath(__file__))
    mt_dir = os.path.realpath(os.path.join(this_script_dir, os.pardir))

    settings_file = '/tmp/mt_settings.py'

    # store settings for later restoration
    cmd.do(f'save_settings {settings_file}')

    cmd.run(os.path.join(mt_dir, 'config_files', 'rendering.py'))

    # TODO: actually make a movie

    # restore previous settings
    cmd.do(f'run {settings_file}')


cmd.extend('mt_movie', mt_movie)
