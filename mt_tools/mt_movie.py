from pymol import cmd, stored
import os

from mt_tools.config_files import rendering

def mt_movie(duration=5, ):
    """
    loads nice settings for rendering and records a movie of the trajectory
    before running, position the camera to have a good point of view

    arguments
    duration: duration in seconds of the final movie

    """
    duration = int(duration)
#    this_script_dir = os.path.dirname(os.path.realpath(__file__))
#    print(this_script_dir)
#    mt_dir = os.path.realpath(os.path.join(this_script_dir, os.pardir))
    rendering.config()

    mt_dir = '/home/brisvag/git/mt_tools'

    settings_file = '/tmp/mt_settings.py'

    # store settings for later restoration
#    cmd.do(f'save_settings {settings_file}')

    cmd.run(os.path.join(mt_dir, 'config_files', 'rendering.py'))

    cmd.mset(f'1x{duration*30}')
    cmd.mview('store', 1, state=1)
    states = cmd.count_states()
    cmd.mview('store', duration*30, state=states)
    cmd.mview('reinterpolate')

    # restore previous settings
#    cmd.do(f'run {settings_file}')


def load():
    cmd.extend('mt_movie', mt_movie)
