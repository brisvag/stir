from pymol import cmd, stored

# local imports
from mt_tools import config


def mt_movie(movie_type, duration=5):
    """
    loads nice settings for rendering and records a movie of the trajectory
    before running, position the camera to have a good point of view

    arguments
    movie_type: one of traj|matrix
                - traj: records the whole trajectory
                - matrix: a whole revolution around the z axis
    duration: duration in seconds of the final movie
    """
    duration = int(duration)

    #    settings_file = '/tmp/mt_settings.py'
    #
    #    # store settings for later restoration
    #    cmd.do(f'save_settings {settings_file}')

    config.rendering()
    cmd.sync()

    frames = duration*30
    cmd.mset(f'1x{frames}')

    if movie_type == 'traj':
        cmd.mview('store', 1, state=1)
        cmd.sync()
        states = cmd.count_states()
        cmd.mview('store', frames, state=states)
        cmd.sync()
        cmd.mview('reinterpolate')
        cmd.sync()
    if movie_type == 'matrix':
        cmd.mview('store', 1)
        cmd.sync()
        for i in range(4):
            frame = ((frames//4)*(i+1))
            cmd.frame(frame)
            cmd.sync()
            cmd.rotate('z', angle=90, camera=0)
            cmd.sync()
            cmd.mview('store', 0)
            cmd.sync()
        cmd.mview('reinterpolate')
        cmd.sync()

    # restore previous settings
#    cmd.do(f'run {settings_file}')


def load():
    cmd.extend('mt_movie', mt_movie)


if __name__ == 'pymol':
    load()
