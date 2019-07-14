# Copyright 2019-2019 the stir authors. See copying.md for legal info.

from pymol import cmd, movie

# local imports
from . import config


def mt_movie(movie_type, savefile=None, duration=5, mode='ray', width=1920, height=1080):
    """
DESCRIPTION

    loads nice settings for rendering and records a movie of the trajectory.
    Before running, position the camera to have a good point of view.

USAGE

    mt_movie movie_type [, savefile [, duration [, mode [, width [, height]]]]]

ARGUMENTS

    movie_type = one of:
        - traj: records the whole trajectory
        - bullettime: a whole revolution around the z axis
    savefile = file name to use to save the movie, without extension. If not given,
               the movie won't be save but only shown (default=None)
    duration = duration in seconds of the final movie (default=5)
    ray = set to 'draw' to disable ray tracing (default='ray')
    width = width of the movie in pixel (default=1920)
    height = height of the movie in pixel (default=1080)
    """
    duration = int(duration)

#    # store settings for later restoration
#    store_settings.load()
#    cmd.sync()
#    settings_file = '/tmp/mt_settings.py'
#    cmd.do(f'store_settings {settings_file}')

    # load nice settings for rendering
    config.rendering()
    cmd.sync()

    # initialize movie frames
    cmd.mdelete()
    cmd.mview('reset')
    frames = duration*30
    cmd.mset(f'1x{frames}')

    # make a still movie of the whole trajectory
    if movie_type == 'traj':
        cmd.mview('store', 1, state=1)
        cmd.sync()
        states = cmd.count_states()
        cmd.mview('store', frames, state=states)
        cmd.sync()
        cmd.mview('reinterpolate')
        cmd.sync()

    # make a movie rotating 360 around the z axis
    if movie_type == 'bullettime':
        cmd.mview('store', 1)
        cmd.mview('store', frames)
        # TODO: moving objects seems a bit hacky. Any better solution?
        for obj in cmd.get_object_list():
            cmd.mview('store', object=obj)
            cmd.sync()
        for i in range(4):
            frame = int(frames * (i+1)/4)
            cmd.frame(frame)
            cmd.sync()
            for obj in cmd.get_object_list():
                cmd.rotate('z', angle=-90, camera=0, object=obj)
                cmd.sync()
                cmd.mview('store', object=obj)
                cmd.sync()
                cmd.mview('reinterpolate', object=obj, power=1)
                cmd.sync()
        cmd.mview('reinterpolate')
        cmd.sync()
    else:
        raise ValueError(f'{movie_type} is not a valid argument. See `help mt_movie`.')

    # save if required, otherwise just play it
    if savefile:
        cmd.viewport(width, height)
        movie.produce(f'{savefile}.mp4', mode=mode, encoder='ffmpeg', quality=100)
    else:
        cmd.mplay()

    # restore previous settings
#    cmd.do(f'run {settings_file}')


def load():
    cmd.extend('mt_movie', mt_movie)


if __name__ == 'pymol':
    load()
