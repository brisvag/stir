# Copyright 2019-2019 the stir authors. See copying.md for legal info.

"""
commodity functions to make taking pictures and making movies easier and faster
"""

from pymol import cmd, movie

# local imports
from . import config


def cheese(render_type, savefile=None, duration=5, mode='ray', width=1920, height=1080):
    """
DESCRIPTION

    loads nice settings for rendering and  take a picture or records a movie of the trajectory.
    Before running, position the camera to have a good point of view.

USAGE

    cheese render_type [, savefile [, duration [, mode [, width [, height]]]]]

ARGUMENTS

    render_type = one of:
        - snap: take a snapshot of the current view
        - traj: record the whole trajectory
        - bullettime: make a whole revolution around the z axis
    savefile = file name to use to save the movie, without extension. If not given,
               the movie won't be save but only shown (default=None)
    duration = duration in seconds of the final movie (default=5)
    mode = set to 'draw' to disable ray tracing (default='ray')
    width = width of the movie in pixel (default=1920)
    height = height of the movie in pixel (default=1080)
    """
    duration = int(duration)

#    # store settings for later restoration
#    store_settings.load()
#    cmd.sync()
#    settings_file = '/tmp/stir_settings.py'
#    cmd.do(f'store_settings {settings_file}')

    # load nice settings for rendering
    config.rendering()
    cmd.sync()

    def init_movie():
        cmd.mdelete()
        cmd.mview('reset')
        frames = duration*30
        cmd.mset(f'1x{frames}')
        return frames

    # if requested, take a picture
    if render_type == 'snap':
        cmd.ray(width=width, height=height)
    # make a still movie of the whole trajectory
    elif render_type == 'traj':
        frames = init_movie()
        cmd.mview('store', 1, state=1)
        cmd.sync()
        states = cmd.count_states()
        cmd.mview('store', frames, state=states)
        cmd.sync()
        cmd.mview('reinterpolate')
        cmd.sync()
    # make a movie rotating 360 around the z axis
    elif render_type == 'bullettime':
        frames = init_movie()
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
        raise ValueError(f'{render_type} is not a valid argument. See `help cheese`.')

    # save if required, otherwise just play it
    if savefile:
        cmd.viewport(width, height)
        movie.produce(f'{savefile}.mp4', mode=mode, encoder='ffmpeg', quality=100)
    else:
        cmd.mplay()

    # restore previous settings
#    cmd.do(f'run {settings_file}')


def load():
    """
    adds cheese command to pymol
    """
    cmd.extend('cheese', cheese)


if __name__ == 'pymol':
    load()
