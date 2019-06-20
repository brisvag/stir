from pymol import cmd, stored, movie

# local imports
from mt_tools import config
from mt_tools.utils import store_settings


def mt_movie(movie_type, savefile=None, duration=5, mode='ray'):
    """
    loads nice settings for rendering and records a movie of the trajectory
    before running, position the camera to have a good point of view

    arguments
    savefile: file name that will be used by the move, without extension
    movie_type: one of traj|matrix
                - traj: records the whole trajectory
                - matrix: a whole revolution around the z axis
    duration: duration in seconds of the final movie (default=5)
    ray: set to 'draw' to disable ray tracing (default='ray')
    """
    duration = int(duration)

#    # store settings for later restoration
#    store_settings.load()
#    cmd.sync()
#    settings_file = '/tmp/mt_settings.py'
#    cmd.do(f'store_settings {settings_file}')

    config.rendering()
    cmd.sync()

    cmd.mview('reset')
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
        cmd.mview('store', frames)
        # TODO: moving objects seems a bit hacky. Any better solution?
        for obj in cmd.get_object_list():
            cmd.mview('store', object=obj)
            cmd.sync()
        for i in range(4):
            frame = (frames * (i+1)/4)
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

    if savefile:
        movie.produce(savefile, mode=mode)
    else:
        cmd.mplay()

    # restore previous settings
#    cmd.do(f'run {settings_file}')


def load():
    cmd.extend('mt_movie', mt_movie)


if __name__ == 'pymol':
    load()
