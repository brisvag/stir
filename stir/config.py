# Copyright 2019-2019 the stir authors. See copying.md for legal info.

"""
configuration functions for pymol containing a collection of settings for a specific purpose
"""

from pymol import cmd, util
import psutil


def pymolrc():
    """
    startup settings, always applicable
    """
    # multithreading
    cmd.set('max_threads', psutil.cpu_count())

    # color space
    cmd.space('pymol')

    # environment and color
    cmd.bg_color('white')
    cmd.set('depth_cue', 'off')
    cmd.set('fog', 'off')

    # nicer visualization
    cmd.set('valence', 'on')
    cmd.set('cartoon_side_chain_helper', 'on')
    cmd.set('cartoon_highlight_color', 'grey50')
    cmd.set('sphere_mode', '9')

    # make sure trajectories are only loaded once defer_builds_mode is 3
    def load_traj(*args, **kwargs):
        cmd.set('defer_builds_mode', 3)
        cmd.sync()
        cmd.load_traj(*args, **kwargs)
    cmd.extend('load_traj', load_traj)


def trajectory():
    """
    optimization settings for trajectory visualization
    """
    # improves performance drastically when loading a lot of states
    cmd.set('defer_builds_mode', 3)


def rendering():
    """
    optimal settings for movie (or picture) rendering
    """

    # ray trace movie options
    cmd.set('ray_trace_frames', 'off')  # can be turned on at conversion
    cmd.set('movie_loop', 'on')
    cmd.set('opaque_background', 'off')     # useless for movies, but not figures
    cmd.set('movie_auto_interpolate', 'off')
    cmd.set('movie_auto_store', 'off')

    # crank up quality and down performance
    util.performance(0)
    cmd.set('hash_max', 1000)
    cmd.set('line_smooth', 'on')
    cmd.set('antialias', 2)

    # lighting both in- and outside of a surface
    cmd.set('two_sided_lighting', 'on')

    # shadows and light options
    cmd.set('ray_shadow', 'on')
    cmd.set('light_count', 4)

    # TODO: hard to optimize these settings. Need more testing
    cmd.set('ambient', 0.1)
    cmd.set('direct', 0.1)
    cmd.set('reflect', 1.3)
    cmd.set('shininess', 10.0)
    cmd.set('specular', 0.3)
    cmd.set('spec_count', 1)
    cmd.set('power', 1.0)
    cmd.set('specular_intensity', 0.1)
    cmd.set('spec_direct', 0.2)

    cmd.set('ray_shadow_decay_factor', 0.10000)
    cmd.set('ray_shadow_decay_range', 2.00000)

    # TODO: ambient occlusion is better than shadows when rendering surfaces
    #cmd.set('ambient_occlusion_mode', 2)
    #cmd.set('ambient_occlusion_scale', 13)
    #cmd.set('ambient_occlusion_smooth', 20)

    cmd.set('surface_quality', 2)
    cmd.set('sphere_quality', 4)

    cmd.set('movie_panel', 'on')

    # multiplayer transparency
    cmd.set('transparency_mode', 1)

    # labels
    cmd.set('label_size', 30)
    cmd.set('label_color', 'black')
    cmd.set('label_outline_color', 'white')
    cmd.set('label_font_id', 7)
    cmd.set('label_shadow_mode', 1)


if __name__ == 'pymol':
    cmd.extend('stir_config_pymolrc', pymolrc)
    cmd.extend('stir_config_trajectory', trajectory)
    cmd.extend('stir_config_rendering', rendering)
