from pymol import cmd
from pymol import util


# ray trace movie options
cmd.set('ray_trace_frames', 'on')
cmd.set('movie_loop', 'off')
cmd.set('opaque_background', 'off')
cmd.set('movie_auto_interpolate', 'off')

# crank up performance
util.performance(0)
cmd.set('hash_max', 1000)
cmd.set('line_smooth', 'on')

cmd.set('depth_cue', 'off')

# lighting both in- and outside of a surface
cmd.set('two_sided_lighting', 'on')

# shadows and light options
cmd.set('ray_shadow', 'on')
cmd.set('light_count', 4)

# TODO: hard to optimize these settings. Need more testing
cmd.set('ambient', 0.25000)
cmd.set('direct', 0.25000)
cmd.set('reflect', 0.80000)
cmd.set('shininess', 10.00000)
cmd.set('specular', 'off')
cmd.set('spec_count', 1)
cmd.set('power', 1.00000)
cmd.set('specular_intensity', 0.10000)
cmd.set('spec_direct', 0.10000)

cmd.set('ray_shadow_decay_factor', 0.10000)
cmd.set('ray_shadow_decay_range', 2.00000)

# ambient occlusion is better than shadows
cmd.set('ambient_occlusion_mode', 2)
cmd.set('ambient_occlusion_scale', 13)
cmd.set('ambient_occlusion_smooth', 20)

cmd.set('surface_quality', 2)
cmd.set('sphere_quality', 4)

cmd.set('movie_panel', 'on')
