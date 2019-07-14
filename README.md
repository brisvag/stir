# Simple Trajectory Inspection and Representation

**STIR** is a wrapper for [PyMOL](https://github.com/schrodinger/pymol-open-source) that provides a
collection of tools for the visualization of trajectories.

Its main focus are [Martini coarse-grained systems](https://cgmartini.nl) 

# Installation

You will first need to install [garnish](https://github.com/mjtadema/garnish), since handling of git dependencies
is apparently a mess.

Once you've done that, you can just use pip:
```
pip install git+git://github.com/brisvag/stir.git#egg=stir
```
pip will automatically add the `stir` command to your `PATH`.

# Usage

Use `stir` to open pymol and automatically run all the default tools. Check out the help with:
```
stir -h
```

# Tools

## stir
Wrapper for pymol. Loads a system and provides access to the other tools.
```
stir system.gro topol.top md.xtc
```
Initializes with the repo's `pymolrc` and automatically runs `garnish` and `mt_nice`. 
Also loads the other tools, ready to be used with `mt_toolname` from within pymol. 

---

## garnish
Draws bonds and elastic network for coarse-grained systems.
```
garnish selection [, tpr_file|top_file]
```
Source: https://github.com/mjtadema/garnish

## mt_nice
Provides a series of function for easy selection and visualization.

#### mt_nice
Wraps other subtools to select, color and show molecules nicely.
```
mt_nice [style [, selection]]
```

#### mt_sele
To automatically create (or delete) commonly used selections, such as `lip` for lipids and `BB` for backbone beads:
```
mt_sele [delete]
```

#### mt_color 
Color all atoms in the selection based on a common identification (chain id, residue id...).
```
mt_color resi|chain|name|... [, selection]
```

## mt_supercell
Shows periodic images. To show a 3x3 grid in the x,y plane:
```
mt_supercell 3,3,1
```
To show first 2 neighbouring cells in the z direction:
```
mt_supercell 1,1,5
```
WARNING: Pymol creates actual copies of the system, so this is an expensive command!

Source: https://github.com/speleo3/pymol-psico/blob/master/psico/xtal.py

## mt_movie
Makes a nice ray-traced movie of the trajectory. A few movie templates are available (see `help mt_movie`). Try:
```
mt_movie bullettime, test_movie, width=640, height=480
```
