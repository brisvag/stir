# Martini Pymol Tools

This repository is a collection of tools that will make your life easier
when using pymol to visualize Martini coarse-grained trajectories.

# Installation

Just clone the repo:
```
pip install git+git://github.com/martinitoolspymol/mtools.git#egg=martinitoolspymol
```

# Usage

To open pymol and automatically run all the default tools, run:
```
python -m mtools
```

To make it easier to launch it from anywhere, add an alias in your `~/.bashrc`:
```
alias mt_pymol="python -m mtools"
```

Alternatively, you can run the individual tools from pymol (or add them to you `~/.pymolrc.pml`)
and then use their functions from within pymol as normal.

# Tools

## mt_pymol
Wrapper for pymol and most of the other tools.
```
mt_pymol system.gro topol.top md.xtc
```
Initializes with the repo's `pymolrc` and automatically runs `garnish` and `mt_sele`. 
Also loads `mt_supercell` and `mt_movie` to be used when needed.

## garnish
Draws bonds and elastic network for coarse-grained systems.
```
garnish selection [, tpr_file|top_file]
```
Source: https://github.com/mjtadema/garnish

## mt_nice
Provides a series of function for easy selection and visualization.

### mt_nice
Wraps other subtools to select, colors and shows molecules nicely.
```
mt_nice [style [, selection]]
```

### mt_sele
To automatically create (or delete) commonly used selections, such as `lip` for lipids and `BB` for backbone beads:
```
mt_sele [delete]
```

### mt_color 
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
WARNING: Pymol creates aactual copies of the system, so this is an expensive command!

Source: https://github.com/speleo3/pymol-psico/blob/master/psico/xtal.py

## mt_movie
Makes a nice ray-traced movie of the trajectory. A few movie templates are available (see `help mt_movie`). Try:
```
mt_movie bullettime, test_movie, width=640, height=480
```
