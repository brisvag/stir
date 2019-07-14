# Simple Trajectory InspectoR

-- "_**Stirred**, not shaken._"

**STIR** is a wrapper for [PyMOL](https://github.com/schrodinger/pymol-open-source) that provides a
collection of tools for the visualization of [gromacs](www.gromacs.org) molecular dynamics trajectories.

Its main focus are [Martini coarse-grained systems](www.cgmartini.nl).

# Installation

Grab the latest PyMOL from [here](https://github.com/schrodinger/pymol-open-source). You will also
have to install [garnish](https://github.com/mjtadema/garnish).

Once you've done that, you can just use pip:
```bash
pip install git+git://github.com/brisvag/stir.git#egg=stir
```
pip will automatically add the `stir` command to your `PATH`.

# How to `stir` your Martini

Use `stir` to open PyMOL and automatically run all the default tools. Check out the help with:
```bash
stir -h
```

For example, to load a typical system:
```bash
stir system.gro topol.top md.xtc
```
`stir` initializes PyMOL with its own `pymolrc` and automatically runs `garnish` and `nice`.
It also loads the other tools, ready to be used from within PyMOL. 

---

# Tools

## garnish
Draws bonds and elastic network for coarse-grained systems.
```
garnish selection [, tpr_file|top_file]
```
Source: https://github.com/mjtadema/garnish

## nice
Provides a series of functions for easy selection and visualization.

#### nice
Wraps other subtools to select, color and show molecules nicely.
```
nice [style [, selection]]
```

#### nicesele
To automatically create (or delete) commonly used selections, such as `lip` for lipids and `BB` for backbone beads:
```
nicesele [delete]
```

#### nicecolor 
Color all atoms in the selection with the same random color based on a common identifier (chain id, residue name...).
To see a list of all the options, check out [PyMOL's iterate command](https://pymolwiki.org/index.php/Iterate).
```
nicecolor resi|chain|name|... [, selection]
```

## supercell
Shows periodic images. To show a 3x3 grid in the x,y plane:
```
supercell 3,3,1
```
To show first 2 neighbouring cells in the z direction:
```
supercell 1,1,5
```
Source: https://github.com/speleo3/pymol-psico/blob/master/psico/xtal.py

## cheese
Makes a nice, ray-traced picture or movie of the system. A few options are available (see `help cheese`). Try:
```
cheese snap
```
to take a quick ray-traced picture. For a movie, try:
```
cheese bullettime, test_movie, width=640, height=480
```
