# Martini Pymol Tools

This repository is a collection of tools that will make your life easier
when using pymol to visualize Martini coarse-grained trajectories.

# Installation

To download the repo, use:
```
git clone --recurse-submodules
```

# Usage

To open pymol and automatically run all the default tools, run:
```
./mt_tools/mt_pymol.py
```

To make it easier to launch it from anywhere, add an alias in your `~/.bashrc`:
```
alias mt_pymol="/path/to/mt_pymol.py"
```

Alternatively, you can run the individual tools from pymol (or add them to you `~/.pymolrc.pml`)
and then use their functions from within pymol as normal.

# Tools

## mt_pymol
Wrapper for most of the other tools.
```
mt_pymol system.gro topol.top md.xtc
```

## pycg_bonds
Draws bonds and elastic network for coarse-grained systems.
```
cg_bonds selection [, tpr_file|top_file]
```
Source: [https://github.com/mjtadema/pycg_bonds]

## mt_sele
Automatically creates useful named selections, such as `lip` for lipids and `BB` for backbone beads.
```
mt_sele
mt_desele
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

Source: [https://github.com/speleo3/pymol-psico/blob/master/psico/xtal.py]
