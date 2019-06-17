#!/usr/bin/env sh

# shell wrapper for mt_pymol.py

if command -v python3 > /dev/null 2>&1; then
    mt_pymol.py $@
else
    echo "You don't seem to have python3 installed!"
    echo "You need python3 to run mt_pymol."
fi
