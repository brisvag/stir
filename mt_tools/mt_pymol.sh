#!/usr/bin/env sh

# shell wrapper for mt_pymol.py

dir_path="`dirname \"$0\"`"

if command -v python3 > /dev/null 2>&1; then
    $dir_path/mt_pymol.py $@
else
    echo "You don't seem to have python3 installed!"
    echo "You need python3 to run mt_pymol."
fi
