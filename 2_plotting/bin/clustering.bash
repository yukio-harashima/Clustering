#!/bin/bash
#-----
echo "cluster_unified.py start (step 1/4)"
cluster_unified.py
echo "cluster_unified.py done"

echo "gensnap2_KDE.py and clmeca.py  start (step 2/4)"
gensnap2_KDE.py
clmeca.py
echo "gensnap2_KDE.py and clmeca.py  done"

echo "ffi_MMcl4.bash  start (step 3/4)"
ffi_MMcl4.bash

echo "ffi_MMcl4.bash  done"

# echo "ffi_PMcl4.bash start (step 4/4)"
# ffi_PMcl4.bash
# echo "ffi_PMcl4.bash done"

echo "all step done! check result."