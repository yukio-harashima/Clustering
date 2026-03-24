#!/bin/bash
#-----
ffi_solution_v2.bash        >  .tmpout
echo ffi_solution.png
ffi_MapMeca.bash         >  .tmpout
echo ffi_MapMeca.png
ffi_RSD.bash             >  .tmpout
echo ffi_RSD.png
ffi_WaveComp.bash        >  .tmpout
echo ffi_WaveComp.png
rm -rf .tmpout
