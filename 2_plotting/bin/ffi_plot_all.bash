#!/bin/bash
#-----
ffi_solution.bash        >  .tmpout
echo ffi_solution.png
ffi_MapMeca.bash         >  .tmpout
echo ffi_MapMeca.png
ffi_MapX.bash p          >  .tmpout
echo ffi_MapX.pdf
ffi_RupStrDip.bash       >  .tmpout
echo ffi_RupStrDip.pdf
ffi_RSD.bash             >  .tmpout
echo ffi_RSD.png
ffi_NDC.bash             >  .tmpout
echo ffi_NDC.pdf
rm -rf .tmpout
#-----
ffi_MapMecaSnap.bash     >  .tmpout
echo ffi_MapMecaSnap.png
ffi_PlaneMecaSnap.bash   >  .tmpout
echo ffi_PlaneMecaSnap.pdf
ffi_MapXsnap.bash        >  .tmpout
echo ffi_MapXsnap.pdf
#-----
ffi_WaveComp.bash        >  .tmpout
echo ffi_WaveComp.png
ffi_AziEquiStaMap.bash   >  .tmpout
echo ffi_AziEquiStaMap.pdf
rm -rf .tmpout
