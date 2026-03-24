#!/bin/bash

clno=5
mecasize=1.5
# info=($(gmt gmtinfo clmeca_ave.dat -C))

# echo `info`

# xmin=$(echo "${info[2]}-0.2" | bc); xmax=$(echo "${info[3]}+0.2" | bc)
# ymin=$(echo "${info[0]}-0.1" | bc); ymax=$(echo "${info[1]}+0.1" | bc)
# xmin=$(echo "scale=2;${xmin}/1." | bc); xmax=$(echo "scale=2;${xmax}/1." | bc)
# ymin=$(echo "scale=2;${ymin}/1." | bc); ymax=$(echo "scale=2;${ymax}/1." | bc)
r=139.7/140.3/35.7/36.3

gmt begin meca_sample${clno} png
    # gmt makecpt -Chot -T0.5/6.5/1 -I --COLOR_NAN=white --COLOR_BACKGROUND=white --COLOR_FOREGROUND=white
    gmt makecpt -Cvik -T0.5/5.5/1 -I --COLOR_NAN=white --COLOR_BACKGROUND=white --COLOR_FOREGROUND=white
    gmt basemap -JM12 -R${r} -BWSNE
    # gmt coast -Df -W0.25 -LJBR+jTR+o0/1+c20+w500+f
    # gmt meca "meca.txt" -Sa3.5 -Gred
    awk -F'\t' '{if ( $1=='${clno}') print 140, 36,$1,$3,$4,$5,$6,$7,$8,22}' clave_mod.dat | gmt meca -Sm${mecasize}c+l -T0 -C
gmt end show
