#!/bin/bash
echo 'Running SNR script'
rm detangles.txt
touch detangles.txt
detAnglefile=detangles.txt
a=$(asy -V spherelite.asy)
thetadet=${a[0]}
phidet=${a[1]}
psidet=${a[2]}
echo "$thetadet" >> "$detAnglefile"
echo "$phidet" >> "$detAnglefile"
echo "$psidet" >> "$detAnglefile"
sed '/^$/d' $detAnglefile > $detAnglefile.output # remove possible blank lines
mv $detAnglefile.output $detAnglefile
python plot.py
