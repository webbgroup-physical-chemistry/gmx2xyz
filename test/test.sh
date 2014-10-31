#! /bin/bash

cd ..
rm -v gmx2xyz
make
cd -

tndx=ndx.ndx
echo "[ test ]" > t.ndx
for i in `seq 1 150`; do echo $i >> t.ndx; done
find . -name "*.xyz" -print | xargs rm
if [ -z $1 ] ; then
time ../gmx2xyz -s tpr.tpr -f gro.gro -p top.top -n $tndx -a ~ritchie/software/tinker/params/CNC_amoebabio09.prm -e 15
else
for i in `seq 0 5`; do
find . -name "*.xyz" -print | xargs rm
time ../gmx2xyz -s tpr.tpr -f gro.gro -p top.top -n $tndx -a ~ritchie/software/tinker/params/CNC_amoebabio09.prm  &> j
rm j
done
fi
