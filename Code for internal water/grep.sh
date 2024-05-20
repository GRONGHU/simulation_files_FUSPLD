#!/bin/bash  
dtime=20
rm internal_water_resid.npy
chmod 777 EDTSurf 
for ((i=0;i<=100000;i=i+1))
do
 rm mdtest_nowion.pdb
 rm mdtest_all.gro
 rm SURFACE_OUTPUT-cav.pdb
 rm SURFACE_OUTPUT.asa
 rm SURFACE_OUTPUT.ply
 let j=dtime*i
 echo $j
 echo 1 1 | gmx trjconv -f md1.xtc -s md1.tpr -o mdtest_nowion.pdb -pbc mol -center -b ${j} -e ${j}
 echo 1 0 | gmx trjconv -f md1.xtc -s md1.tpr -o mdtest_all.gro -pbc mol -center -b ${j} -e ${j}
 ./EDTSurf -i mdtest_nowion.pdb -s 2 -h 2 -p 1.4 -o SURFACE_OUTPUT
 python select_water.py

done
