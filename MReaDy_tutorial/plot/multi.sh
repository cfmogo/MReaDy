#!/bin/bash
if [ "$1" == "" ]; then
      echo "Please give a number after multi.sh"
      exit 0
fi

#go to folder with data from MReaDy
cd ../bin/out/
pwd
#work on files
#filter energies to file geo*.dat
grep -in 'Step' geo$1.xyz > geo$1.dat
ls geo*.dat
read

echo "load '../../plot/time_changes.gnu'" | gnuplot
read

echo "load '../../plot/multi.gnu'" | gnuplot
read

 
cp *.eps ../../graphs/
cd -
