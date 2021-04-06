#!/bin/bash

# Folder naming conventions: r$(radius)_l$(lipid_count)
function run
{
lipid_count_list=$2
radius_list=$3

for rad in $radius_list
do 
  for lcount in $lipid_count_list
  do
    folder="r${rad}_l${lcount}"
    echo $folder

    cd $folder
      psubmit que precycle_mdrun_lammps options
    cd ..
  done  
done
}

r_run=100000 #steps should work fine
lipid_count=1000
radius=10

run $r_run "$lipid_count" "$radius"


