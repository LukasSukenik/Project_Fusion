#!/bin/bash

# Folder naming conventions: r$(radius)_l$(lipid_count)
function create
{
lipid_count_list=$2
radius_list=$3

for rad in $radius_list
do 
  for lcount in $lipid_count_list
  do
    folder="r${rad}_l${lcount}"
    echo $folder

    rcount=$(($lcount / 2))
    
    mkdir $folder
    
    cp prototype/* $folder/

    cd $folder
      ./prep.sh $1 $lcount $rcount $rad $4
    cd ..
  done  
done
}

r_run=100000 #steps should work fine
lipid_count=1000
radius=10
box=30

create $r_run "$lipid_count" "$radius" $box


