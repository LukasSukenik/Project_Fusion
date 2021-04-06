#!/bin/bash

# Folder naming conventions: r${radius}_i${interaction}_c${prolatenes}_s${scale}
function create
{
run=$1
rad=$2
inter=$3
prol=$4
sca=$5
boxx=$6

folder="r${rad}_i${inter}_c${prol}_s${sca}"

echo $folder

mkdir $folder

cp prototype/* $folder/

# Copy endfile (data.end) from vesicle equilibration and rename to data.equi

cd $folder
  ./prep.sh $run ${inter} ${prol} ${sca} ${boxx}
cd ..

}

r_run=50000 # 5 M steps seems reasonable - will finish in ~8h on 4 cores
radius=10
interaction=3
scale=2
prolatenes=7 # CANT use 1.0 -> 1.001
box=45

create $r_run $radius $interaction $prolatenes $scale $box



