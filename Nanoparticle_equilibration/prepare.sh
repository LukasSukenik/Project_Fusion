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
head=$7

folder="r${rad}_h${head}_i${inter}_c${prol}_s${sca}"

echo $folder

mkdir $folder

cp prototype/* $folder/

# Copy endfile (data.end) from vesicle equilibration and rename to data.equi
cp ../Equilibration/r${rad}\_h${head}\_*/data.end $folder/data.equi

cd $folder
  ./prep.sh $run ${inter} ${prol} ${sca} ${boxx} ${head}
cd ..

}

r_run=500000 
radius=10
interaction=2
prolatenes=3 # CANT use 1.0 -> 1.001
scale=3
box=30
head_size=0.95

create $r_run $radius $interaction $prolatenes $scale $box $head_size



