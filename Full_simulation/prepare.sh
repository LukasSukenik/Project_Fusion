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

folder="r${rad}_h${head}_4i${inter}_c${prol}_s${sca}_long"

echo $folder

mkdir $folder/

cp prototype/* $folder/

# Copy endfile (data.end) from vesicle equilibration and rename to data.equi
cp ../Nanoparticle_equilibration/$folder/data.equi $folder/
cp ../Nanoparticle_equilibration/$folder/data.end $folder/data.nano_equi

cd $folder
  ./prep.sh $run ${inter} ${prol} ${sca} ${boxx} ${head}
cd ..

}

r_run=5000000 # 5 M steps seems reasonable - will finish in ~8h on 4 cores
radius=10
interaction=2
prolatenes=4 # CANT use 1.0 -> 1.001
scale=2
box=55
head_size=0.95

create $r_run $radius $interaction $prolatenes $scale $box $head_size



