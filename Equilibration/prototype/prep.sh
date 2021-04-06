#!/bin/bash

# input: number_of_steps lipid_count receptor_count radius

function gen_in
{
  r_run=$3
  r_damp=100.0
  r_timestep=0.01

  cp $1 $2
  sed -i -e 's/r_damp/'"$r_damp"'/g' $2
  sed -i -e 's/r_random/'"$RANDOM"'/g' $2
  sed -i -e 's/r_timestep/'"$r_timestep"'/g' $2
  sed -i -e 's/r_run/'"$r_run"'/g' $2
}

function gen_ves
{
  lipid_count=$3
  receptor_count=$4 
  radius=$5

  cp $1 $2
  sed -i -e 's/sed_lipids/'"$lipid_count"'/g' $2
  sed -i -e 's/sed_receptors/'"$receptor_count"'/g' $2
  sed -i -e 's/sed_radius/'"$radius"'/g' $2
}

gen_in "in.prescript_equi" "in.equi" $1
gen_ves prescript_vesicle vesicle $2 $3 $4
./gen_membrane vesicle


