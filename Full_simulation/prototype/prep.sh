#!/bin/bash

# input: number_of_steps interaction_strength nano_size nano_ratio ligand_count

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



function gen_ff
{
  interaction=$3

  cp $1 $2
  sed -i -e 's/sed_interaction/'"$interaction"'/g' $2
}

function gen_file
{
  boxx=$3

  cp $1 $2
  sed -i -e 's/sed_box/'"$boxx"'/g' $2
}


run_steps=$1
strength=$2
box=$5

gen_file load_equi_nano_prescript load_equi_nano $box
gen_file load_ves2_prescript load_ves2 $box
./ico load_equi_nano load_ves2 > data.start
gen_in "in.prescript_production" "in.production" $run_steps
gen_in "in.prescript_restart" "in.restart" $run_steps
gen_ff "force_field_prescript" "force_field" $strength

sed -i -e 's/mass_1/'"1"'/g' data.start
sed -i -e 's/mass_2/'"1"'/g' data.start
sed -i -e 's/mass_3/'"1"'/g' data.start
sed -i -e 's/mass_4/'"1"'/g' data.start
sed -i -e 's/mass_5/'"1"'/g' data.start


