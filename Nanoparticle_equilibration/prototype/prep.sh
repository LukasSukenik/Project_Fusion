#!/bin/bash

# input: number_of_steps interaction_strength nano_size nano_ratio ligand_count

function gen_in
{
  r_run=$3
  r_damp=100.0
  r_timestep=0.01
  r_radius=$(($4-3))

  cp $1 $2
  sed -i -e 's/r_damp/'"$r_damp"'/g' $2
  sed -i -e 's/r_random/'"$RANDOM"'/g' $2
  sed -i -e 's/r_timestep/'"$r_timestep"'/g' $2
  sed -i -e 's/r_run/'"$r_run"'/g' $2
  sed -i -e 's/r_radius/'"$r_radius"'/g' $2
}

function gen_nano
{
  rscale=$3
  c=$4
  lig_dens=$5
  nano_dens=$6
  boxx=$7
  
  S_sphere=`echo $rscale | awk '{pi=atan2(0, -1); print 4*pi*$1*$1;}'`
  equi_rad=$rscale # equatorial radius
  polar_rad=`echo $rscale $c | awk '{print $1*$2}'` 
  S_actual=`echo $equi_rad $polar_rad | awk 'function asin(x) { return atan2(x, sqrt(1-x*x)) }
{pi=atan2(0, -1); a=$1; c=$2;
print 2*pi*a*a+((2*pi*a*c*c)/sqrt(c*c-a*a))*asin((sqrt(c*c-a*a))/(c))}'`

  if [[ $c -eq 1 ]]; then
    echo "Prolate/oblate ratio c can't be 1, use 1.001"
    exit
  else
    nano=`echo $S_actual $nano_dens | awk '{print $1*$2}'`
    lig=`echo $S_actual $lig_dens | awk '{print $1*$2}'`
  fi
  
  pos_z=`grep -B 99999 Velocities data.equi | grep -A 99999 Atoms | awk 'BEGIN{max=-999}{if(max<$7 && length($7) != 0) max=$7} END{print max}'`
  pos_z=`echo $pos_z $equi_rad | awk '{print ($1 + $2 +1) }'`

  cp $1 $2
  sed -i -e 's/sed_nano/'"$nano"'/g' $2
  sed -i -e 's/sed_scale/'"$rscale"'/g' $2
  sed -i -e 's/sed_c/'"$c"'/g' $2
  sed -i -e 's/sed_lig/'"$lig"'/g' $2
  sed -i -e 's/sed_z/'"$pos_z"'/g' $2
  sed -i -e 's/sed_box/'"$boxx"'/g' $2
}

function gen_ves
{
  boxx=$3

  cp $1 $2
  sed -i -e 's/sed_box/'"$boxx"'/g' $2
}

function gen_ff
{
  interaction=$3

  cp $1 $2
  sed -i -e 's/sed_interaction/'"$interaction"'/g' $2
}

# constants
num_lig_area=1 # number of nanoparticle ligands per area
num_nano_area=3

# variables $run ${inter} ${prol} ${sca}
run_steps=$1
strength=$2
c_param=$3
scale=$4
box=$5

gen_ves load_ves_prescript load_ves $box
gen_nano "nanoparticle_prescript" "nanoparticle" $scale $c_param $num_lig_area $num_nano_area $box
./ico load_ves nanoparticle > data.ves_nano
gen_in "in.prescript_production" "in.production" $run_steps $box
gen_ff "force_field_prescript" "force_field" $strength

sed -i -e 's/mass_1/'"1"'/g' data.ves_nano
sed -i -e 's/mass_2/'"1"'/g' data.ves_nano
sed -i -e 's/mass_3/'"1"'/g' data.ves_nano
sed -i -e 's/mass_4/'"1"'/g' data.ves_nano
sed -i -e 's/mass_5/'"1"'/g' data.ves_nano
sed -i -e 's/mass_6/'"1"'/g' data.ves_nano

