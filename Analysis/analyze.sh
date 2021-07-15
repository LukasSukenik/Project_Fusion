#!/bin/bash

#
# Run this script in storage folder of simulation
#

#
# Base on naming scheme of used precycle script
#
xtcname=res.0001.xtc
logname=res.0001.log

#
# Extract lammps thermo_style output from *.log file
#
function gen_thermo
{
egrep '[0-9]{2}[.][0-9]{5}' res.0001.log | head -n 50 > thermo
for name in `ls -v res.*.end`
do
  #
  # Check fusion state from structure file
  #
  state=`./check_fused.sh $name`
  if [[ $state != *"fused"* ]]
  then
    egrep '[0-9]{2}[.][0-9]{5}' "res.`echo $name | cut -d"." -f2`.log" >> thermo
    break
  fi
done
}

#
# Get time of fusion (past stalk phase) from thermo output
#
# columns: step temp pe ke c_pair_lig_rec c_pair_ves1_rec c_pair_ves2_rec c_gyr c_gyr1 c_gyr2 c_cm1[1]
#
function get_fusion_time
{
  #
  # Because of performance we read end to start.
  #
  tac thermo | awk '
    {
  
    }
    END{
      print fusion_time
    }
  '
}



#
# Get time when liposome hydrophobic content  first mixes
#
function get_hydrophob_contact_time
{
  awk '
    END{
      print contact_time
    }    
  ' thermo
}



#
# Specify folders, default all
#
if [ "$1" = "" ]; then
  inn="r*"
else
  inn=$1
fi

echo "Selected folders: " \`ls -d ${inn}\`
echo "Selected folders: " `ls -d ${inn}` 1>&2

for i in `ls -d ${inn}`
do
  cd $i/storage
    cp ../.././check_fused.sh .    
    gen_thermo
    fus_time=`get_fusion_time `
    hydrofob_contact=``
    echo "$i $hydrofob_contact $fus_time"
    echo -n "$i, " 
  cd ../..
done
