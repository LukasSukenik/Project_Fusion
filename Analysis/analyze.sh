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
    echo -n "$i, " 
  cd ../..
done
