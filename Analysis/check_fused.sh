#!/bin/bash

input_file=$1

awk '
BEGIN{
  true=1
  false=0
  fused=false
}

{
  #
  # Determine fused vesicles from lammps data file
  #
}

END{
  if( fused == true)
  {
    print fused
  } 
  else
  {
    print meh 
  }
  
} ' $input_file
