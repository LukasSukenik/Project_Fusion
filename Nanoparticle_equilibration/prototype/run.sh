#!/usr/bin/bash

restart=0
mpirun -np 8 ./lmp_mpi -in in.production

./get_angle.sh

check=`awk '{print $1}' OUT`
echo $check
while [ "$check" != "DONE" ];do
	restart=$(($restart+1))
	echo $restart	
	mv data $restart.log
	sed  's/sed_rest/'"$restart"'/' in.restart >in.restartok
	mpirun -np 8 ./lmp_mpi -in in.restartok
	./get_angle.sh
	check=`awk '{print $1}' OUT`
	
done


