#!/bin/sh
INSTALLDIR=/home/b1033128/Documents/CSM
INSTALLDIR_TOPSY=/home/b1033128/CSM
INSTALLDIR_HOME=/Users/hayleymoore/Documents/Uni/Documents/CSM

if [ -f $INSTALLDIR/params.py.dist ]; then
	DIR=$INSTALLDIR
fi
if [ -f $INSTALLDIR_TOPSY/params.py.dist ]; then
	DIR=$INSTALLDIR_TOPSY
	module load python/2.7.9
	module load hdf5
	module load topsy-openmpi 
fi
if [ -f $INSTALLDIR_HOME/params.py.dist ]; then
	DIR=$INSTALLDIR_HOME
fi

HERE=$(pwd)
for wall in Square Circular; do
	cd $wall
	for pred in pred*; do
		cd $pred
		for ni in group*; do
			cd $ni
			echo $(pwd)
			python -b $DIR/getStatistics.py
			python -b $DIR/getBeta.py
			cd ..
		done
		cd ..
	done
	cd $HERE
done
