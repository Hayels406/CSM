#!/bin/sh
INSTALLDIR=/home/b1033128/Documents/CSM
INSTALLDIR_TOPSY=/home/b1033128/CSM
INSTALLDIR_HOME=/Users/hayleymoore/Documents/Uni/Documents/CSM

if [ -f $INSTALLDIR/params.py.dist ]; then
	DIR=$INSTALLDIR
fi
if [ -f $INSTALLDIR_TOPSY/params.py.dist ]; then
	DIR=$INSTALLDIR_TOPSY
fi
if [ -f $INSTALLDIR_HOME/params.py.dist ]; then
	DIR=$INSTALLDIR_HOME
fi

HERE=$(pwd)
for wall in Square* Circular*; do
	cd $wall
	for noise in noise*; do
		cd noise
		for pred in pred*; do
			cd $pred
			for ni in group[0-9]*; do
				cd $ni
				echo $(pwd)
				#python -b $DIR/getStatistics.py
				#python -b $DIR/getAlphaBeta.py
				cd ..
			done
			#python -b $DIR/plotStats.py
			#python -b $DIR/plotAlpha.py
			#python -b $DIR/plotBeta.py
			#python -b $DIR/plotGroupPredation.py
			cd ..
		done
		cd ..
	done
	for pred in pred*; do
		cd $pred
		for ni in group[0-9]*; do
			cd $ni
			echo $(pwd)
			#python -b $DIR/getStatistics.py
			#python -b $DIR/getAlphaBeta.py
			cd ..
		done
		#python -b $DIR/plotStats.py
		#python -b $DIR/plotAlpha.py
		#python -b $DIR/plotBeta.py
		#python -b $DIR/plotGroupPredation.py
		cd ..
	done
	cd $HERE
done
