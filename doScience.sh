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
if [ -d $HERE/Circular ]; then
	for wall in Circular; do
		cd $wall
		if [ -f $wall/noise1.0 ]; then
			for noise in noise*; do
				cd $noise
				for pred in pred*; do
					cd $pred
					for ni in group[0-9]*; do
						cd $ni
						echo $(pwd)
						python -b $DIR/postProcessing/getStatistics.py
						python -b $DIR/postProcessing/getAlphaBeta.py
						cd ..
					done
					python -b $DIR/plotting/plotting/plotStats.py
					python -b $DIR/plotting/plotting/plotAlpha.py
					python -b $DIR/plotting/plotting/plotBeta.py
					python -b $DIR/plotting/plotting/plotGroupPredation.py
					cd ..
				done
				cd ..
			done
		fi
		for pred in pred*; do
			cd $pred
			for ni in group[0-9]*; do
				cd $ni
				echo $(pwd)
				python -b $DIR/postProcessing/getStatistics.py
				python -b $DIR/postProcessing/getAlphaBeta.py
				cd ..
			done
			python -b $DIR/plotting/plotStats.py
			python -b $DIR/plotting/plotAlpha.py
			python -b $DIR/plotting/plotBeta.py
			python -b $DIR/plotting/plotGroupPredation.py
			cd ..
		done
		cd $HERE
	done
fi


for wall in Square* ; do
	cd $wall
	if [ -d noise1.0 ]; then
		for noise in noise*; do
			cd $noise
			for pred in pred*; do
				cd $pred
				for ni in group[0-9]*; do
					cd $ni
					echo $(pwd)
					if [ -f data-*-1.h5 ]; then
						python -b $DIR/postProcessing/getStatistics.py
						python -b $DIR/postProcessing/getAlphaBeta.py
					fi
					cd ..
				done
				python -b $DIR/plotting/plotGroupPredation.py
				cd ..
			done
			if [ -f $noise/$pred/$ni/output/beta]; then
				python -b $DIR/plotting/plotStats.py
				python -b $DIR/plotting/plotAlpha.py
				python -b $DIR/plotting/plotBeta.py
			fi
			cd ..
		done
	else
	for pred in pred*; do
		cd $pred
		for ni in group[0-9]*; do
			cd $ni
			echo $(pwd)
			if [ -f data-*-1.h5 ]; then
				python -b $DIR/postProcessing/getStatistics.py
				python -b $DIR/postProcessing/getAlphaBeta.py
			fi
			cd ..
		done
		python -b $DIR/plotting/plotGroupPredation.py
		cd ..
	done
	python -b $DIR/plotting/plotStats.py
	python -b $DIR/plotting/plotAlpha.py
	python -b $DIR/plotting/plotBeta.py
	cd $HERE
fi
done
