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

rm -f *.h5
rm -f *.pyc

python -b ~/CSM/CSM.py $1

if [ "$1" = "1" ]; then
	mkdir -p frames
	rm -f frames/*
	python $DIR/plotFrames.py
	p=`pwd`
	vid=`echo ${p##*/}`
	ffmpeg -framerate 10 -pattern_type glob -i 'frames/*.png' -vcodec mpeg4 -qscale:v 2 frames/${vid}.mp4

	MAC=$INSTALLDIR_HOME
	if [ -f $MAC/params.py.dist ]; then
	  open frames/${vid}.mp4
	fi

	mkdir -p plots
	mkdir -p output

	python $DIR/plotSheepAcceleration.py
	python $DIR/plotSheepVelocity.py
	python $DIR/plotDogAcceleration.py
	python $DIR/plotDogVelocity.py
	python $DIR/plotTimeStep.py
	python $DIR/plotPredation.py

	if [ -f $MAC/params.py.dist ]; then
	  terminal-notifier -title "CSM" -message "Done with simulation" -subtitle ${vid} -sound 'glass'
	  open frames/${vid}.mp4
	fi

fi
