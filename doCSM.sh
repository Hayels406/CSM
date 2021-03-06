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

rm -f *-$1.h5
rm -f *.pyc

python -b $DIR/CSM.py $1

if [ "$1" = "1" ]; then
	echo 'Starting frames'
	mkdir -p frames
	rm -f frames/*
	python $DIR/plotting/plotFrames.py
	p=`pwd`
	vid=`echo ${p##*/}`
	ffmpeg -framerate 10 -pattern_type glob -i 'frames/*.png' -vcodec mpeg4 -qscale:v 2 frames/${vid}.mp4

	MAC=$INSTALLDIR_HOME
	if [ -f $MAC/params.py.dist ]; then
	  open frames/${vid}.mp4
	fi

	mkdir -p plots
	mkdir -p output
	echo 'Starting plots'
	python $DIR/plotting/plotSheepAcceleration.py
	python $DIR/plotting/plotSheepVelocity.py
	python $DIR/plotting/plotDogAcceleration.py
	python $DIR/plotting/plotDogVelocity.py

	if [ -f $MAC/params.py.dist ]; then
	  terminal-notifier -title "CSM" -message "Done with simulation" -subtitle ${vid} -sound 'glass'
	  open frames/${vid}.mp4
	fi

fi
echo 'Done'
