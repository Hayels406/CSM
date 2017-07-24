#!/bin/sh
INSTALLDIR=/home/b1033128/Documents/CSM
INSTALLDIR_TOPSY=/home/b1033128/CSM
INSTALLDIR_HOME=/Users/hayleymoore/Documents/Uni/Documents/CSM


rm -f *.h5
python CSM.py
mkdir -p frames
rm -f frames/*
if [ -f $INSTALLDIR/params.py.dist ]; then
  python $INSTALLDIR/plotFrames.py
  p=`pwd`
  vid=`echo ${p##*/}`
  ffmpeg -framerate 5 -pattern_type glob -i 'frames/*.png' -vcodec mpeg4 -qscale:v 2 frames/${vid}.mp4
  mkdir -p plots
  python $INSTALLDIR/plotSheepAcceleration.py
  python $INSTALLDIR/plotSheepVelocity.py
  python $INSTALLDIR/plotDogAcceleration.py
  python $INSTALLDIR/plotDogVelocity.py
  python $INSTALLDIR/plotTimeStep.py
  python $INSTALLDIR/convexHull.py

fi

if [ -f $INSTALLDIR_HOME/params.py.dist ]; then
  python $INSTALLDIR_HOME/plotFrames.py
  p=`pwd`
  vid=`echo ${p##*/}`
  ffmpeg -framerate 5 -pattern_type glob -i 'frames/*.png' -vcodec mpeg4 -qscale:v 2 frames/${vid}.mp4
  mkdir -p plots
  python $INSTALLDIR_HOME/plotSheepAcceleration.py
  python $INSTALLDIR_HOME/plotSheepVelocity.py
  python $INSTALLDIR_HOME/plotDogAcceleration.py
  python $INSTALLDIR_HOME/plotDogVelocity.py
  python $INSTALLDIR_HOME/plotTimeStep.py
  python $INSTALLDIR_HOME/convexHull.py

fi

if [ -f $INSTALLDIR_TOPSY/params.py.dist ]; then
  python $INSTALLDIR_TOPSY/plotFrames.py
  p=`pwd`
  vid=`echo ${p##*/}`
  ffmpeg -framerate 5 -pattern_type glob -i 'frames/*.png' -vcodec mpeg4 -qscale:v 2 frames/${vid}.mp4
  mkdir -p plots
  python $INSTALLDIR_TOPSY/plotSheepAcceleration.py
  python $INSTALLDIR_TOPSY/plotSheepVelocity.py
  python $INSTALLDIR_TOPSY/plotDogAcceleration.py
  python $INSTALLDIR_TOPSY/plotDogVelocity.py
  python $INSTALLDIR_TOPSY/plotTimeStep.py
  python $INSTALLDIR_TOPSY/convexHull.py

fi
