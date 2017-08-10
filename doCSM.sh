#!/bin/sh
rm -f *.h5
rm -f *.pyc

python -b CSM.py

mkdir -p frames
rm -f frames/*
python ~/Documents/Uni/Documents/CSM/plotFrames.py
p=`pwd`
vid=`echo ${p##*/}`
ffmpeg -framerate 10 -pattern_type glob -i 'frames/*.png' -vcodec mpeg4 -qscale:v 2 frames/${vid}.mp4

mkdir -p plots
mkdir -p output
python ~/Documents/Uni/Documents/CSM/plotSheepAcceleration.py
python ~/Documents/Uni/Documents/CSM/plotSheepVelocity.py
python ~/Documents/Uni/Documents/CSM/plotDogAcceleration.py
python ~/Documents/Uni/Documents/CSM/plotDogVelocity.py
python ~/Documents/Uni/Documents/CSM/plotTimeStep.py
python ~/Documents/Uni/Documents/CSM/convexHull.py
python ~/Documents/Uni/Documents/CSM/plotVoronoiDistances.py
