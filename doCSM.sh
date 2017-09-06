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

MAC=/Users/hayleymoore/Documents/Uni/Documents/CSM
if [ -f $MAC/params.py.dist ]; then
  open frames/${vid}.mp4
fi

mkdir -p plots
mkdir -p output
python ~/Documents/Uni/Documents/CSM/plotSheepAcceleration.py
python ~/Documents/Uni/Documents/CSM/plotSheepVelocity.py
python ~/Documents/Uni/Documents/CSM/plotDogAcceleration.py
python ~/Documents/Uni/Documents/CSM/plotDogVelocity.py
python ~/Documents/Uni/Documents/CSM/plotTimeStep.py
python ~/Documents/Uni/Documents/CSM/plotConvexHull.py
python ~/Documents/Uni/Documents/CSM/plotVoronoiDistances.py
python ~/Documents/Uni/Documents/CSM/findClusters.py

if [ -f $MAC/params.py.dist ]; then
  terminal-notifier -title "CSM" -message "Done with simulation" -subtitle ${vid} -sound 'glass'
  open frames/${vid}.mp4
fi
