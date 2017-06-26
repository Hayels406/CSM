#!/bin/sh
python CSM.py
mkdir -p frames
rm -f frames/*
python ~/Documents/CSM/plotFrames.py
p=`pwd`
vid=`echo ${p##*/}`
ffmpeg -framerate 20 -pattern_type glob -i 'frames/*.png' -vcodec mpeg4 -qscale:v 2 frames/${vid}.mp4
mkdir -p plots
python ~/Documents/CSM/plotSheepAcceleration.py
python ~/Documents/CSM/plotSheepVelocity.py
python ~/Documents/CSM/plotDogAcceleration.py
python ~/Documents/CSM/plotDogVelocity.py
python ~/Documents/CSM/plotTimeStep.py
