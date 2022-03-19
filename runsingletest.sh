#!/bin/sh
cd Orient2D
gcc -o timetestsingle timetestsingle.c -lm
./timetestsingle
cd ..
cd Orient3D
gcc -o timetestsingle timetestsingle.c -lm
./timetestsingle
cd ..
cd InCircle
gcc -o timetestsingle timetestsingle.c -lm
./timetestsingle
cd ..
cd InSphere
gcc -o timetestsingle timetestsingle.c -lm
./timetestsingle
cd ..
cd ConvexHullArea
gcc -o timetestsingle timetestsingle.c -lm
./timetestsingle
cd ..
cd Intersection2D
gcc -o timetestsingle timetestsingle.c -lm
./timetestsingle
cd ..
cd Intersection3D
gcc -o timetestsingle timetestsingle.c -lm
./timetestsingle
cd ..
cd Polynomial
gcc -o timetestsingle timetestsingle.c -lm
./timetestsingle
cd ..