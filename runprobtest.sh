#!/bin/sh
cd Orient2D
gcc -o probtest probtest.c -lm
./probtest
cd ..
cd Orient3D
gcc -o probtest probtest.c -lm
./probtest
cd ..
cd InCircle
gcc -o probtest probtest.c -lm
./probtest
cd ..
cd InSphere
gcc -o probtest probtest.c -lm
./probtest
cd ..
cd ConvexHullArea
gcc -o probtest probtest.c -lm
./probtest
cd ..
cd Intersection2D
gcc -o probtest probtest.c -lm
./probtest
cd ..
cd Intersection3D
gcc -o probtest probtest.c -lm
./probtest
cd ..
cd Polynomial
gcc -o probtest probtest.c -lm
./probtest
cd ..
