#!/bin/sh
cd Orient2D
cc -o runtest runtest.c -lmpfr -lgmp -lm
./runtest
cd ..
cd Orient3D
cc -o runtest runtest.c -lmpfr -lgmp -lm
./runtest
cd ..
cd InCircle
cc -o runtest runtest.c -lmpfr -lgmp -lm
./runtest
cd ..
cd InSphere
cc -o runtest runtest.c -lmpfr -lgmp -lm
./runtest
cd ..
cd ConvexHullArea
cc -o runtest runtest.c -lmpfr -lgmp -lm
./runtest
cd ..
cd Intersection2D
cc -o runtest runtest.c -lmpfr -lgmp -lm
./runtest
cd ..
cd Intersection3D
cc -o runtest runtest.c -lmpfr -lgmp -lm
./runtest
cd ..
cd Polynomial
cc -o runtest runtest.c -lmpfr -lgmp -lm
./runtest
cd ..