# -*- coding: utf-8 -*-
"""
Created on Mon Sep 14 11:11:40 2020

@author: thanh
"""

import FPSynLib as FPS
from sympy import *


numver=10
p=0

ver=[]
vec=[]

for i in range(numver):
    xt = symbols('x'+str(i))
    yt = symbols('y'+str(i))
    ver.append((xt,yt))

(x0,y0)=ver[numver-1]

for i in range(numver-1):
    (a,b)=ver[i]
    vec.append((a-x0,b-y0))

for i in range(numver-2):
    (a,b)=vec[i]
    (c,d)=vec[i+1]
    p=p+a*d-c*b    

exp=0.5*p

errexp=FPS.FPSynthesis(exp,"ConvexHullArea")

