# -*- coding: utf-8 -*-
"""
Created on Mon Sep 14 11:11:40 2020

@author: thanh
"""

import FPSynLib as FPS
import copy
import math
from sympy import *
import itertools
import random

x, y, z, t, m,n,eps = symbols('x y z t m n eps')
a,b,c,d,e,f,g,h =symbols('a, b c d e f g h')

t=x
p=1+x
for i in range(9):
    j = i+2
    t= t*x*(1/j)
    p=p+t

##Running---------------------------###



p=FPS.FPSynthesis(p,"Polynomial")

