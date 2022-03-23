# -*- coding: utf-8 -*-
"""
Created on Mon Sep 14 11:11:40 2020

@author: thanh
"""

import FPSynLib as FPS
from sympy import *


x = symbols('x')

t=x
p=1+x
for i in range(9):
    j = i+2
    t= t*x*(1/j)
    p=p+t

##Running---------------------------###



errorexp=FPS.FPSynthesis(p,"Polynomial")

