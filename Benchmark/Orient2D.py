# -*- coding: utf-8 -*-
"""
Created on Mon Sep 14 11:11:40 2020

@author: thanh
"""

import FPSynLib as FPS
from sympy import *




ax,ay,bx,by,cx,cy =symbols('ax ay bx by cx cy')

orient2D=(ax-cx)*(by-cy)-(bx-cx)*(ay-cy)

##Running---------------------------###


exp=orient2D
errorexp=FPS.FPSynthesis(exp,"Orient2D")

