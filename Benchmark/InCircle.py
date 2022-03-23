# -*- coding: utf-8 -*-
"""
Created on Mon Sep 14 11:11:40 2020

@author: thanh
"""

import FPSynLib as FPS
from sympy import *


ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez =symbols('ax ay az bx by bz cx cy cz dx dy dz ex ey ez')



incira = ((ax-dx)**2 + (ay-dy)**2)*((bx-dx)*(cy-dy)-(by-dy)*(cx-dx)) 
incirb = ((bx-dx)**2 + (by-dy)**2)*((cx-dx)*(ay-dy)-(cy-dy)*(ax-dx)) 
incirc = ((cx-dx)**2 + (cy-dy)**2)*((ax-dx)*(by-dy)-(ay-dy)*(bx-dx)) 
incircle = incira + incirb + incirc


##Running---------------------------###

funcname = "InCircle"
exp=incircle
errexp=FPS.FPSynthesis(exp,funcname)

