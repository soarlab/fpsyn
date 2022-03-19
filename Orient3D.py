# -*- coding: utf-8 -*-
"""
Created on Mon Sep 14 11:11:40 2020

@author: thanh
"""

import FPSynLib as FPS
from sympy import *


ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz =symbols('ax ay az bx by bz cx cy cz dx dy dz')


orient3Da=(ax-dx)*((by-dy)*(cz-dz)-(bz-dz)*(cy-dy))
orient3Db=(bx-dx)*((cy-dy)*(az-dz)-(cz-dz)*(ay-dy))
orient3Dc=(cx-dx)*((ay-dy)*(bz-dz)-(az-dz)*(by-dy))
orient3D=orient3Da+orient3Db+orient3Dc


##Running---------------------------###


p=orient3D
p=FPS.FPSynthesis(p,"Orient3D")

