# -*- coding: utf-8 -*-
"""
Created on Mon Sep 14 11:11:40 2020

@author: thanh
"""

import FPSynLib as FPS
from sympy import *



ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez =symbols('ax ay az bx by bz cx cy cz dx dy dz ex ey ez')
    

aex = ax - ex
bex = bx - ex
cex = cx - ex
dex = dx - ex

aey = ay - ey
bey = by - ey
cey = cy - ey
dey = dy - ey

aez = az - ez
bez = bz - ez
cez = cz - ez
dez = dz - ez


aexbey = aex * bey
bexaey = bex * aey
ab = aexbey - bexaey
bexcey = bex * cey
cexbey = cex * bey
bc = bexcey - cexbey
cexdey = cex * dey
dexcey = dex * cey
cd = cexdey - dexcey
dexaey = dex * aey
aexdey = aex * dey
da = dexaey - aexdey

aexcey = aex * cey
cexaey = cex * aey
ac = aexcey - cexaey
bexdey = bex * dey
dexbey = dex * bey
bd = bexdey - dexbey

abc = aez * bc - bez * ac + cez * ab
bcd = bez * cd - cez * bd + dez * bc
cda = cez * da + dez * ac + aez * cd
dab = dez * ab + aez * bd + bez * da

alift = aex * aex + aey * aey + aez * aez
blift = bex * bex + bey * bey + bez * bez
clift = cex * cex + cey * cey + cez * cez
dlift = dex * dex + dey * dey + dez * dez

insphere = (dlift * abc - clift * dab) + (blift * cda - alift * bcd)




##Running---------------------------###


exp=insphere
errorexp=FPS.FPSynthesis(p,"InSphere")

