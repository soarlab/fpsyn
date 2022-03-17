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

p=x*y*z+m*n*t
p1=(x*y+x*x+x*y*z)*(z*z+y)*(x+y)+m*n+t
p2=(x+y+z)*(x+t)+y*m
p2=expand(p2)


ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez =symbols('ax ay az bx by bz cx cy cz dx dy dz ex ey ez')



incira = ((ax-dx)**2 + (ay-dy)**2)*((bx-dx)*(cy-dy)-(by-dy)*(cx-dx)) 
incirb = ((bx-dx)**2 + (by-dy)**2)*((cx-dx)*(ay-dy)-(cy-dy)*(ax-dx)) 
incirc = ((cx-dx)**2 + (cy-dy)**2)*((ax-dx)*(by-dy)-(ay-dy)*(bx-dx)) 
incircle = incira + incirb + incirc


    

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


p=insphere
p=FPS.FPSynthesis(p,"InSphere")

