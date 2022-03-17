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

orient2D=(ax-cx)*(by-cy)-(bx-cx)*(ay-cy)
orient3Da=(ax-dx)*((by-dy)*(cz-dz)-(bz-dz)*(cy-dy))
orient3Db=(bx-dx)*((cy-dy)*(az-dz)-(cz-dz)*(ay-dy))
orient3Dc=(cx-dx)*((ay-dy)*(bz-dz)-(az-dz)*(by-dy))
orient3D=orient3Da+orient3Db+orient3Dc

incira = ((ax-dx)**2 + (ay-dy)**2)*((bx-dx)*(cy-dy)-(by-dy)*(cx-dx)) 
incirb = ((bx-dx)**2 + (by-dy)**2)*((cx-dx)*(ay-dy)-(cy-dy)*(ax-dx)) 
incirc = ((cx-dx)**2 + (cy-dy)**2)*((ax-dx)*(by-dy)-(ay-dy)*(bx-dx)) 
incircle = incira + incirb + incirc

insphererow1=[ax-ex,ay-ey,az-ez,(ax-ex)**2+(ay-ey)**2+(az-ez)**2]
insphererow2=[bx-ex,by-ey,bz-ez,(bx-ex)**2+(by-ey)**2+(bz-ez)**2]
insphererow3=[cx-ex,cy-ey,cz-ez,(cx-ex)**2+(cy-ey)**2+(cz-ez)**2]
insphererow4=[dx-ex,dy-ey,dz-ez,(dx-ex)**2+(dy-ey)**2+(dz-ez)**2]
Mat = [insphererow1,insphererow2,insphererow3,insphererow4]
tmat=[[]]
for i in range(4):
    trow=[]
    for j in range(4):
        m=symbols('t'+str(i)+str(j))
        trow.append(m)
    tmat.append(trow)

inspherematrix=Matrix([insphererow1,insphererow2,insphererow3,insphererow4])
tmatrix=Matrix(tmat)
insphere=tmatrix.det()
for i in range(4):
    for j in range(4):
        m=symbols('t'+str(i)+str(j))
        insphere=insphere.subs(m,Mat[i][j])
    



exp=(a**6+b*b)+c*d+(a**6+b*b)*c*d


##Running---------------------------###


p=exp
p=FPS.FPSynthesis(p,"Example")

