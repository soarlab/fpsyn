# -*- coding: utf-8 -*-
"""
Created on Mon Sep 14 11:11:40 2020

@author: thanh
"""

import FPSynDivLib as FPSDiv
import copy
import math
from sympy import *
import itertools
import random

x, y, z, t, m,n,eps = symbols('x y z t m n eps')
ax1,ay1,az1,ax2,ay2,az2,ax3,ay3,az3 =symbols('ax1 ay1 az1 ax2 ay2 az2 ax3 ay3 az3')
bx1,by1,bz1,bx2,by2,bz2,bx3,by3,bz3 =symbols('bx1 by1 bz1 bx2 by2 bz2 bx3 by3 bz3')
cx1,cy1,cz1,cx2,cy2,cz2,cx3,cy3,cz3 =symbols('cx1 cy1 cz1 cx2 cy2 cz2 cx3 cy3 cz3')


plane1x=(by1-ay1)*(cz1-az1)-(cy1-ay1)*(bz1-az1)
plane1y=(cx1-ax1)*(bz1-az1)-(bx1-ax1)*(cz1-az1)
plane1z=(bx1-ax1)*(cy1-ay1)-(cx1-ax1)*(by1-ay1)
plane1t=ax1*plane1x+ay1*plane1y+az1*plane1z

plane2x=(by2-ay2)*(cz2-az2)-(cy2-ay2)*(bz2-az2)
plane2y=(cx2-ax2)*(bz2-az2)-(bx2-ax2)*(cz2-az2)
plane2z=(bx2-ax2)*(cy2-ay2)-(cx2-ax2)*(by2-ay2)
plane2t=ax2*plane2x+ay2*plane2y+az2*plane2z

plane3x=(by3-ay3)*(cz3-az3)-(cy3-ay3)*(bz3-az3)
plane3y=(cx3-ax3)*(bz3-az3)-(bx3-ax3)*(cz3-az3)
plane3z=(bx3-ax3)*(cy3-ay3)-(cx3-ax3)*(by3-ay3)
plane3t=ax3*plane3x+ay3*plane3y+az3*plane3z

den1=plane1x*(plane2y*plane3z-plane3y*plane2z)
den2=plane1y*(plane2x*plane3z-plane3x*plane2z)
den3=plane1z*(plane2x*plane3y-plane3x*plane2y)

den=den1-den2+den3

num1=plane1x*(plane2y*plane3t-plane3y*plane2t)
num2=plane1y*(plane2x*plane3t-plane3x*plane2t)
num3=plane1t*(plane2x*plane3y-plane3x*plane2y)

num=num1-num2+num3



##Running---------------------------###
#div=Function('div')
#p=div(num,den)
#gg=FPSDiv.orgexpedit(p)

errp=FPSDiv.FPSynthesisDiv(num,den,"Intersection3D")

