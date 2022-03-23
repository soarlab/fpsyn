#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Apr 19 10:56:26 2021

@author: thanh
"""


import copy
import math
from sympy import *
import itertools
import random
import IALib as IA

def makesingletestfile_nomanual(funcname):
    filename=funcname+'/'+funcname+ "_timetestsingle.txt"
    file=open(filename,'w')
    headfile=open('header.txt','r')
    codelines=headfile.readlines()
    for o in codelines:
        file.write(o)
    headfile.close()
    
    file.write("//FPSyn Function\n")    
    fpsynfile=open(funcname+'/'+funcname+ "_program.txt",'r')
    codelines=fpsynfile.readlines()
    fpsynfunc=codelines[0]
    varlist = fpsynfunc[fpsynfunc.find('(')+1:fpsynfunc.find(')')]
    varlist = varlist.split(',')
    varlist = varlist[:len(varlist)-2]
    print(varlist)
    for o in codelines:
        file.write(o)
    fpsynfile.close()
    
    file.write("//IA Function\n")    
    IAfile=open(funcname+'/'+funcname+ "_IA.txt",'r')
    codelines=IAfile.readlines()
    IAfunc=codelines[0]
    for o in codelines:
        file.write(o)
    IAfile.close()    
    
    file.write("//Main test Function\n")    
    bodyfile=open("timesinglemain_nomanual.txt",'r')
    codelines=bodyfile.readlines()
    for o in codelines:
        file.write(o)
    bodyfile.close()       
    
    file.close()


makesingletestfile_nomanual('Polynomial')