#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Oct 23 10:02:30 2021

@author: thanh
"""
import numpy as np
import csv
def getinfo(funcname):
    filelist = ["SingleReport.txt", "ProbReport.txt", "RuntestReport.txt"]
    num=[]
    for f in filelist:
        filename=funcname+'/'+f
        file=open(filename,'r')
        codelines=file.readlines()
        for l in codelines:
            spt = l[:len(l)-1].split(' ')
            last = spt[len(spt)-1]
            if len(last)>0:
                num.append(last)
            else:
                num.append(spt[len(spt)-2])
            #num.append(float(spt[len(spt)-1]))
        file.close()
    return num

def getstd_nomanual(funcname,savename):
    filelist = [savename]
    IAList=np.zeros(10000)
    FPSList=np.zeros(10000)
    out = []
    for f in filelist:
        filename=funcname+'/'+f
        file=open(filename,'r')
        codelines=file.readlines()
        for t in range(10000):
            l=codelines[2*t]
            FPSList[t] = float(l[:len(l)-1])
            l=codelines[2*t+1]
            IAList[t] = float(l[:len(l)-1])
        file.close()    
    return [np.std(FPSList)/100, np.std(IAList)/100]


def getstd_nomanual2(funcname,savename):
    filelist = [savename]
    IAList=np.zeros(10000)
    FPSList=np.zeros(10000)
    for f in filelist:
        filename=funcname+'/'+f
        file=open(filename,'r')
        codelines=file.readlines()
        for t in range(10000):
            l=codelines[2*t]
            FPSList[t] = float(l[:len(l)-1])
            l=codelines[2*t+1]
            IAList[t] = float(l[:len(l)-1])
        file.close()    
    return IAList

def getstd_manual(funcname,savename):
    filelist = [savename]
    IAList=np.zeros(10000)
    FPSList=np.zeros(10000)
    ManList=np.zeros(10000)
    for f in filelist:
        filename=funcname+'/'+f
        file=open(filename,'r')
        codelines=file.readlines()
        for t in range(10000):
            l=codelines[3*t]
            FPSList[t] = float(l[:len(l)-1])
            l=codelines[3*t+1]
            ManList[t] = float(l[:len(l)-1])
            l=codelines[3*t+2]
            IAList[t] = float(l[:len(l)-1])
        file.close()    
    return [np.std(FPSList)/100, np.std(ManList)/100, np.std(IAList)/100]

manuallist=['Orient2D','Orient3D','InCircle','InSphere']
nomanuallist=['Polynomial','ConvexHullArea','Intersection2D','Intersection3D']
alllist = ['Orient2D','Orient3D','InCircle','InSphere','Polynomial','ConvexHullArea','Intersection2D','Intersection3D']

meanstdsingle=[]
meanstd=[]
prob =[]


ff = getstd_nomanual2('Orient2D',"SingleSave.txt")

for f in manuallist:
    info = getinfo(f)
    std1 = getstd_manual(f,"SingleSave.txt")
    std2 = getstd_manual(f,"RuntestSave.txt")
    l = [float(info[0]),std1[0],float(info[1]),std1[1],float(info[2]),std1[2]]
    meanstdsingle.append(l)
    prob.append(info[3:12])
    l = [float(info[12]),std2[0],float(info[13]),std2[1],float(info[14]),std2[2]]
    meanstd.append(l)
    
for f in nomanuallist:
    info = getinfo(f)
    std1 = getstd_nomanual(f,"SingleSave.txt")
    std2 = getstd_nomanual(f,"RuntestSave.txt")
    l = [float(info[0]),std1[0],0,0,float(info[1]),std1[1]]
    meanstdsingle.append(l)
    prob.append([info[2]]+[0]+info[3:7]+[0,0,0])
    l = [float(info[7]),std2[0],0,0,float(info[8]),std2[1]]
    meanstd.append(l)
    
meanstdsingle = np.array(meanstdsingle)
meanstd= np.array(meanstd)
prob = np.array(prob)


singleruntime = np.delete(meanstdsingle,1,1)
singleruntime = np.delete(singleruntime,2,1)
singleruntime = np.delete(singleruntime,3,1)
singleruntime = 1000000*singleruntime


fullruntime = np.delete(meanstd,1,1)
fullruntime = np.delete(fullruntime,2,1)
fullruntime = np.delete(fullruntime,3,1)
fullruntime = 1000000*fullruntime





    