#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Apr 19 10:56:26 2021

@author: thanh
"""



from sympy import *
import IALib as IA
import MPFRGen as MPFRGen
import InterReduce as IR

numrun=1000



def makesingletestfile_nomanual(funcname,randrange):
    IA.IAprogram(funcname)
    filename=funcname+'/'+"timetestsingle.c"
    file=open(filename,'w')
    headfile=open('header.txt','r')
    codelines=headfile.readlines()
    for o in codelines:
        file.write(o)
    headfile.close()
    
    file.write("//FPSyn Function\n")    
    fpsynfile=open(funcname+'/'+funcname+ "_FPSyn.txt",'r')
    codelines=fpsynfile.readlines()
    fpsynfunc=codelines[0][5:len(codelines[0])-2]
    varlist = fpsynfunc[fpsynfunc.find('(')+1:fpsynfunc.find(')')]
    varlist = varlist.split(',')
    varlist = varlist[:len(varlist)-2]
    varlist = [v[5:] for v in varlist]

    strvl=''
    for v in varlist:
        strvl=strvl+v+','
    strvl=strvl[:len(strvl)-1]    
    
    for o in codelines:
        file.write(o)
    fpsynfile.close()
    
    file.write("//IA Function\n")    
    IAfile=open(funcname+'/'+funcname+ "_IA.txt",'r')
    codelines=IAfile.readlines()
    IAfunc=codelines[0][5:len(codelines[0])-2]
    for o in codelines:
        file.write(o)
    IAfile.close()    
    
    file.write("//Main test Function\n")    
    bodyfile=open("timesinglemain_nomanual.txt",'r')
    codelines=bodyfile.readlines()
    for o in codelines:
        if "range=" in o:
            file.write('int range = '+str(randrange)+';\n')
        elif 'int i' in o:
            file.write('REAL '+strvl+';\n')
            file.write(o)
        elif 'openfile' in o:
            file.write('file=fopen(\"'+funcname+'_singletimetestreport.txt\",\"w");\n')
        elif 'varlist=randfl(range)' in o:
            for v in varlist:
                file.write(v+"=randfl(range);\n")
        elif 'FPSynx10' in o:
            for i in range(numrun):
                file.write("FPSyn_"+funcname+"("+strvl+',&res,&err);\n')
        elif 'IAx10' in o:
            for i in range(numrun):
                file.write("IA_"+funcname+"("+strvl+',&res,&herr,&lerr);\n')
        else:                
            file.write(o)
    bodyfile.close()       
    
    file.close()


def makesingletestfile_manual(funcname,randrange):
    IA.IAprogram(funcname)
    filename=funcname+'/'+"timetestsingle.c"
    file=open(filename,'w')
    headfile=open('header.txt','r')
    codelines=headfile.readlines()
    for o in codelines:
        file.write(o)
    headfile.close()
    
    file.write("//FPSyn Function\n")    
    fpsynfile=open(funcname+'/'+funcname+ "_FPSyn.txt",'r')
    codelines=fpsynfile.readlines()
    fpsynfunc=codelines[0][5:len(codelines[0])-2]
    varlist = fpsynfunc[fpsynfunc.find('(')+1:fpsynfunc.find(')')]
    varlist = varlist.split(',')
    varlist = varlist[:len(varlist)-2]
    varlist = [v[5:] for v in varlist]

    strvl=''
    for v in varlist:
        strvl=strvl+v+','
    strvl=strvl[:len(strvl)-1]    
    
    for o in codelines:
        file.write(o)
    fpsynfile.close()

    file.write("//Manual Function\n")    
    Manualfile=open(funcname+'/'+funcname+ "_Manual.txt",'r')
    codelines=Manualfile.readlines()
    for o in codelines:
        file.write(o)
    Manualfile.close()    
    
    file.write("//IA Function\n")    
    IAfile=open(funcname+'/'+funcname+ "_IA.txt",'r')
    codelines=IAfile.readlines()
    for o in codelines:
        file.write(o)
    IAfile.close()    
    
    file.write("//Main test Function\n")    
    bodyfile=open("timesinglemain_manual.txt",'r')
    codelines=bodyfile.readlines()
    for o in codelines:
        if "range=" in o:
            file.write('int range = '+str(randrange)+';\n')
        elif 'int i' in o:
            file.write('REAL '+strvl+';\n')
            file.write(o)
        elif 'openfile' in o:
            file.write('file=fopen(\"'+funcname+'_singletimetestreport.txt\",\"w");\n')
        elif 'varlist=randfl(range)' in o:
            for v in varlist:
                file.write(v+"=randfl(range);\n")
        elif 'FPSynx10' in o:
            for i in range(numrun):
                file.write("FPSyn_"+funcname+"("+strvl+',&res,&err);\n')
        elif 'Manualx10' in o:
            for i in range(numrun):
                file.write("Manual_"+funcname+"("+strvl+',&res,&err);\n')                
        elif 'IAx10' in o:
            for i in range(numrun):
                file.write("IA_"+funcname+"("+strvl+',&res,&herr,&lerr);\n')
        else:                
            file.write(o)
    bodyfile.close()       
    
    file.close()

def makeprobtestfile_nomanual(funcname,randrange):
    IA.IAprogram(funcname)
    filename=funcname+'/'+"probtest.c"
    file=open(filename,'w')
    headfile=open('header.txt','r')
    codelines=headfile.readlines()
    for o in codelines:
        file.write(o)
    headfile.close()
    
    file.write("//FPSyn Function\n")    
    fpsynfile=open(funcname+'/'+funcname+ "_FPSyn.txt",'r')
    codelines=fpsynfile.readlines()
    fpsynfunc=codelines[0][5:len(codelines[0])-2]
    varlist = fpsynfunc[fpsynfunc.find('(')+1:fpsynfunc.find(')')]
    varlist = varlist.split(',')
    varlist = varlist[:len(varlist)-2]
    varlist = [v[5:] for v in varlist]

    strvl=''
    for v in varlist:
        strvl=strvl+v+','
    strvl=strvl[:len(strvl)-1]    
    
    for o in codelines:
        file.write(o)
    fpsynfile.close()



    file.write("//IA Function\n")    
    IAfile=open(funcname+'/'+funcname+ "_IA.txt",'r')
    codelines=IAfile.readlines()
    for o in codelines:
        file.write(o)
    IAfile.close()    
    
    file.write("//Main test Function\n")    
    bodyfile=open("probmain_nomanual.txt",'r')
    codelines=bodyfile.readlines()
    for o in codelines:
        if "range=" in o:
            file.write('int range = '+str(randrange)+';\n')
        elif 'int i' in o:
            file.write('REAL '+strvl+';\n')
            file.write(o)
        elif 'openreport' in o:
            file.write('freport=fopen(\"'+'ProbReport.txt\",\"w");\n')
        elif 'opensave' in o:
            file.write('fsave=fopen(\"'+'ProbSave.txt\",\"w");\n')
        elif 'varlist=randfl(range)' in o:
            for v in varlist:
                file.write(v+"=randfl(range);\n")
        elif 'FPSynFunc;' in o:
            file.write("FPSyn_"+funcname+"("+strvl+',&res,&err);\n')
        elif 'IAFunc;' in o:
            file.write("IA_"+funcname+"("+strvl+',&res,&herr,&lerr);\n')            
        else:                
            file.write(o)
    bodyfile.close()       
    
    file.close()


def makeprobtestfile_manual(funcname,randrange):
    IA.IAprogram(funcname)
    filename=funcname+'/'+"probtest.c"
    file=open(filename,'w')
    headfile=open('header.txt','r')
    codelines=headfile.readlines()
    for o in codelines:
        file.write(o)
    headfile.close()
    
    file.write("//FPSyn Function\n")    
    fpsynfile=open(funcname+'/'+funcname+ "_FPSyn.txt",'r')
    codelines=fpsynfile.readlines()
    fpsynfunc=codelines[0][5:len(codelines[0])-2]
    varlist = fpsynfunc[fpsynfunc.find('(')+1:fpsynfunc.find(')')]
    varlist = varlist.split(',')
    varlist = varlist[:len(varlist)-2]
    varlist = [v[5:] for v in varlist]

    strvl=''
    for v in varlist:
        strvl=strvl+v+','
    strvl=strvl[:len(strvl)-1]    
    
    for o in codelines:
        file.write(o)
    fpsynfile.close()

    file.write("//Manual Function\n")    
    Manualfile=open(funcname+'/'+funcname+ "_Manual.txt",'r')
    codelines=Manualfile.readlines()
    for o in codelines:
        file.write(o)
    Manualfile.close()    
    
    file.write("//IA Function\n")    
    IAfile=open(funcname+'/'+funcname+ "_IA.txt",'r')
    codelines=IAfile.readlines()
    for o in codelines:
        file.write(o)
    IAfile.close()    
    
    file.write("//Main test Function\n")    
    bodyfile=open("probmain_manual.txt",'r')
    codelines=bodyfile.readlines()
    for o in codelines:
        if "range=" in o:
            file.write('int range = '+str(randrange)+';\n')
        elif 'int i' in o:
            file.write('REAL '+strvl+';\n')
            file.write(o)
        elif 'openreport' in o:
            file.write('freport=fopen(\"'+'ProbReport.txt\",\"w");\n')
        elif 'opensave' in o:
            file.write('fsave=fopen(\"'+'ProbSave.txt\",\"w");\n')
        elif 'varlist=randfl(range)' in o:
            for v in varlist:
                file.write(v+"=randfl(range);\n")
        elif 'FPSynFunc;' in o:
            file.write("FPSyn_"+funcname+"("+strvl+',&res,&err);\n')
        elif 'IAFunc;' in o:
            file.write("IA_"+funcname+"("+strvl+',&res,&herr,&lerr);\n')
        elif 'ManualFunc;' in o:
            file.write("Manual_"+funcname+"("+strvl+',&res,&merr);\n')                      
        else:                
            file.write(o)
    bodyfile.close()       
    
    file.close()
    

def makeruntestfile_nomanual(funcname,randrange):
    MPFRGen.Exactprogram(funcname)
    filename=funcname+'/'+"runtest.c"
    file=open(filename,'w')
    headfile=open('header.txt','r')
    codelines=headfile.readlines()
    for o in codelines:
        file.write(o)
    headfile.close()
    
    file.write("//FPSyn Function\n")    
    fpsynfile=open(funcname+'/'+funcname+ "_FPSyn.txt",'r')
    codelines=fpsynfile.readlines()
    fpsynfunc=codelines[0][5:len(codelines[0])-2]
    varlist = fpsynfunc[fpsynfunc.find('(')+1:fpsynfunc.find(')')]
    varlist = varlist.split(',')
    varlist = varlist[:len(varlist)-2]
    varlist = [v[5:] for v in varlist]

    strvl=''
    for v in varlist:
        strvl=strvl+v+','
    strvl=strvl[:len(strvl)-1]    
    
    for o in codelines:
        file.write(o)
    fpsynfile.close()

    file.write("//IA Function\n")    
    IAfile=open(funcname+'/'+funcname+ "_IA.txt",'r')
    codelines=IAfile.readlines()
    for o in codelines:
        file.write(o)
    IAfile.close()    
    
    file.write("//Exact Function\n")    
    Exactfile=open(funcname+'/'+funcname+ "_exact.txt",'r')
    codelines=Exactfile.readlines()
    for o in codelines:
        file.write(o)
    Exactfile.close()   
    
    file.write("//FPSyn Adapt Function\n")
    proto = codelines[0]
    file.write(proto.replace("Exact","FPSyn_Adapt"))
    file.write("REAL res, err;\n")
    file.write("FPSyn_"+funcname+"("+strvl+',&res,&err);\n')
    file.write("if ((res<err)&&(res>-err)) "+ "Exact_"+funcname+"("+strvl+');\n')
    file.write("}\n")
    
    file.write("//IA Adapt Function\n")
    proto = codelines[0]
    file.write(proto.replace("Exact","IA_Adapt"))
    file.write("REAL res, lerr, herr;\n")
    file.write("IA_"+funcname+"("+strvl+',&res,&herr,&lerr);\n')
    file.write("if ((0<herr)&&(0>lerr)) "+ "Exact_"+funcname+"("+strvl+');\n')
    file.write("}\n")    
    
    
    file.write("//Main test Function\n")    
    bodyfile=open("runmain_nomanual.txt",'r')
    codelines=bodyfile.readlines()
    for o in codelines:
        if "range=" in o:
            file.write('int range = '+str(randrange)+';\n')
        elif 'int i' in o:
            file.write('REAL '+strvl+';\n')
            file.write(o)
        elif 'openfile' in o:
            file.write('file=fopen(\"'+funcname+'_singletimetestreport.txt\",\"w");\n')
        elif 'varlist=randfl(range)' in o:
            for v in varlist:
                file.write(v+"=randfl(range);\n")
        elif 'FPSynx10' in o:
            for i in range(numrun):
                file.write("FPSyn_Adapt_"+funcname+"("+strvl+');\n')
        elif 'IAx10' in o:
            for i in range(numrun):
                file.write("IA_Adapt_"+funcname+"("+strvl+');\n')
        else:                
            file.write(o)
    bodyfile.close()       
    
    file.close()



def makeruntestfile_manual(funcname,randrange):
    MPFRGen.Exactprogram(funcname)
    filename=funcname+'/'+"runtest.c"
    file=open(filename,'w')
    headfile=open('header.txt','r')
    codelines=headfile.readlines()
    for o in codelines:
        file.write(o)
    headfile.close()
    
    file.write("//FPSyn Function\n")    
    fpsynfile=open(funcname+'/'+funcname+ "_FPSyn.txt",'r')
    codelines=fpsynfile.readlines()
    fpsynfunc=codelines[0][5:len(codelines[0])-2]
    varlist = fpsynfunc[fpsynfunc.find('(')+1:fpsynfunc.find(')')]
    varlist = varlist.split(',')
    varlist = varlist[:len(varlist)-2]
    varlist = [v[5:] for v in varlist]

    strvl=''
    for v in varlist:
        strvl=strvl+v+','
    strvl=strvl[:len(strvl)-1]    
    
    for o in codelines:
        file.write(o)
    fpsynfile.close()

    file.write("//Manual Function\n")    
    Manualfile=open(funcname+'/'+funcname+ "_Manual.txt",'r')
    codelines=Manualfile.readlines()
    for o in codelines:
        file.write(o)
    Manualfile.close()    
    
    file.write("//IA Function\n")    
    IAfile=open(funcname+'/'+funcname+ "_IA.txt",'r')
    codelines=IAfile.readlines()
    for o in codelines:
        file.write(o)
    IAfile.close()    

    file.write("//Exact Function\n")    
    Exactfile=open(funcname+'/'+funcname+ "_exact.txt",'r')
    codelines=Exactfile.readlines()
    for o in codelines:
        file.write(o)
    Exactfile.close()   

    file.write("//FPSyn Adapt Function\n")
    proto = codelines[0]
    file.write(proto.replace("Exact","FPSyn_Adapt"))
    file.write("REAL res, err;\n")
    file.write("FPSyn_"+funcname+"("+strvl+',&res,&err);\n')
    file.write("if ((res<err)&&(res>-err)) "+ "Exact_"+funcname+"("+strvl+');\n')
    file.write("}\n")

    file.write("//Manual Adapt Function\n")
    proto = codelines[0]
    file.write(proto.replace("Exact","Manual_Adapt"))
    file.write("REAL res, err;\n")
    file.write("Manual_"+funcname+"("+strvl+',&res,&err);\n')
    file.write("if ((res<err)&&(res>-err)) "+ "Exact_"+funcname+"("+strvl+');\n')
    file.write("}\n")

    
    file.write("//IA Adapt Function\n")
    proto = codelines[0]
    file.write(proto.replace("Exact","IA_Adapt"))
    file.write("REAL res, lerr, herr;\n")
    file.write("IA_"+funcname+"("+strvl+',&res,&herr,&lerr);\n')
    file.write("if ((0<herr)&&(0>lerr)) "+ "Exact_"+funcname+"("+strvl+');\n')
    file.write("}\n")    
    
    file.write("//Main test Function\n")    
    bodyfile=open("runmain_manual.txt",'r')
    codelines=bodyfile.readlines()
    for o in codelines:
        if "range=" in o:
            file.write('int range = '+str(randrange)+';\n')
        elif 'int i' in o:
            file.write('REAL '+strvl+';\n')
            file.write(o)
        elif 'openfile' in o:
            file.write('file=fopen(\"'+funcname+'_singletimetestreport.txt\",\"w");\n')
        elif 'varlist=randfl(range)' in o:
            for v in varlist:
                file.write(v+"=randfl(range);\n")
        elif 'FPSynx10' in o:
            for i in range(numrun):
                file.write("FPSyn_Adapt_"+funcname+"("+strvl+');\n')
        elif 'Manualx10' in o:
            for i in range(numrun):
                file.write("Manual_Adapt_"+funcname+"("+strvl+');\n')                
        elif 'IAx10' in o:
            for i in range(numrun):
                file.write("IA_Adapt_"+funcname+"("+strvl+');\n')
        else:                
            file.write(o)
    bodyfile.close()       
    
    file.close()
    

for e in ['Orient2D','Orient3D','InCircle','InSphere','Polynomial','ConvexHullArea']:
    IR.reduceinter(e)
rangedic = {'Orient2D':126,'Orient3D':126,'InCircle':126,'InSphere':126,'Polynomial':126,'ConvexHullArea':126,'Intersection2D':126,'Intersection3D':30}
manuallist=['Orient2D','Orient3D','InCircle','InSphere']
nomanuallist=['Polynomial','ConvexHullArea','Intersection2D','Intersection3D']
for file in manuallist:
    makesingletestfile_manual(file,rangedic[file])
    makeprobtestfile_manual(file,rangedic[file])
    makeruntestfile_manual(file,rangedic[file])
for file in nomanuallist:
    makesingletestfile_nomanual(file,rangedic[file])
    makeprobtestfile_nomanual(file,rangedic[file])    
    makeruntestfile_nomanual(file,rangedic[file])