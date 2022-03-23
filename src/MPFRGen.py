# -*- coding: utf-8 -*-
"""
Created on Thu Nov  5 15:08:52 2020

@author: thanh
"""
from sympy import *

def breakcodeline(cl,varlist):
    left=cl[:cl.find('=')]
    op='n'
    out=("","","","")
    if '+' in cl:
        op='+'
        term1=cl[cl.find('=')+1:cl.find('+')]
        term2=cl[cl.find('+')+1:cl.find(';')]
        
    elif '*' in cl:
        op='*'
        term1=cl[cl.find('=')+1:cl.find('*')]
        term2=cl[cl.find('*')+1:cl.find(';')]
         
    elif '/' in cl:
        op='/'
        term1=cl[cl.find('=')+1:cl.find('/')]
        term2=cl[cl.find('/')+1:cl.find(';')]
        
    elif '-' in cl:
        op='-'
        term1=cl[cl.find('=')+1:cl.find('-')]
        term2=cl[cl.find('-')+1:cl.find(';')]
        if term1=='':
            term1='0'
    
    if not isnum(term1):
        term1="MP"+term1
    else:
        if '.' not in term1:
            term1 = term1 + ".0"
    if not isnum(term2):
        term2="MP"+term2
    else:
        if '.' not in term2:
            term2 = term2 + ".0"
    left="MP"+left[1:]    
    out=(left,op,term1,term2)        
    return out

def readinput(filename, varlist):
    file=open(filename,'r')
    codelines=file.readlines()
    code=[]
    out=[]
    i=0
    for c in codelines:
        c=c[4:]
        if (i==0) and ('=' in c):
            code.append(c)
            i=1
            continue
        if (i==1) and ('=' in c):
            code.append(c)
        if (i==1) and ('=' not in c):
            i=2
    for c in code:
        out.append(breakcodeline(c,varlist))
    return out

def readinputraw(filename):
    file=open(filename,'r')
    codelines=file.readlines()
    proto = codelines[0]
    proto=proto.replace("FPSyn", "Exact")
    proto=proto.replace(",REAL* result, REAL* error_bound", "")
    out=[]
    out.append(proto)
    code=[]

    i=0
    for c in codelines:
        if (i==0) and ('=' in c):
            code.append(c)
            i=1
            continue
        if (i==1) and ('=' in c):
            code.append(c)
        if (i==1) and ('=' not in c):
            i=2
    for c in code:
        out.append(c)
    return out




def isnum(s):
    if '.' in s:
        return True
    if '-' == s[0]:
        return True
    if s.isnumeric(): 
        return True
    return False

def generatecode(stlist, varlist):
    out=[]
    #n=52*(2**(len(stlist)))
    n = 1024
    for v in varlist:
        out.append("mpfr_t MP"+ v)
        out.append("mpfr_init2(MP" + v + "," + str(n)+")")
        
    
    rndm = "MPFR_RNDN"
    for st in stlist:
        
        left,op,term1,term2=st
        out.append("mpfr_t "+ left)
        out.append("mpfr_init2(" + left+ "," + str(n)+")")

        if op=='+':
            if isnum(term1):
                out.append("mpfr_add_d(" + left + "," + term2 + "," + term1 + "," + rndm + ")" )
            elif isnum(term2):
                out.append("mpfr_add_d(" + left + "," + term1 + "," + term2 + "," + rndm + ")" )
            else:
                out.append("mpfr_add(" + left + "," + term1 + "," + term2 + "," + rndm + ")" )

        if op=='-':
            if isnum(term1):
                out.append("mpfr_d_sub(" + left + "," + term1 + "," + term2 + "," + rndm + ")" )
            elif isnum(term2):
                out.append("mpfr_sub_d(" + left + "," + term1 + "," + term2 + "," + rndm + ")" )
            else:
                out.append("mpfr_sub(" + left + "," + term1 + "," + term2 + "," + rndm + ")" )

        if op=='*':
            if isnum(term1):
                out.append("mpfr_mul_d(" + left + "," + term2 + "," + term1 + "," + rndm + ")" )
            elif isnum(term2):
                out.append("mpfr_mul_d(" + left + "," + term2 + "," + term1 + "," + rndm + ")" )
            else:
                out.append("mpfr_mul(" + left + "," + term1 + "," + term2 + "," + rndm + ")" )
           
        elif op=='/':
            if isnum(term1):
                out.append("mpfr_d_div(" + left + "," + term1 + "," + term2 + "," + rndm + ")" )
            elif isnum(term2):
                out.append("mpfr_div_d(" + left + "," + term1 + "," + term2 + "," + rndm + ")" )
            else:
                out.append("mpfr_div(" + left + "," + term1 + "," + term2 + "," + rndm + ")" )
    
    for v in varlist:
        out.append("mpfr_clear(MP"+ v+ ")")
    for st in stlist:
        left,op,term1,term2=st
        out.append("mpfr_clear(" + left + ")")

    return out


def writecode(s,file):
    for e in s:
        if ('=' in e) and ('int' not in e) and ("if" not in e):
            e=e.replace(' ','')
            if (e[0]!='*') and ("if" not in e) and ("else" not in e) and ("df" not in e):
                e="REAL "+e
            if ("df" in e):
                e=e[2:]
        if(str(e)!='')  :
            file.write(str(e))
            if (e!="}") and ("{" not in e):
                file.write(';')
            file.write('\n')

def Exactprogram(funcname):        
    filename=funcname+'/'+funcname+ "_program.txt"
    file=open(funcname+'/'+funcname+"_exact.txt",'w')
    out=readinputraw(filename)
    #for o in out:
    file.write(out[0])
    
    pfile=open(funcname+'/'+funcname+ "_program.txt",'r')    
    codelines=pfile.readlines()
    fpsynfunc=codelines[0][5:len(codelines[0])-2]
    varlist = fpsynfunc[fpsynfunc.find('(')+1:fpsynfunc.find(')')]
    varlist = varlist.split(',')
    varlist = varlist[:len(varlist)-2]
    varlist = [v[5:] for v in varlist]
    pfile.close()
    
    out=generatecode(readinput(filename,varlist), varlist)
    writecode(out,file)
    
    file.write('}\n')
    file.close()


#IAprogram('Intersection2D')






            

        