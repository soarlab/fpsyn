# -*- coding: utf-8 -*-
"""
Created on Thu Nov  5 15:08:52 2020

@author: thanh
"""
from sympy import *

def breakcodeline(cl):
    left=cl[:cl.find('=')]
    op='n'
    out=("","","","")
    if ((('A' in cl) or ('B' in cl)) or ('E' in cl)) or ('?' in cl):
        op="none"
        term1=cl[cl.find('=')+1:cl.find(';')]
        for i in ['+','-','*','/','?']:
            term1 = term1.replace(i, ' '+i+' ')
        term2=""
        out=(left,op,term1,term2)  
    elif '+' in cl:
        op='+'
        term1=cl[cl.find('=')+1:cl.find('+')]
        term2=cl[cl.find('+')+1:cl.find(';')]
        out=(left,op,term1,term2)
    elif '*' in cl[1:]:
        op='*'
        if cl[0]!="*":
            term1=cl[cl.find('=')+1:cl.find('*')]
            term2=cl[cl.find('*')+1:cl.find(';')]
            out=(left,op,term1,term2)
        else:
            cl=cl[1:]
            term1=cl[cl.find('=')+1:cl.find('*')]
            term2=cl[cl.find('*')+1:cl.find(';')]
            out=(left,op,term1,term2)            
    elif '/' in cl:
        op='/'
        term1=cl[cl.find('=')+1:cl.find('/')]
        term2=cl[cl.find('/')+1:cl.find(';')]
        out=(left,op,term1,term2)
    elif '-' in cl:
        op='-'
        term1=cl[cl.find('=')+1:cl.find('-')]
        term2=cl[cl.find('-')+1:cl.find(';')]
        if term1=='':
            term1='0'
        out=(left,op,term1,term2)
    elif 'fabs' in cl:
        op='fabs'
        term1=cl[cl.find('(')+1:cl.find(')')]
        term2=""
        out=(left+"f",op,term1,term2)
    else:
        op="none"
        term1=cl[cl.find('=')+1:cl.find(';')]
        term2=""
        out=(left,op,term1,term2)  

    return out

def readinput(filename):
    file=open(filename,'r')
    codelines=file.readlines()
    file.close()
    code=[]
    out=[]
    for c in codelines:
        if "REAL" in c:
            c=c[5:]
        if '=' in c:
            code.append(c)
    for c in code:
        out.append(breakcodeline(c))
    
    out2=[]
    ab = []
    for (left,op,term1,term2) in out:
        if op=="fabs":
            ab.append(term1)
    
    for (left,op,term1,term2) in out:
        if "tem" in left:
            if (("tem" not in term1) and ("t" in term1)) and (term1 in ab):
                term1=term1+"f"
            if (("tem" not in term2) and ("t" in term2)) and (term2 in ab):
                term2=term2+"f"
        if left in ['A','B','E','*error_bound']:
            for t in ab:
                term1 = term1.replace(t+" ",t+"f ")
        out2.append((left,op,term1,term2))        
        
    return out2

def readinputraw(filename):
    file=open(filename,'r')
    codelines=file.readlines()
    file.close()
    proto = codelines[0]
    proto=proto.replace("FPSyn", "IA")
    proto=proto.replace("error_bound", "hbound, REAL* lbound")
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



def leftdic(out):
    dic = {}
    for (left,op,t1,t2) in out:
        if op=="fabs":
            dic[left] = " fabs( "+t1+" ) "
        elif op =="none":
            dic[left] =" "+ t1+" "
        else:          
            dic[left] = " " + t1 + " " + op + " " + t2 + " "
    return dic

def reducedicstep(dic):   
    unique = [] 
    for left in dic:
        c = 0
        for st in dic:
            c = c + dic[st].count(left + " ")
        if (c == 1) and (left not in ['A','B','E']):
            unique.append(left)      
    for left in unique:
        for st in dic:
            dic[st] = dic[st].replace(" "+left+" ", " ( "+ dic[left] + " ) ")
        del dic[left]
    return dic
    

def reducedic(dic):
    lb =len(dic)
    la = len(reducedicstep(dic))
    if la==lb:
        return dic
    return reducedic(reducedicstep(dic))
    
def removespace(dic):
    for l in dic:
        dic[l]=dic[l].replace(" ","")
    return dic
    
def reduceinter(funcname):
    filename = funcname + "/"+funcname+"_program.txt"
    file=open(filename,'r')
    proto =file.readlines()[0]
    file.close()
    out = readinput(filename)
    dic = leftdic(out)
    rdic =   removespace(reducedic(dic))  
    filename = funcname + "/"+funcname+"_FPSyn.txt"
    file = open(filename,'w')
    file.write(proto)
    for e in rdic:
        if e[0]=="*":
            st = e + " = " + rdic[e] + ";\n"
        else:
            st = "REAL " + e + " = " + rdic[e] + ";\n"
        st=st.replace("fabs","Absolute")
        file.write(st)
    file.write("}")
    file.close()


funcname = "Intersection2D"
reduceinter(funcname)

#filename = "InCircle/InCircle_program.txt"
#filename = "Orient3D/Orient3D_program.txt"
#
#
#for e in out:
#    print(e)
#
#print("----------------")    
#dic = leftdic(out)
#
#for e in dic:
#    print(e + " = " + dic[e])
#
#print("------------------")



