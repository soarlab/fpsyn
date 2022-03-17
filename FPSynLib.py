# -*- coding: utf-8 -*-
"""
Created on Mon Aug 31 11:08:44 2020

@author: thanh
"""

# -*- coding: utf-8 -*-
"""
Created on Sun Aug  2 21:17:10 2020

@author: thanh
"""



import copy
import math
# -*- coding: utf-8 -*-


from sympy import *
import itertools
import random
count=0

def iappend(s,x):
    global count
    (n,p)=x
    #nonex ={Symbol, Integer, Float, One}
    if (n.func==Symbol) and p==1:
        return n
    if n.is_number:
        return n
    for v in s:
        (c,a,b)=v
        if (a==n) and (b==p):
            return c
    if p>1:
        i=1
        while i<p:
            iappend(s,(n,i))
            i=i*2
    
    m=symbols('tmp'+str(count))        
    s.append((m,n,p))
    count=count+1
    return m

def replacepower(p,s):
    global count
    if p.func==Pow:
        x=p.args[0]
        pw=p.args[1]
        m=symbols('tmp'+str(count))
        #nonex ={Symbol, Integer, Float, One}
        if (x.is_number) or (x.func==Symbol):
            t=iappend(s,(x,pw))
            return t
        else:
            t=replacepower(x,s)
            t1=iappend(s,(t,1))
            t2=iappend(s,(t1,pw))
            return t2

        return m 
    else:
        p1=p
        for i in p.args:
            v=replacepower(i,s)
            p1=p1.subs(i,v)
        k=iappend(s,(p1,1))
        return k


def maxpower(p,i):
    var=p.free_symbols
    mp=1
    s=[]
    for j1 in range(100):
        j=j1+2
        m=symbols('p'+str(i))
        p1=p.subs(i**j,m)
        if p!=p1:
            print(str(i)+" to " + str(j))
            mp=j
    t=2
    m=symbols(str(i))
    while t<=mp:
        s.append(str(i)+str(t)+"="+str(m)+"*"+str(m))
        m=symbols(str(i)+str(t))
        t=t*2
        
    print(s)
            
def binarylist(a):
    s=[]
    while a>1:
        i=1
        j=0
        while i<=a:
            i=i*2
            j=j+1
        a=a-(i//2)
        s.append(2**(j-1))
    if a==1:
        s.append(1)
    return s

def edit(s):
    nopow = []
    for mem in s:
        (a,b,c)=mem
        if c==1:
            nopow.append((a,b))
        else:
            bl= binarylist(c)
            if len(bl)==1:
                t=iappend(s,(b,bl[0]//2))
                nopow.append((a,t*t))
            else:
                #print(mem)
                t=iappend(s,(b,bl[0]))
                j=0
                for i in bl:
                    if j>0:
                        if i==1:
                            t=t*b
                        else:
                            t=t*iappend(s,(b,i))
                    j=j+1
                nopow.append((a,t))
               # print((a,t))
    return nopow
                    
                


def t_cal_step1(p,s):
    t=cse(p, optimizations='basic')
    numvar=len(p.free_symbols)
    p1=t[1][0]
    exp=[]
    for i in t[0]:
        (a,b)=i
        exp.append((a,b))
    #print(t)
    m=symbols('err')
    exp.append((m,p1))
    
    sublist=[]
    for mem in exp:
        (a,b)=mem
        for mem2 in sublist:
            b=b.subs(mem2[0],mem2[1])
        tmp=replacepower(b,s)
        sublist.append((a,tmp))
    
    s = edit(s)
    #numop=0
    for mem in s:
        (a,b)=mem
        #numop=numop+b.count_ops()
    
    return s

def t_cal_step2(exp):
    global count
    s=[]
    count=0
    return t_cal_step1(exp,s)





def writeassignment(asglist):
    s=[]
    for a in asglist:
        (l,r)=a
        if (r.func==Pow):
            t=r.args[0]
            s.append(str(l)+"="+str(t)+"*"+str(t))
        else:
            s.append(str(l)+"=" + str(r))
    return s


#----------------------------------------------##
    
def listbinary(n):
    if n==0:
        return [[]];
    s=listbinary(n-1)
    s0=[]
    for t in s:
        u=[0]+t
        s0.append(u) 
    for t in s:
        u=[1]+t
        s0.append(u)   
    return s0

def listpartition(exp):
    exp=expand(exp)
    if exp.func!=Add:
        return []
    s=[]
    numterm=len(exp.args)
    listbin=listbinary(numterm)
    
    for i in listbin:
        exp1=exp
        for j in range(numterm):
            if i[j]==0:
                exp1=exp1-exp.args[j]
        exp0=exp-exp1
        exp0=expand(exp0)
        fexp0=factor(exp0)
        if ((fexp0!=exp0) or (exp0.func!=Add)) and (exp0!=0):
            s.append([fexp0,expand(exp1)])
    return s


def allformfac(exp):
#    return [exp]
    if exp.func!=Mul:
        return [exp]
    s=[]
    for term in exp.args:
        #print(str(exp)+":   "+str(term))
        listform=allform(term)
        
        #print(str(exp)+":   "+str(term)+ " listform")
        #print(listform)
        s.append(listform)
    #print("finallist")
    #print(s)
    pro = itertools.product(*s)
    listform2=[]
    for e in pro:
        t=e[0]
        for i in range(len(e)-1):
            t=t*e[i+1]
        listform2.append(t)
    
    return listform2
        
    

def allfrompart(part):
    listform=[]
    print("New partition:")
    print(part)
    print("------")
    left=part[0]
    right=part[1]
    leftlist=allformfac(left)
    #print(leftlist)
    rightlist=allform(right)
    listformpart=[x+y for x in leftlist for y in rightlist]
    for k in listformpart:
        if not(k in listform):
            listform.append(k)
            print(k)
    
    


def allform(exp):
    #print(exp)
    exp=expand(exp)
    listform=[exp]
    listpart=listpartition(exp)
    for part in listpart:
        left=part[0]
        right=part[1]
        leftlist=allformfac(left)
        #print(leftlist)
        rightlist=allform(right)
        listformpart=[x+y for x in leftlist for y in rightlist]
        for k in listformpart:
            if not(k in listform):
                listform.append(k)
                #print(k)
    
    return listform
        

def t_allform(exp,c):
    #print(exp)
    exp=expand(exp)
    listform=[exp]
    listpart=listpartition(exp)
    for part in listpart:
        left=part[0]
        right=part[1]
        leftlist=t_allformfac(left,c+1)
        #print(leftlist)
        rightlist=t_allform(right,c+1)
        listformpart=[x+y for x in leftlist for y in rightlist]
        for k in listformpart:
            if not(k in listform):
                listform.append(k)
                if c==0:
                    print(k)
    
    return listform


def t_allformfac(exp,c):
#    return [exp]
    if exp.func!=Mul:
        return [exp]
    s=[]
    for term in exp.args:
        #print(str(exp)+":   "+str(term))
        listform=t_allform(term,c+1)
        
        #print(str(exp)+":   "+str(term)+ " listform")
        #print(listform)
        s.append(listform)
    #print("finallist")
    #print(s)
    pro = itertools.product(*s)
    listform2=[]
    for e in pro:
        t=e[0]
        for i in range(len(e)-1):
            t=t*e[i+1]
        listform2.append(t)
    
    return listform2
        


        
#------------------------------------------------------##

def editinter(inter,prefix):
    dic={}
    c=0
    out=[]
    for e in inter:
        var,exp =e
        dic[var]=symbols(prefix+str(c))
        c=c+1
        for d in dic:
            exp=exp.subs(d,dic[d])
        out.append((dic[var],exp))
    return out
        
def fastsubs(s,e):
    e1=e
    for v in s:
        (a,b)=v
        e1=e1.subs(b,a)
    return e1
    
    
    
def get_intermediate(exp,prefix='t'):
    inter = t_cal_step2(exp)
    #prtlist(inter)
    s=[]
    count=0
    for e in inter:
        t,ex = e
        ex=fastsubs(s,ex)

        if len(ex.args)<=2:
            s.append((t,ex))
        else:
            
            if ex.func==Add:
                count=0
                ex0=ex.args[0]
                for a in range(len(ex.args)-1):
                    i=a+1
                    if i<len(ex.args)-1:
                        m=symbols(str(t)+"_"+str(count))
                    else:
                        m=t
                    ex0=ex0+ex.args[i]
                    s.append((m,ex0))
                    ex0=m
                    count=count+1
                    
            if (ex.func==Mul) and (ex.args[0].is_number == False):
                count=0
                ex0=ex.args[0]
                for a in range(len(ex.args)-1):
                    i=a+1
                    if i<len(ex.args)-1:
                        m=symbols(str(t)+"_"+str(count))
                    else:
                        m=t
                    ex0=ex0*ex.args[i]
                    s.append((m,ex0))
                    ex0=m
                    count=count+1 
            elif (ex.func==Mul):
                count=0
                ex0=ex.args[1]
                for a in range(len(ex.args)-2):
                    i=a+2
                    
                    m=symbols(str(t)+"_"+str(count))
                    
                    ex0=ex0*ex.args[i]
                    s.append((m,ex0))
                    ex0=m
                    count=count+1
                ex0=ex0*ex.args[0]
                m=t
                s.append((m,ex0))
                    

    return editinter(s,prefix)

def countops(exp,prefix='c'):
    inter=get_intermediate(exp,prefix)
    return inter, len(inter)


def absops(exp,inter,signlist):
    pair=[]
    noabs=[]
    out=[]
    for i in inter:
        var,ex=i
        if (ex.func==Mul) and (ex.args[0]==-1):
            pair.append(var)
            pair.append(ex.args[1])
    for i in inter:
        var,ex=i
        if (signlist[var]==1):
            noabs.append(var)
            if (signlist[var]==-1) and (var in pair):
                noabs.append(var)    

    #print(noabs)    
    c=0
    for t in exp.free_symbols:
        if t not in noabs:
            out.append((t,'fabs('+str(t)+')'))
            c=c+1
    return out

def final_countops(exp, asg,inter):
    sl=signlist(inter)
    absop=absops(exp,inter,sl)
    out=[]
    for e in absop:
        out.append(e)
    for e in asg:
        out.append(e)
    numop=len(out)
    return out, numop

def naieve_err(exp):
    inter, numop=count_ops_minus(exp)
    c=0
    for e in inter:      
        if '+' in e:
            c=c+4
        if '*' in e:
            c=c+8
    return c+numop

def editminuspower(s):
    out=[]
    for st in s:
        if st[len(st)-1]=='\n':
            st=st[:len(st)-1]
        st1=st.replace('+ -', '- ')
        l=len(st)
        if st[l-3:]=='**2':
            start=st.find('=')+1
            end=l-3
            term=st[start:end]
            st1=st[:start]+term+" * "+term
        if ('=-' in st) and ('*' not in st) and ('+' in st):
            sterm1=st.find('-')
            sterm2=st.find('+')
            term1=st[sterm1:sterm2]
            term2=st[sterm2+1:]
            st1=st[:sterm1]+term2+term1
        out.append(st1)
    return out

def count_ops_minus(exp):
    inter=get_intermediate(exp)
    dic={}
    rdic={}
    reuse={}
    for line in inter:
        var,exp=line
        reuse[var]=0
        
    s=[] #removable
    for line in inter:
        var,exp=line
        if (exp.func==Mul) and (exp.args[0]==-1):
            dic[var]=exp.args[1]
            rdic[exp.args[1]]=var
            s.append(var)
    for line in inter:
        var,exp=line
        if exp.func==Mul:
            if (exp.args[0] in s) and (exp.args[1] not in s):
                s.remove(exp.args[0])
            if (exp.args[1] in s) and (exp.args[0] not in s):
                s.remove(exp.args[1]) 
    for line in inter:
        var,exp=line
        if exp.func==Add:
            if (exp.args[0] in s) and (exp.args[1] in s):
                s.remove(exp.args[0])
            

    assignlist=[]
    for e in inter:
        var, exp =e
        if var not in s:
            st=str(exp)

            if exp.args[0] in s:
                st=st.replace(str(exp.args[0]),"-"+str(dic[exp.args[0]]))
            if exp.args[1] in s:
                st=st.replace(str(exp.args[1]),"-"+str(dic[exp.args[1]]))  
            st=str(var)+"="+st
            assignlist.append(st)
    
    numop=len(assignlist)
    return assignlist, numop


# def get_error(inter):
#     dic={}
#     out=[]
#     eps=symbols('eps')
#     for e in inter:
#         var, exp = e
#         fs=exp.free_symbols
#         for s in fs:
#             dic[s]=0
#     for e in inter:
#         var, exp = e
#         if exp.func==Add:
#             left=exp.args[0]
#             right=exp.args[1]            
#             if left.is_number:
#                 dic[left]=0
#             dic[var]=dic[left] + dic[right]+var*eps
            
#         if exp.func==Mul:
#             left=exp.args[0]
#             right=exp.args[1]
#             if left.is_number:
#                 if left==-1:
#                     dic[var]=dic[right]
#                 else:
#                     dic[var]=left*dic[right]+var*eps
#             else:    
#                 dic[var]=right*dic[left] + left*dic[right]+ dic[left]*dic[right]+ var*eps
            
#         if exp.func==Pow:
#             ex0=exp.args[0]
#             err=dic[ex0]
#             dic[var]=err*err+2*ex0*err+var*eps
            
#         out.append((var,dic[var]))
#         print((var,dic[var]))
#     return out


def get_error_sub(inter,subsmul,subsadd,epval,thrhld):
    dic={}
    out=[]
    #eps=symbols('eps')
    eps=epval
    for e in inter:
        var, exp = e
        fs=exp.free_symbols
        for s in fs:
            dic[s]=0
    for e in inter:
        var, exp = e
        if exp.func==Add:
            left=exp.args[0]
            right=exp.args[1]            
            if left.is_number:
                dic[left]=0
            dic[var]=dic[left] + dic[right]+var*eps
            
        if exp.func==Mul:
            left=exp.args[0]
            right=exp.args[1]
            if left.is_number:
                if left==-1:
                    dic[var]=dic[right]
                else:
                    if left<0:
                        temp=-left
                    else:
                        temp=left
                    dic[var]=temp*dic[right]+var*eps
            else:    
                dic[var]=right*dic[left] + left*dic[right]+ dic[left]*dic[right]+ var*eps
            
        if exp.func==Pow:
            ex0=exp.args[0]
            err=dic[ex0]
            dic[var]=err*err+2*ex0*err+var*eps
        numop=0    
        print()
        print(var)
        oexp=dic[var]
        adjust=oexp
        subs1=oexp
        if dic[var]!=0:          
            subs1, adjust, subs2, dic[var], numop =subsexp_adjustcoeff(oexp,subsmul,subsadd,thrhld)
        
        out.append((var,oexp,subs1,adjust,dic[var],numop))

        print(dic[var])
        print(str(numop))
    return out

def signlist(inter):
    dic={}
    for e in inter:
        var, exp = e
        fs=exp.free_symbols
        for s in fs:
            dic[s]=0
            
    for e in inter:
        var, exp = e
        if exp.func==Pow:
            dic[var]=1
        if exp.func==Mul:
            ex0=exp.args[0]
            ex1=exp.args[1]
            if ex0.is_number:
                if ex0 > 0:
                    dic[var]=dic[ex1]
                elif ex0 < 0:
                    dic[var]=-dic[ex1]
                else:
                    dic[var]=0
            else:
                dic[var]=dic[ex0]*dic[ex1]
        if exp.func==Add:
            ex0=exp.args[0]
            ex1=exp.args[1]
            
            if ex0.is_number:
                if (ex0 >= 0) and (dic[ex1]>0):
                    dic[var]=1
                elif (ex0 <= 0) and (dic[ex1]<0):
                    dic[var]=-1
                else:
                    dic[var]=0
            else:
                if (dic[ex0] > 0) and (dic[ex1]>0):
                    dic[var]=1
                elif (dic[ex0] < 0) and (dic[ex1]<0):
                    dic[var]=-1
                else:
                    dic[var]=0  
    return dic

def get_minus(inter):
    minus={}
    for e in inter:
        var, exp = e
        if (exp.func==Mul) and (exp.args[0]==-1):
            minus[var]=exp.args[1]
            minus[exp.args[1]]=var
    return minus

def subslist(inter,eps):
    outmul=[]
    outadd=[]
    sign = signlist(inter)
    #print(sign)
    #eps=symbols('eps')
    epsp2 = (1+2*eps)
    minus={}
    for e in inter:
        var, exp = e
        if (exp.func==Mul) and (exp.args[0]==-1):
            minus[var]=exp.args[1]
            minus[exp.args[1]]=var
           
    #prtlist(minus)
    for e in inter:
        var, exp = e
        if (exp.func==Mul) and (not(exp.args[0].is_number)):
            if (exp.args[0] in minus) and (exp.args[1] in minus):
                outmul.append((epsp2*var,exp))
                outmul.append((epsp2*var,exp.args[0]*minus[exp.args[1]]))
                outmul.append((epsp2*var,exp.args[1]*minus[exp.args[0]]))
                outmul.append((epsp2*var,minus[exp.args[1]]*minus[exp.args[0]]))
            elif (exp.args[0] in minus):
                outmul.append((epsp2*var,exp))
                outmul.append((epsp2*var,exp.args[1]*minus[exp.args[0]]))
            elif (exp.args[1] in minus):
                outmul.append((epsp2*var,exp))
                outmul.append((epsp2*var,exp.args[0]*minus[exp.args[1]]))                
            else:
                outmul.append((epsp2*var,exp))

        if ((exp.func==Mul) and (exp.args[0]!=-1)) and (exp.args[0].is_number):
            if (exp.args[1] in minus):
                outmul.append((epsp2*var*(1/abs(exp.args[0])),exp.args[1]))
                outmul.append((epsp2*var*(1/abs(exp.args[0])),minus[exp.args[1]]))                
            else:
                outmul.append((epsp2*var*(1/abs(exp.args[0])),exp.args[1]))

        if exp.func==Pow:
            if (exp.args[0] in minus):
                outmul.append((epsp2*var,exp))
                outmul.append((epsp2*var,exp.args[1]*minus[exp.args[0]]))                
            else:
                outmul.append((epsp2*var,exp))

        if exp.func==Add:
            if (exp.args[0] in minus) and (exp.args[1] in minus):
                outadd.append((var,epsp2*exp))
                outadd.append((var,epsp2*(exp.args[0]+minus[exp.args[1]])))
                outadd.append((var,epsp2*(exp.args[1]+minus[exp.args[0]])))
                outadd.append((var,epsp2*(minus[exp.args[1]]+minus[exp.args[0]])))
            elif (exp.args[0] in minus):
                outadd.append((var,epsp2*exp))
                outadd.append((var,epsp2*(exp.args[1]+minus[exp.args[0]])))
            elif (exp.args[1] in minus):
                outadd.append((var,epsp2*exp))
                outadd.append((var,epsp2*(exp.args[0]+minus[exp.args[1]])))                
            else:
                outadd.append((var,epsp2*exp))
            
            si0=0
            si1=0
            if exp.args[0].is_number:
                if exp.args[0]>0:
                    si0=1
                else:
                    si0=-1
            else:
                si0=sign[exp.args[0]]
          
            if exp.args[1].is_number:
                if exp.args[1]>0:
                    si1=1
                else:
                    si1=-1
            else:
                si1=sign[exp.args[1]]   
                 
            if si0*si1==1:
                outadd.append((exp.args[0],epsp2*var-exp.args[1]))
    return outmul, outadd

def revsubslist(inter,eps):
    outmul=[]
    outadd=[]
    sign = signlist(inter)
    #eps=symbols('eps')
    epsp2 = (1+2*eps)
    minus={}
    for e in inter:
        var, exp = e
        if (exp.func==Mul) and (exp.args[0]==-1):
            minus[var]=exp.args[1]
            minus[exp.args[1]]=var
           
    #prtlist(minus)
    for e in inter:
        var, exp = e
        if exp.func==Mul:
            if (exp.args[0] in minus) and (exp.args[1] in minus):
                outmul.append((var,epsp2*exp))
                outmul.append((var,epsp2*exp.args[0]*minus[exp.args[1]]))
                outmul.append((var,epsp2*exp.args[1]*minus[exp.args[0]]))
                outmul.append((var,epsp2*minus[exp.args[1]]*minus[exp.args[0]]))
            elif (exp.args[0] in minus):
                outmul.append((var,epsp2*exp))
                outmul.append((var,epsp2*exp.args[1]*minus[exp.args[0]]))
            elif (exp.args[1] in minus):
                outmul.append((var,epsp2*exp))
                outmul.append((var,epsp2*exp.args[0]*minus[exp.args[1]]))                
            else:
                outmul.append((var,epsp2*exp))

        if exp.func==Pow:
            if (exp.args[0] in minus):
                outmul.append((var,epsp2*exp))
                outmul.append((var,epsp2*exp.args[1]*minus[exp.args[0]]))                
            else:
                outmul.append((var,epsp2*exp))

        if exp.func==Add:
            si0=0
            si1=0
            if exp.args[0].is_number:
                if exp.args[0]>0:
                    si0=1
                else:
                    si0=-1
            else:
                si0=sign[exp.args[0]]
          
            if exp.args[1].is_number:
                if exp.args[1]>0:
                    si1=1
                else:
                    si1=-1
            else:
                si1=sign[exp.args[1]]   
                 
            if si0*si1==1:
                outmul.append((var,epsp2*(exp.args[0]+exp.args[1])))
    return outmul



def subsexp_mul(exp,subs):
    #exp=expand(exp)
    #print(exp)
    asg,numop=countops(exp)
    tnumop=0
    exp1=exp
    while tnumop<numop:
        asg,numop=countops(exp1)
        symbollist=exp1.free_symbols
        for e in subs:
            (a,b)=e
            if (len(b.args)==0): 
                if (b not in symbollist):
                    continue           
            elif (b.args[0] not in symbollist) and (b.args[1] not in symbollist)  :                               
                continue
            exp2=exp1.subs(b,a)
            asg,t1=countops(exp1)
            asg,t2=countops(exp2)
            if t2 < t1:
                exp1 = exp2
        asg,tnumop=countops(exp1)
    return exp1
        
def subsexp(exp,subsmul,subsadd):
    exp1=subsexp_mul(exp,subsmul)
    asg,numop=countops(exp1)
    tnumop=0
    c=0
    while tnumop<numop:
        asg,numop=countops(exp1)
        tnumop=numop
        exp2=exp1
        for e in subsadd:
            (a,b)=e
            exp2=exp1.subs(a,b)
            exp2=expand(exp2)
            exp2=subsexp_mul(exp2,subsmul)            
            asg,tnumop2=countops(exp2)
            #print(tnumop2)
            if tnumop2<tnumop:
                exp1=exp2
                tnumop=tnumop2
        asg,tnumop=countops(exp1)
        c=c+1
        #print(c)
    return exp1               

def subsexp_noex(exp,subsmul,subsadd):
    if exp.is_number:
        return exp
    exp1=subsexp_mul(exp,subsmul)
    asg,numop=countops(exp1)
    tnumop=0
    c=0
    while tnumop<numop:
        asg,numop=countops(exp1)
        tnumop=numop
        exp2=exp1
        symbollist=exp1.free_symbols
        for e in subsadd:
            (a,b)=e
            if a not in symbollist:
                continue
            exp2=exp1.subs(a,b)
            exp2=subsexp_mul(exp2,subsmul)
            asg,tnumop2=countops(exp2)
            if tnumop2<tnumop:
                exp1=exp2
                tnumop=tnumop2
        asg,tnumop=countops(exp1)
        c=c+1
    return exp1, tnumop 

def t_subsexp_noex(exp,subsmul,subsadd):
    if exp.is_number:
        return exp
    exp1=exp
    asg,numop=countops(exp1)
    tnumop=0
    c=0
    while tnumop<numop:
        asg,numop=countops(exp1)
        tnumop=numop
        exp2=exp1
        symbollist=exp1.free_symbols
        for e in subsadd:
            (a,b)=e
            if a not in symbollist:
                continue
            exp2=exp1.subs(a,b)
            exp2=subsexp_mul(exp2,subsmul)
            asg,tnumop2=countops(exp2)
            if tnumop2<tnumop:
                exp1=exp2
                tnumop=tnumop2
        asg,tnumop=countops(exp1)
        c=c+1
    return exp1, tnumop 

def doublesubsexp(exp,subsmul,subsadd,revsubs):
    if exp.is_number:
        return exp
    exp1=subsexp_mul(exp,subsmul)
    asg,numop=countops(exp1)
    tnumop=0
    while tnumop<numop:
        asg,numop=countops(exp1)
        tnumop=numop
        exp2=exp1
        symbollist=exp1.free_symbols
        count=0
        for e in revsubs:
            print(str(count)+ "    " + str(len(revsubs)))
            count=count+1
            (a,b)=e
            if a not in symbollist:
                continue
            exp2=exp1.subs(a,b)
            exp2,tnumop2=t_subsexp_noex(exp2,subsmul,subsadd)
            #asg,tnumop2=countops(exp2)
            if tnumop2<tnumop:
                exp1=exp2
                tnumop=tnumop2
        asg,tnumop=countops(exp1)
    return exp1, tnumop     
    


def subsexp_adjustcoeff(exp,subsmul,subsadd,thrhld):
    print(exp)
    subs1,numop1=subsexp_noex(exp,subsmul,subsadd)
    print(subs1)
    adjust = adjust_coeff_deep(subs1,thrhld)
    print(adjust)
    subs2,numop2=subsexp_noex(adjust,subsmul,subsadd)
    if (numop2<numop1) and (adjust!=subs2):
        return subs1, adjust, subs2, subs2, numop2
    else:
        return subs1, adjust, subs2, subs1, numop1
    
def subsexp_adjustcoeff_final(exp,subsmul,subsadd,thrhld):
    print(exp)
    subs1,numop1=subsexp_noex(exp,subsmul,subsadd)
    print(subs1)
    adjust = adjust_coeff_deep(subs1,thrhld)
    print(adjust)
    subs2,numop2=subsexp_noex(adjust,subsmul,subsadd)
    if numop2<=numop1:
        return subs1, adjust, subs2, subs2, numop2
    else:
        return subs1, adjust, subs2, subs1, numop1


def get_coeff(exp):
    s=[]
    for e in exp.args:
        if len(e.args)>0:
            if e.args[0].is_number:
                s.append(e.args[0])
    return s

def erase_coeff(exp):
    texp=1
    if len(exp.args)==2:
        return exp.args[1]
    for c in range(len(exp.args)-1):
        texp=texp*exp.args[c+1]
    return texp

def adjust_coeff_dic(s, thrhld):
    dic={}
    while len(s)>0:
        m=max(s)
        s1=[]
        for e in s:      
            if m-e<=thrhld:
                dic[e]=m
                s1.append(e)
        for u in s1:
            s.remove(u)
    return dic
        
def adjust_coeff(exp,thrhld):
    #print('Aj Co ')
    #print('origin: ' + str(exp))
    exp1=exp
    s=get_coeff(exp)
    #print('set co: ')
    #print(s)
    dic=adjust_coeff_dic(s, thrhld)
    #print('dic')
    #print(dic)
    for e in exp.args:
        if len(e.args)>0:
            if e.args[0].is_number:
                exp1=exp1-e+dic[e.args[0]]*erase_coeff(e)
    #print('out: '+str(exp1))
    return exp1    

def adjust_coeff_deep(exp,thrhld):
    if exp.is_number:
        return exp
    oexp=exp
    if exp.func==Add:
        oexp=0
        for e in exp.args:
            oexp=oexp+simplify(e)
        oexp=adjust_coeff(oexp,thrhld)
        ex1=0
        for e in oexp.args:
            ex1=ex1+adjust_coeff_deep(e,thrhld)
        oexp=ex1
    if exp.func==Mul:
        oexp=simplify(exp)
        ex1=1
        for e in oexp.args:
            ex1=ex1*adjust_coeff_deep(e,thrhld)
        oexp=ex1        
    if exp.func==Pow:
        #oexp=factor(exp)
        ex1=adjust_coeff_deep(exp.args[0],thrhld)
        oexp=ex1**exp.args[1]
    return simplify(oexp)
            
def adjust_coeff_deep_nofac(exp,thrhld):
    if exp.is_number:
        return exp
    oexp=exp
    if exp.func==Add:
        oexp=adjust_coeff(oexp,thrhld)
        ex1=0
        for e in oexp.args:
            ex1=ex1+adjust_coeff_deep_nofac(e,thrhld)
        oexp=ex1
    if exp.func==Mul:
        ex1=1
        for e in oexp.args:
            ex1=ex1*adjust_coeff_deep_nofac(e,thrhld)
        oexp=ex1        
    if exp.func==Pow:
        #oexp=factor(exp)
        ex1=adjust_coeff_deep_nofac(exp.args[0],thrhld)
        oexp=ex1**exp.args[1]
    return simplify(oexp)

def prtlist(s):
    for e in s:
        print(e)

def inter_des(inter):
    out=[]
    for e in inter:
        var, exp = e
        fexp=exp
        sexp=exp
        for e1 in inter:
            var1, exp1 =e1
            sexp=sexp.subs(var1,exp1)
        while sexp!=fexp:
            fexp = sexp
            for e1 in out:
                var1, exp1, subsexp =e1
                sexp=sexp.subs(var1,subsexp)
        out.append((var,exp,sexp))   
        print((var,exp,sexp))   
    return out
            

def writelist(s,file):
    for e in s:
        file.write(str(e)+'\n')  

def writecode(s,file):
    for e in s:
        if ('=' in e) and ('int' not in e):
            e=e.replace(' ','')
            if ('fabs' not in e) and (e[0]!='*'):
                e="REAL "+e
        if(str(e)!=''):
            file.write(str(e))
            file.write(';')
            file.write('\n')
        
def writeassignmentfile(s,file):
    for e in s:
        (l,r)=e
        file.write(str(l)+'='+ str(r)+'\n')

def writesubsfile(s,file):
    for e in s:
        (l,r)=e
        file.write(str(l)+' -> '+ str(r)+'\n')

def write_errorlist(s,inter,file):
    c=0
    for e in s:
        var, exp, subs1, adjust, subs2, numop =e
        v,ex=inter[c]
        c=c+1
        file.write("Inter step: " +str(var)+" = "+ str(ex)+"   # Operations:" + str(numop)+'\n')
        file.write("Before 1st subs exp: \n")
        file.write(str(exp)+'\n')
        file.write("After 1st subs exp: \n")
        file.write(str(subs1)+'\n') 
        file.write("After adjust coeff exp: \n")
        file.write(str(adjust)+'\n')
        file.write("After 2nd subs exp: \n")
        file.write(str(subs2)+'\n\n')      


                
            
def get_inter_dic(desinter,varlist,value):
    l=len(varlist)
    out={}
    for v in range(l):
        out[varlist[v]]=abs(value[v])
    for i in desinter:
        inter,texp,exp=i
        ex1=exp
        for v in range(l):
            ex1=ex1.subs(varlist[v],value[v])
        out[inter]=abs(ex1)
    return out

def evalexp(exp,desinter,varlist,value):
    interdic=get_inter_dic(desinter,varlist,value)
    out=exp
    for i in interdic:
        out=out.subs(i,interdic[i])
    return out
            
def randomdata(lenlist,numsam,minval,maxval,filename):
    random.seed()
    out=[]
    for i in range(numsam):
        row=[]
        for j in range(lenlist):
            x=random.uniform(minval,maxval)
            row.append(x)
        out.append(row)
    file=open(filename,'w')
    writelist(out,file)
    file.close()
    return out

def compare_func(funclist, funcname, varlist ,data, desinter, filename):
    numfunc=len(funclist)
    out=[]
    for i in range(numfunc):
        out.append(0)
    file=open(filename,'w')
    c=0
    for row in data:
        file.write("Row "+str(c)+":\n")
        c=c+1
        print(c)
        for i in range(numfunc):
            v=evalexp(funclist[i],desinter,varlist,row)
            out[i]=out[i]+v
            file.write(funcname[i]+":  "+str(v)+ ",         ")
        file.write('\n')
    file.close()
    for i in out:
        i=i/numfunc
        
    return out


#------------------------------------------------------##            

def firstitem(s,e):
    for i in range(len(s)):
        if i==2:
            return i

def randlist(l,p):
    s=[]
    q=0
    for i in range(l):
        k = random.randint(0,q)
        s.append(k)
        q=min(max(s)+1,p-1)
    return s

def goodpartition(s):
    for t in range(len(s)):
        exp=s[t]
        fexp=factor(exp)
        if ((fexp!=exp) or (exp.func!=Add)) and (exp!=0):
            return t
    return -1

def goodpartition2(s):
    c=0
    l=len(s)-1
    for t in range(len(s)):
        exp=s[t]
        fexp=factor(exp)
        if (fexp!=exp) :
            c=c+1
        if exp==0:
            l=l-1
    if c==l:
        return c
    return -1
    
def randpartition(exp,numpar):
    exp=expand(exp)
    if exp.func!=Add:
        return []
    s=[]
    numterm=len(exp.args)
    assign = randlist(numterm,numpar)

    for i in range(numpar):
        
        exp1=exp
        for j in range(numterm):
            if assign[j]!=i:
                exp1=exp1-exp.args[j]
        s.append(exp1)
    
    return s



def editoplist(oplist,prefix='tem'):
    out=[]
    appear={}
    subs={}
    count=0
    for op in oplist:
        (var,exp)=op
        texp=0
        while texp!=exp:
            texp=exp
            for s in subs:
                exp=exp.subs(s,subs[s])
                
        if exp in appear:
            subs[var]=appear[exp]
        else:
            appear[exp]=var
            out.append((var,exp))
    return editinter(out,prefix)
            
        
    

def fitness(exp,assign,numpar):
    par=partition(exp,assign,numpar)
    asglist=[]
    count=0
    tosum=[]
    for e in par:
        asg, no = countops(e,'p'+str(count))
        count=count+1
        for t in asg:
            asglist.append(t)
        if len(asg)!=0:
            #print(str(e))
            #print(asg)
            (a,b)=asg[len(asg)-1]
            tosum.append(a)
    last=tosum[0]        
    for c in range(len(tosum)-1):
        i=c+1
        m=symbols('psum'+str(c))
        asglist.append((m,last+tosum[i]))
        last=m            
        
    asglist=editoplist(asglist)   
    numop=len(asglist)
    return numop

def prtfitness(exp,assign,numpar):
    par=partition(exp,assign,numpar)
    asglist=[]
    count=0
    tosum=[]
    for e in par:
        asg, no = countops(e,'p'+str(count))
        count=count+1
        for t in asg:
            asglist.append(t)
        if len(asg)!=0:
            #print(str(e))
            #print(asg)
            (a,b)=asg[len(asg)-1]
            tosum.append(a)
    last=tosum[0]        
    for c in range(len(tosum)-1):
        i=c+1
        m=symbols('psum'+str(c))
        asglist.append((m,last+tosum[i]))
        last=m            
        
    asglist=editoplist(asglist)   
    numop=len(asglist)
    return asglist


def countops_partition(exp,assign,numpar):
    par=partition(exp,assign,numpar)
    asglist=[]
    parstring = ''
    count=0
    tosum=[]
    for e in par:
        parstring=parstring+'+('+str(simplify(e))+')\n'
        asg, no = countops(e,'p'+str(count))
        count=count+1
        for t in asg:
            asglist.append(t)
        if len(asg)!=0:
            (a,b)=asg[len(asg)-1]
            tosum.append(a)
    
    last=tosum[0]        
    for c in range(len(tosum)-1):
        i=c+1
        m=symbols('psum'+str(c))
        asglist.append((m,last+tosum[i]))
        last=m                  
            
    asglist=editoplist(asglist)   
    numop=len(asglist)
    parstring=parstring[1:]
    return parstring, asglist, numop             

def partition(exp,assign,numpar):
    if exp.func!=Add:
        return []
    s=[]
    numterm=len(exp.args)

    for i in range(numpar):
        
        exp1=exp
        for j in range(numterm):
            if assign[j]!=i:
                exp1=exp1-exp.args[j]
        s.append(exp1)
    
    return s   



def somepartitions(exp,numpar,numsam, maxtry=1000):
    pl=[]
    i=0
    while (len(pl)<numsam) and (i<maxtry):
        s =randpartition(exp,numpar)
        if s not in pl:
            pl.append(s)
        i=i+1
    return pl

def getneighborhood(assign,numpar):
    s=[]
    for i in range(len(assign)):
        for j in range(numpar):
            if j!=assign[i]:
                t=assign[:i]+[j]+assign[i+1:]
                s.append(t)
    return s
 
def bestpartition(listpar):
    numop=2000
    bestexp=0
    for par in listpar:
        exp,no = fitness(par)
        if no<numop:
            numop=no
            bestexp=exp
    return bestexp, numop

def tabusearchpar(start,exp,numpar,numsearch,listsize,cyclelen,secondcyclelen,file):
    sbest = start
    bestcand = start
    l=len(exp.args)
    bestfit = fitness(exp,sbest,numpar)
    tabulist=[]
    tabulist.insert(0,start)
    count=0
    cycle=0
    secondcycle=0
    lastbestfitcand=10000
    while count<numsearch:
        if cycle == cyclelen:
            secondcycle=secondcycle+1
            if secondcycle==secondcyclelen:
                numpar=numpar+1
                secondcycle=0
            bestcand=randlist(l,numpar)
            cycle=0
            
        sneighbor=getneighborhood(bestcand,numpar)
        for e in sneighbor:
            if e in tabulist:
                sneighbor.remove(e)
        bestcand=sneighbor[0]
        bestfitcand=fitness(exp,bestcand,numpar)
        for candidate in sneighbor:
            if (fitness(exp,candidate,numpar)<bestfitcand):
                bestcand=candidate
                bestfitcand = fitness(exp,bestcand,numpar)
        print(partition(exp,bestcand,numpar))
        #print(prtfitness(exp,bestcand,numpar))
        file.write(str(bestcand))
        if bestfitcand < bestfit:
            sbest=bestcand
            bestfit=bestfitcand
        if bestfitcand>=lastbestfitcand:
            cycle=cycle+1
        tabulist.insert(0,bestcand)
        if len(tabulist)>listsize:
            tabulist.pop()
        count=count+1
        lastbestfitcand=bestfitcand
        #prtlist(tabulist)
        print(str(count)+"  "+ str(cycle) + "   "+ str(secondcycle) + "   " + str(numpar) + "   "+ str(bestfitcand)+"   "+ str(bestfit))
        file.write("\nStep: " + str(count)+" RCycle: "+ str(cycle) + " NPCycle: "+ str(secondcycle) + " #Par: "+ str(numpar) + " BFCandicate: " + str(bestfitcand)+" BestFound: "+ str(bestfit)+"\n")
    return sbest, numpar

def tabusearch(exp,numpar,numsearch,listsize,cyclelen,secondcyclelen,file):
    if exp.func!=Add:
        asg,numop=countops(exp,'tem')
        return str(exp), asg, numop
    random.seed()    
    l=len(exp.args)
    k=randlist(l,numpar)
    result, numpar=tabusearchpar(k,exp,numpar,numsearch,listsize,cyclelen,secondcyclelen,file)
    print(result)
    bestexp, assign, bestnumop = countops_partition(exp,result,numpar)
    print('result:')
    print(bestexp)
    prtlist(assign)
    print(bestnumop) 
    return bestexp, assign, bestnumop

            
    


def FPSynthesis(exp,funcname):
    epval=2**(-53)
    thrhld=6*epval
    inter=get_intermediate(exp,'t')  
    numinter=len(inter)    
    subsmul,subsadd = subslist(inter,epval)
    revsubs = revsubslist(inter,epval)
    #minus = get_minus(inter)
    reportfilename=funcname+'/' + funcname+"_report.txt"
    outputfilename=funcname+'/' + funcname+"_program.txt"
    file=open(reportfilename,'w')
    outfile=open("temp.txt",'w')
    
    print('expression: ')
    print(exp)
    file.write('Expression: \n')
    file.write(str(exp))
    
    print('\n\nInter-steps:\n')
    desinter=inter_des(inter)
    prtlist(desinter)
    prtlist(inter)
    file.write('\n\nInter-steps:\n')
    writelist(desinter,file)
   
    file.write('\nMultiplication substitutions:\n')
    #writelist(subsmul,file)
    for e in subsmul:
        (a,b)=e
        if '-' not in str(b):
            file.write(str(b)+' -> ' + str(a) + '\n')
            
    file.write('\nAddition substitutions:\n')
    writesubsfile(subsadd,file)
    file.write('\nReverse substitutions:\n')
    writesubsfile(revsubs,file)
    
    errlist=get_error_sub(inter,subsmul,subsadd,epval,thrhld)
    t,oexp,subs1,adj, errexp,no = errlist[numinter-1]
    #print('\nError-list:')
    #prtlist(errlist)
    file.write('\nError-list:\n')
    write_errorlist(errlist,inter,file)
    
    print('\nDouble substitution:')
    errexp,no =doublesubsexp(errexp,subsmul,subsadd,revsubs)
    errexp=simplify(errexp)
    print(errexp)
    file.write('\nDouble substitution:\n')
    file.write(str(errexp)+'\n')
    file.write("#Operations: " + str(no) + "\n\n")
    print('\n')
    #subs1, adjust, subs2, errexp, no=subsexp_adjustcoeff_final(errexp,subsmul,subsadd,thrhld)
    errexp = adjust_coeff_deep(errexp,thrhld)
    print('\nFinal adjusting coefficient:')
    print(errexp)
    file.write('\nFinal adjusting coefficient:\n')
    file.write(str(errexp)+'\n')
    file.write("#Operations: " + str(no) + "\n\n")    
    
    numpar=2  
    numsearch=200
    listsize=40
    cyclelen=5
    secondcyclelen=5
    file.write('\nTabusearch parameters: ')
    file.write('\nNumber  of partitions: ' + str(numpar))
    file.write('\nNumber  of search steps: ' + str(numsearch))    
    file.write('\nList size: ' + str(listsize))  
    file.write('\nRandom cycle length: ' + str(cyclelen))
    file.write('\nNumpar cycle length: ' + str(secondcyclelen))
    file.write('\n\nTabusearch steps: \n')
    #errexp=expand(errexp)
    bestexp, tabuassign, bestnumop=tabusearch(errexp,numpar,numsearch,listsize,cyclelen,secondcyclelen,file)
    
    asg,numop=final_countops(errexp,tabuassign,inter) 
    print('\n\nFinal Error-bound expression: '+bestexp)
    file.write('\n\nFinal Error-bound expression: '+bestexp)
    print('\nNumber of operations: '+str(numop))
    file.write('\nNumber of operations: '+str(numop))
    print('\nCalculation steps: \n')
    file.write('\nCalculation steps: \n')
    prtlist(asg)
    writeassignmentfile(asg,file)
    file.close()
    
    
    outfile.write('\n//Calculation steps for input: \n')
    #outfile.write('\nint roundmode = fegetround()\n')
    assignlist, numop=count_ops_minus(exp)
    writelist(assignlist,outfile)
    outfile.write('\n//Calculation steps for Error-bound: \n')
    #outfile.write('\nfesetround(FE_UPWARD) \n')
    writeassignmentfile(asg,outfile)    
    inputres=assignlist[len(assignlist)-1]
    pos=inputres.find('=')
    inputres=inputres[:pos]
    (erbres,b) = asg[len(asg)-1]
    outfile.write('\n*result ='+ str(inputres)+'\n*error_bound = 1.0000000001* ' +str(erbres))
    #outfile.write('\nfesetround(roundmode) \n')
   
    outfile.close()
    outfile=open("temp.txt",'r')
    codelines=outfile.readlines()
    outfile.close()
    codelines=editminuspower(codelines)
    outfile=open(outputfilename,'w')
    funcpro = "void FPSyn_"+ funcname +"("
    freesym = []
    for fs in exp.free_symbols:
        freesym.append(str(fs))
    freesym.sort()
    for fsym in freesym:
        funcpro=funcpro+"REAL " + fsym + ","
    funcpro=funcpro+"REAL* result, REAL* error_bound){\n"
    outfile.write(funcpro)
    writecode(codelines,outfile)
    outfile.write("}\n")
    outfile.close()
    return errexp



























