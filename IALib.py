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
    if '+' in cl:
        op='+'
        term1=cl[cl.find('=')+1:cl.find('+')]
        term2=cl[cl.find('+')+1:cl.find(';')]
        out=(left,op,term1,term2)
    elif '*' in cl:
        op='*'
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
    return out

def readinput(filename):
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
        out.append(breakcodeline(c))
    return out

def readinputraw(filename):
    file=open(filename,'r')
    codelines=file.readlines()
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


def get_inter_bound(stlist):
    out=[]
    var={}
    for st in stlist:
        left,op,term1,term2=st
        if 't' not in term1:
            var[term1]=1
        if 't' not in term2:
            var[term2]=1

    for st in stlist:
        left,op,term1,term2=st
        if term1 in var:
            low1=term1
            high1=term1
        else:
            low1=term1+"_l"
            high1=term1+"_h"
            
        if term2 in var:
            low2=term2
            high2=term2
        else:
            low2=term2+"_l"
            high2=term2+"_h"
            
        low=left+"_l"
        high=left+"_h"
        
        if op=='+':
            if (term1 in var) and (term2 in var):
                lowst=low+"="+left+"-eps*"+left
                highst=high+"="+left+"+eps*"+left
            else:
                lowst=low+"="+low1+"+"+low2+"-eps*"+left
                highst=high+"="+high1+"+"+high2+"+eps*"+left               
              
        if op=='-':
            if (term1 in var) and (term2 in var):
                lowst=low+"="+left+"-eps*"+left
                highst=high+"="+left+"+eps*"+left
            else:
                lowst=low+"="+low1+"-"+high2+"-eps*"+left
                highst=high+"="+high1+"-"+low2+"+eps*"+left   
        
                

        if op=='*':
            if (term1 in var) and (term2 in var):
                lowst=low+"="+left+"-eps*"+left
                highst=high+"="+left+"+eps*"+left
            else:
                x1y1=low1+"*"+low2
                x1y2=low1+"*"+high2
                x2y1=high1+"*"+low2
                x2y2=high1+"*"+high2
                max2f="max("+x1y1+","+x1y2+")"
                max2l="max("+x2y1+","+x2y2+")"
                max4="max("+max2f+","+max2l+")"

                min2f="min("+x1y1+","+x1y2+")"
                min2l="min("+x2y1+","+x2y2+")"
                min4="min("+min2f+","+min2l+")"
                
                lowst=low+"="+min4+"-eps*"+left
                highst=high+"="+max4+"+eps*"+left  
        
        out.append(lowst)
        out.append(highst)
    return out

def isnum(s):
    if '.' in s:
        return True
    if '-' == s[0]:
        return True
    if s.isnumeric(): 
        return True
    return False

def get_inter_bound_binary(stlist):
    out=[]
    var={}
    defin=[]
    for st in stlist:
        left,op,term1,term2=st
        if 't' not in term1:
            var[term1]=1
        if 't' not in term2:
            var[term2]=1
    out.append("eps = pow(2,-53)")
    out.append('int inferr=0')
    out.append("int roundmode = fegetround()")
    out.append("fesetround(FE_UPWARD)")
    for st in stlist:
        
        left,op,term1,term2=st
        if term1 in var:
            low1=term1
            high1=term1
        else:
            low1=term1+"_l"
            high1=term1+"_h"
            
        if term2 in var:
            low2=term2
            high2=term2
        else:
            low2=term2+"_l"
            high2=term2+"_h"
            
        low=left+"_l"
        neglow=left+"_l_n"
        high=left+"_h"
        
        if (op in ["+","-","*","/"]) and (term1!='0'):            
            out.append(left+"_e_1=fabs("+left+")")
            out.append(left+"_e=eps*"+left+"_e_1")
        if op=='+':
            if (term1 in var) and (term2 in var):
                highst=high+"="+left+"+"+left+"_e"
                out.append(highst)
                
                neglowst=neglow+"="+left+"_e"+"-"+left
                out.append(neglowst)   
                
                #lowst=low+"=-"+neglowst
                #out.append(lowst) 
            else:
                neglow1=low1+"_n"
                neglow2=low2+"_n"
                if term1 in var:
                    if isnum(term1):
                        neglow1="-"+low1
                    else:
                        if low1 not in defin:
                            negt1st=low1+"_n=-"+low1
                            out.append(negt1st)
                            defin.append(low1)
                if term2 in var:
                    if isnum(term2):
                        neglow2="-"+low2
                    else:
                        if low2 not in defin:
                            negt2st=low2+"_n=-"+low2
                            out.append(negt2st)
                            defin.append(low2)
                
                highst1=high+"_1="+high1+"+"+high2
                out.append(highst1)
                highst=high+"="+high+"_1+"+ left+"_e"
                out.append(highst)      
                
             
                neglowst1=neglow+"_1="+neglow1+"+"+neglow2
                out.append(neglowst1)
                neglowst=neglow+"="+neglow+"_1+"+ left+"_e"
                out.append(neglowst) 
                            
                #lowst=low+"=-"+neglow
                #out.append(lowst) 
                
        if op=='-':
            if term1=='0':
                if (low2 not in defin) and (term2 in var):
                    negt2st=low2+"_n=-"+low2
                    out.append(negt2st)
                    defin.append(low2)
                out.append(high+"="+low2+"_n")    
                out.append(neglow+"="+high2)  
            elif (term1 in var) and (term2 in var):
                highst=high+"="+left+"+"+left+"_e"
                out.append(highst)
                
                neglowst=neglow+"="+left+"_e"+"-"+left
                out.append(neglowst)   
                
                #lowst=low+"=-"+neglow
                #out.append(lowst) 
            else:
                neglow1=low1+"_n"
                neglow2=low2+"_n"
                if term1 in var:
                    if isnum(term1):
                        neglow1="-"+low1
                    else:
                        if low1 not in defin:
                            negt1st=low1+"_n=-"+low1
                            out.append(negt1st)
                            defin.append(low1)
                if term2 in var:
                    if isnum(term2):
                        neglow2="-"+low2
                    else:
                        if low2 not in defin:
                            negt2st=low2+"_n=-"+low2
                            out.append(negt2st)
                            defin.append(low2) 
                    
                highst1=high+"_1="+high1+"+"+neglow2
                out.append(highst1)
                highst=high+"="+high+"_1+"+ left+"_e"
                out.append(highst)      
                
                neglowst1=neglow+"_1="+high2+"+"+neglow1
                out.append(neglowst1)
                neglowst=neglow+"="+neglow+"_1+"+ left+"_e"
                out.append(neglowst)    
                         
                #lowst=low+"=-"+neglow
                #out.append(lowst)                                        
                

        if op=='*':
            if (term1 in var) and (term2 in var):
                highst=high+"="+left+"+"+left+"_e"
                out.append(highst)
                
                neglowst=neglow+"="+left+"_e"+"-"+left
                out.append(neglowst) 
                
               
            else:
                neglow1=low1+"_n"
                neglow2=low2+"_n"
                if term1 in var:
                    if isnum(term1):
                        neglow1="-"+low1
                    else:
                        if low1 not in defin:
                            negt1st=low1+"_n=-"+low1
                            out.append(negt1st)
                            defin.append(low1)
                if term2 in var:
                    if isnum(term2):
                        neglow2="-"+low2
                    else:
                        if low2 not in defin:
                            negt2st=low2+"_n=-"+low2
                            out.append(negt2st)
                            defin.append(low2)  
                
                if isnum(term1):
                    if '-' != term1[0]:
                        highst1=high+"_1="+high1+"*"+high2
                        out.append(highst1)
                        highst=high+"="+high+"_1+"+ left+"_e"
                        out.append(highst)      
                        
                        neglowst1=neglow+"_1="+high1+"*"+neglow2
                        out.append(neglowst1)
                        neglowst=neglow+"="+neglow+"_1+"+ left+"_e"
                        out.append(neglowst)     
                    else:
                        highst1=high+"_1="+high1[1:]+"*"+neglow2
                        out.append(highst1)
                        highst=high+"="+high+"_1+"+ left+"_e"
                        out.append(highst)      
                        
                        neglowst1=neglow+"_1="+high1[1:]+"*"+high2
                        out.append(neglowst1)
                        neglowst=neglow+"="+neglow+"_1+"+ left+"_e"
                        out.append(neglowst)     
                else:                                                                                                                                            
                    out.append("REAL "+high+","+neglow)
                    ifst1="if("+high1+"<=0){"
                    out.append(ifst1)
                    ifst2="if("+high2+"<=0){"
                    out.append(ifst2)
                    neghigh1st=high1+"_n=-"+high1
                    out.append(neghigh1st)
                    neglowst1=neglow+"_1="+high1+"_n*"+high2
                    out.append(neglowst1)                
                    highst1=high+"_1="+neglow1+"*"+neglow2
                    out.append(highst1)
                    highst=high+"="+high+"_1+"+ left+"_e"
                    out.append("df"+highst)                    
                    neglowst=neglow+"="+neglow+"_1+"+ left+"_e"
                    out.append("df"+neglowst)
                    out.append("}")
                    
                    ifst3="else if("+neglow2+">=0){"
                    out.append(ifst3)
                    neglowst1=neglow+"_1="+neglow1+"*"+high2
                    out.append(neglowst1)                
                    highst1=high+"_1="+neglow1+"*"+neglow2
                    out.append(highst1)
                    highst=high+"="+high+"_1+"+ left+"_e"
                    out.append("df"+highst)                    
                    neglowst=neglow+"="+neglow+"_1+"+ left+"_e"
                    out.append("df"+neglowst)                    
                    out.append("}")
            
                    ifst4="else{"
                    out.append(ifst4)
                    neglowst1=neglow+"_1="+neglow1+"*"+high2
                    out.append(neglowst1)                                
                    neghigh1st=high1+"_n=-"+high1
                    out.append(neghigh1st)
                    highst1=high+"_1="+high1+"_n*"+neglow2
                    out.append(highst1)
                    highst=high+"="+high+"_1+"+ left+"_e"
                    out.append("df"+highst)                    
                    neglowst=neglow+"="+neglow+"_1+"+ left+"_e"
                    out.append("df"+neglowst)                    
                    out.append("}")
                    
                    out.append("}")
                    
                    ifst5="else if("+neglow1+">=0){"
                    out.append(ifst5)
                    ifst6="if("+high2+"<=0){"
                    out.append(ifst6)
                    neglowst1=neglow+"_1="+high1+"*"+neglow2
                    out.append(neglowst1)                
                    highst1=high+"_1="+neglow1+"*"+neglow2
                    out.append(highst1)
                    highst=high+"="+high+"_1+"+ left+"_e"
                    out.append("df"+highst)                    
                    neglowst=neglow+"="+neglow+"_1+"+ left+"_e"
                    out.append("df"+neglowst)
                    out.append("}")     
                    
                    ifst7="else if("+neglow2+">=0){"
                    out.append(ifst7)
                    neglowst11=neglow+"_1_1="+neglow1+"*"+high2
                    out.append(neglowst11)
                    neglowst12=neglow+"_1_2="+neglow2+"*"+high1
                    out.append(neglowst12)  
                    neglowst1=neglow+"_1=max("+neglow+"_1_1,"+neglow+"_1_2)"
                    out.append(neglowst1)     
                    
                    highst11=high+"_1_1="+neglow1+"*"+neglow2
                    out.append(highst11)
                    highst12=high+"_1_2="+high1+"*"+high2
                    out.append(highst12)  
                    highst1=high+"_1=max("+high+"_1_1,"+high+"_1_2)"
                    out.append(highst1)                   
                    highst=high+"="+high+"_1+"+ left+"_e"
                    out.append("df"+highst)                    
                    neglowst=neglow+"="+neglow+"_1+"+ left+"_e"
                    out.append("df"+neglowst)
                    out.append("}")     
                    
                    ifst8="else{"
                    out.append(ifst8)
                    neglowst1=neglow+"_1="+neglow1+"*"+high2
                    out.append(neglowst1)                                
                    highst1=high+"_1="+high1+"*"+high2
                    out.append(highst1)
                    highst=high+"="+high+"_1+"+ left+"_e"
                    out.append("df"+highst)                    
                    neglowst=neglow+"="+neglow+"_1+"+ left+"_e"
                    out.append("df"+neglowst)                    
                    out.append("}")
                    
                    out.append("}")                
    
                    ifst9="else{"
                    out.append(ifst9)
                    ifst10="if("+high2+"<=0){"
                    out.append(ifst10)
                    neglowst1=neglow+"_1="+neglow2+"*"+high1
                    out.append(neglowst1)                
                    neghigh2st=high2+"_n=-"+high2
                    out.append(neghigh2st)
                    highst1=high+"_1="+high2+"_n*"+neglow1
                    out.append(highst1)
                    highst=high+"="+high+"_1+"+ left+"_e"
                    out.append("df"+highst)                    
                    neglowst=neglow+"="+neglow+"_1+"+ left+"_e"
                    out.append("df"+neglowst)
                    out.append("}")
                    
                    ifst11="else if("+neglow2+">=0){"
                    out.append(ifst11)
                    neglowst1=neglow+"_1="+neglow2+"*"+high1
                    out.append(neglowst1)                
                    highst1=high+"_1="+high1+"*"+high2
                    out.append(highst1)
                    highst=high+"="+high+"_1+"+ left+"_e"
                    out.append("df"+highst)                    
                    neglowst=neglow+"="+neglow+"_1+"+ left+"_e"
                    out.append("df"+neglowst)                    
                    out.append("}")
            
                    ifst12="else{"
                    out.append(ifst12)
                    low1st=low1+"=-"+neglow1
                    out.append(low1st)                                
                    neglowst1=neglow+"_1="+low1+"*"+neglow2
                    out.append(neglowst1)                
                    highst1=high+"_1="+high1+"*"+high2
                    out.append(highst1)
                    highst=high+"="+high+"_1+"+ left+"_e"
                    out.append("df"+highst)                    
                    neglowst=neglow+"="+neglow+"_1+"+ left+"_e"
                    out.append("df"+neglowst)                    
                    out.append("}")
                    
                    out.append("}")


        if op=='/':            
            if (term1 in var) and (term2 in var):
                highst=high+"="+left+"+"+left+"_e"
                out.append(highst)
                
                neglowst=neglow+"="+left+"_e"+"-"+left
                out.append(neglowst) 
                
               
            else:
                neglow1=low1+"_n"
                neglow2=low2+"_n"
                if term1 in var:
                    if isnum(term1):
                        neglow1="-"+low1
                    else:
                        if low1 not in defin:
                            negt1st=low1+"_n=-"+low1
                            out.append(negt1st)
                            defin.append(low1)
                if term2 in var:
                    if isnum(term2):
                        neglow2="-"+low2
                    else:
                        if low2 not in defin:
                            negt2st=low2+"_n=-"+low2
                            out.append(negt2st)
                            defin.append(low2)                     
                else:
                    out.append('if('+ neglow2 + '*' + high2 + '>=0) inferr=1' )
                
                negt2invst=neglow2+"_inv=1/"+low2+'_n'
                out.append(negt2invst)   
                high2invst=high2+"_inv=1/"+high2
                out.append(high2invst)
                temp=high2
                high2=neglow2+"_inv"
                neglow2=temp+"_inv"
                
                if isnum(term1):
                    if '-' != term1[0]:
                        highst1=high+"_1="+high1+"*"+high2
                        out.append(highst1)
                        highst=high+"="+high+"_1+"+ left+"_e"
                        out.append(highst)      
                        
                        neglowst1=neglow+"_1="+high1+"*"+neglow2
                        out.append(neglowst1)
                        neglowst=neglow+"="+neglow+"_1+"+ left+"_e"
                        out.append(neglowst)     
                    else:
                        highst1=high+"_1="+high1[1:]+"*"+neglow2
                        out.append(highst1)
                        highst=high+"="+high+"_1+"+ left+"_e"
                        out.append(highst)      
                        
                        neglowst1=neglow+"_1="+high1[1:]+"*"+high2
                        out.append(neglowst1)
                        neglowst=neglow+"="+neglow+"_1+"+ left+"_e"
                        out.append(neglowst)     
                else:                                                                                                                                            
                    out.append("REAL "+high+","+neglow)
                    ifst1="if("+high1+"<=0){"
                    out.append(ifst1)
                    ifst2="if("+high2+"<=0){"
                    out.append(ifst2)
                    neghigh1st=high1+"_n=-"+high1
                    out.append(neghigh1st)
                    neglowst1=neglow+"_1="+high1+"_n*"+high2
                    out.append(neglowst1)                
                    highst1=high+"_1="+neglow1+"*"+neglow2
                    out.append(highst1)
                    highst=high+"="+high+"_1+"+ left+"_e"
                    out.append("df"+highst)                    
                    neglowst=neglow+"="+neglow+"_1+"+ left+"_e"
                    out.append("df"+neglowst)
                    out.append("}")
                    
                    ifst3="else if("+neglow2+">=0){"
                    out.append(ifst3)
                    neglowst1=neglow+"_1="+neglow1+"*"+high2
                    out.append(neglowst1)                
                    highst1=high+"_1="+neglow1+"*"+neglow2
                    out.append(highst1)
                    highst=high+"="+high+"_1+"+ left+"_e"
                    out.append("df"+highst)                    
                    neglowst=neglow+"="+neglow+"_1+"+ left+"_e"
                    out.append("df"+neglowst)                    
                    out.append("}")
            
                    ifst4="else{"
                    out.append(ifst4)
                    neglowst1=neglow+"_1="+neglow1+"*"+high2
                    out.append(neglowst1)                                
                    neghigh1st=high1+"_n=-"+high1
                    out.append(neghigh1st)
                    highst1=high+"_1="+high1+"_n*"+neglow2
                    out.append(highst1)
                    highst=high+"="+high+"_1+"+ left+"_e"
                    out.append("df"+highst)                    
                    neglowst=neglow+"="+neglow+"_1+"+ left+"_e"
                    out.append("df"+neglowst)                    
                    out.append("}")
                    
                    out.append("}")
                    
                    ifst5="else if("+neglow1+">=0){"
                    out.append(ifst5)
                    ifst6="if("+high2+"<=0){"
                    out.append(ifst6)
                    neglowst1=neglow+"_1="+high1+"*"+neglow2
                    out.append(neglowst1)                
                    highst1=high+"_1="+neglow1+"*"+neglow2
                    out.append(highst1)
                    highst=high+"="+high+"_1+"+ left+"_e"
                    out.append("df"+highst)                    
                    neglowst=neglow+"="+neglow+"_1+"+ left+"_e"
                    out.append("df"+neglowst)
                    out.append("}")     
                    
                    ifst7="else if("+neglow2+">=0){"
                    out.append(ifst7)
                    neglowst11=neglow+"_1_1="+neglow1+"*"+high2
                    out.append(neglowst11)
                    neglowst12=neglow+"_1_2="+neglow2+"*"+high1
                    out.append(neglowst12)  
                    neglowst1=neglow+"_1=max("+neglow+"_1_1,"+neglow+"_1_2)"
                    out.append(neglowst1)     
                    
                    highst11=high+"_1_1="+neglow1+"*"+neglow2
                    out.append(highst11)
                    highst12=high+"_1_2="+high1+"*"+high2
                    out.append(highst12)  
                    highst1=high+"_1=max("+high+"_1_1,"+high+"_1_2)"
                    out.append(highst1)                   
                    highst=high+"="+high+"_1+"+ left+"_e"
                    out.append("df"+highst)                    
                    neglowst=neglow+"="+neglow+"_1+"+ left+"_e"
                    out.append("df"+neglowst)
                    out.append("}")     
                    
                    ifst8="else{"
                    out.append(ifst8)
                    neglowst1=neglow+"_1="+neglow1+"*"+high2
                    out.append(neglowst1)                                
                    highst1=high+"_1="+high1+"*"+high2
                    out.append(highst1)
                    highst=high+"="+high+"_1+"+ left+"_e"
                    out.append("df"+highst)                    
                    neglowst=neglow+"="+neglow+"_1+"+ left+"_e"
                    out.append("df"+neglowst)                    
                    out.append("}")
                    
                    out.append("}")                
    
                    ifst9="else{"
                    out.append(ifst9)
                    ifst10="if("+high2+"<=0){"
                    out.append(ifst10)
                    neglowst1=neglow+"_1="+neglow2+"*"+high1
                    out.append(neglowst1)                
                    neghigh2st=high2+"_n=-"+high2
                    out.append(neghigh2st)
                    highst1=high+"_1="+high2+"_n*"+neglow1
                    out.append(highst1)
                    highst=high+"="+high+"_1+"+ left+"_e"
                    out.append("df"+highst)                    
                    neglowst=neglow+"="+neglow+"_1+"+ left+"_e"
                    out.append("df"+neglowst)
                    out.append("}")
                    
                    ifst11="else if("+neglow2+">=0){"
                    out.append(ifst11)
                    neglowst1=neglow+"_1="+neglow2+"*"+high1
                    out.append(neglowst1)                
                    highst1=high+"_1="+high1+"*"+high2
                    out.append(highst1)
                    highst=high+"="+high+"_1+"+ left+"_e"
                    out.append("df"+highst)                    
                    neglowst=neglow+"="+neglow+"_1+"+ left+"_e"
                    out.append("df"+neglowst)                    
                    out.append("}")
            
                    ifst12="else{"
                    out.append(ifst12)
                    low1st=low1+"=-"+neglow1
                    out.append(low1st)                                
                    neglowst1=neglow+"_1="+low1+"*"+neglow2
                    out.append(neglowst1)                
                    highst1=high+"_1="+high1+"*"+high2
                    out.append(highst1)
                    highst=high+"="+high+"_1+"+ left+"_e"
                    out.append("df"+highst)                    
                    neglowst=neglow+"="+neglow+"_1+"+ left+"_e"
                    out.append("df"+neglowst)                    
                    out.append("}")
                    
                    out.append("}")                

    
    out.append('fesetround(roundmode)')    
    out.append('*result = '+left)
    out.append('if (inferr==0) {')
    out.append('*hbound = '+high)
    out.append('*lbound = -'+neglow)
    out.append('}')

    out.append('else{')
    out.append('*hbound = 1/0')
    out.append('*lbound = 1/0')
    out.append('}')
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

def IAprogram(funcname):        
    filename=funcname+'/'+funcname+ "_program.txt"
    file=open(funcname+'/'+funcname+"_IA.txt",'w')
    out=readinputraw(filename)
    for o in out:
        file.write(o)
    
    file.write("//Calculate bound\n")    
    out=get_inter_bound_binary(readinput(filename))
    writecode(out,file)
    
    file.write('}\n')
    file.close()


#IAprogram('Intersection2D')






            

        