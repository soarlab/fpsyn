# fpsyn
Synthesis of rigorous floating-point predicates

The tool is implemented as a python library in a single file FPSynLib.py

# Input:

To synthesize the program which calculate an input expression in floating-point arithmetic 
together with its error-bound, we create a python program with the following structure and
run it directly with python:

>import FPSynLib as FPS  
>from sympy import * 
>#input variables# =  symbols('#list of variables' name#')  
>expression  = #sympy expression#  
>errorexpression = FPS.FPSynthesis(expression,#expression name#)  


For example:  

>import FPSynLib as FPS  
>from sympy import *  
>ax,ay,bx,by,cx,cy =symbols('ax ay bx by cx cy')  
>orient2D=(ax-cx)*(by-cy)-(bx-cx)*(ay-cy)  
>errorexp=FPS.FPSynthesis(orient2D,"Orient2D")  

# Output

FPSyn will create a folder with the same name as the expression and three files in it  

- The *expressionname* _FPSyn.txt file: The final output function of FPSynth in C language
- The *expressionname* _program.txt file: The output function of FPSynth in C language, in which the right-hand-side of all assignments consist of one operation.
- The *expressionname* _report.txt file: Record the results of all steps that FPSynth synthesizes the output.


