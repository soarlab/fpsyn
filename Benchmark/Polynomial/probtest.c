#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <fenv.h>
#include <gmp.h>
#include <mpfr.h>

#define REAL double

#define max(a,b) \
   ({ __typeof__ (a) _a = (a); \
       __typeof__ (b) _b = (b); \
     _a > _b ? _a : _b; })

#define min(a,b) \
   ({ __typeof__ (a) _a = (a); \
       __typeof__ (b) _b = (b); \
     _a < _b ? _a : _b; })

#define Absolute(a)  ((a) >= 0.0 ? (a) : -(a))

REAL randfl_uniReal(a)
REAL a;
{
    return (REAL)rand()/((REAL)RAND_MAX/a);
}

REAL randfl(a)
int a;
{
    int r = (rand() % a)-(a/2);
	REAL t = (REAL)rand()/((REAL)RAND_MAX);
    return (t - (1/2))*pow(2,r);
}
//FPSyn Function
void FPSyn_Polynomial(REAL x,REAL* result, REAL* error_bound){
REAL t0 = x*x;
REAL t1 = 0.5*t0;
REAL t2 = t0*t0;
REAL t3 = t0*t2;
REAL t4 = 0.00138888888888889*t3;
REAL t5 = t2*t2;
REAL t7 = 2.75573192239859e-6*(t5*x);
REAL t9 = 0.00833333333333333*(t2*x);
REAL t11 = 2.75573192239859e-7*(t0*t5);
REAL t12 = 0.0416666666666667*t2;
REAL t14 = 0.166666666666667*(t0*x);
REAL t15 = 2.48015873015873e-5*t5;
REAL t17 = 0.000198412698412698*(t3*x);
REAL t21 = t14+(t12+(t11+(t1+1)));
REAL t23 = t17+(t15+t21);
REAL t25 = (t23+t4)+t7;
REAL t27 = (t25+t9)+x;
*result = t27;
*error_bound = 1.0000000001*((1.33226762955019e-15*((Absolute(t7))+(t4+((Absolute(t17))+(t11+t15)))))+(6.66133814775094e-16*((Absolute(t9))+((Absolute(t27))+((Absolute(t25))+((Absolute(t23))+((Absolute(t21))+((Absolute(t14))+(t1+t12)))))))));
}//IA Function
void IA_Polynomial(REAL x,REAL* result, REAL* hbound, REAL* lbound){
REAL t0=x*x;
REAL t1=0.5*t0;
REAL t2=t0*t0;
REAL t3=t0*t2;
REAL t4=0.00138888888888889*t3;
REAL t5=t2*t2;
REAL t6=t5*x;
REAL t7=2.75573192239859e-6*t6;
REAL t8=t2*x;
REAL t9=0.00833333333333333*t8;
REAL t10=t0*t5;
REAL t11=2.75573192239859e-7*t10;
REAL t12=0.0416666666666667*t2;
REAL t13=t0*x;
REAL t14=0.166666666666667*t13;
REAL t15=2.48015873015873e-5*t5;
REAL t16=t3*x;
REAL t17=0.000198412698412698*t16;
REAL t18=t1+1;
REAL t19=t11+t18;
REAL t20=t12+t19;
REAL t21=t14+t20;
REAL t22=t15+t21;
REAL t23=t17+t22;
REAL t24=t23+t4;
REAL t25=t24+t7;
REAL t26=t25+t9;
REAL t27=t26+x;
//Calculate bound
REAL eps=pow(2,-53);
int inferr=0;
int roundmode = fegetround();
fesetround(FE_UPWARD);
REAL t0_e_1=fabs(t0);
REAL t0_e=eps*t0_e_1;
REAL t0_h=t0+t0_e;
REAL t0_l_n=t0_e-t0;
REAL t1_e_1=fabs(t1);
REAL t1_e=eps*t1_e_1;
REAL t1_h_1=0.5*t0_h;
REAL t1_h=t1_h_1+t1_e;
REAL t1_l_n_1=0.5*t0_l_n;
REAL t1_l_n=t1_l_n_1+t1_e;
REAL t2_e_1=fabs(t2);
REAL t2_e=eps*t2_e_1;
REAL  t2_h, t2_l_n;
if(t0_h<=0){
if(t0_h<=0){
REAL t0_h_n=-t0_h;
REAL t2_l_n_1=t0_h_n*t0_h;
REAL t2_h_1=t0_l_n*t0_l_n;
t2_h=t2_h_1+t2_e;
t2_l_n=t2_l_n_1+t2_e;
}
else if(t0_l_n>=0){
REAL t2_l_n_1=t0_l_n*t0_h;
REAL t2_h_1=t0_l_n*t0_l_n;
t2_h=t2_h_1+t2_e;
t2_l_n=t2_l_n_1+t2_e;
}
else{
REAL t2_l_n_1=t0_l_n*t0_h;
REAL t0_h_n=-t0_h;
REAL t2_h_1=t0_h_n*t0_l_n;
t2_h=t2_h_1+t2_e;
t2_l_n=t2_l_n_1+t2_e;
}
}
else if(t0_l_n>=0){
if(t0_h<=0){
REAL t2_l_n_1=t0_h*t0_l_n;
REAL t2_h_1=t0_l_n*t0_l_n;
t2_h=t2_h_1+t2_e;
t2_l_n=t2_l_n_1+t2_e;
}
else if(t0_l_n>=0){
REAL t2_l_n_1_1=t0_l_n*t0_h;
REAL t2_l_n_1_2=t0_l_n*t0_h;
REAL t2_l_n_1=max(t2_l_n_1_1,t2_l_n_1_2);
REAL t2_h_1_1=t0_l_n*t0_l_n;
REAL t2_h_1_2=t0_h*t0_h;
REAL t2_h_1=max(t2_h_1_1,t2_h_1_2);
t2_h=t2_h_1+t2_e;
t2_l_n=t2_l_n_1+t2_e;
}
else{
REAL t2_l_n_1=t0_l_n*t0_h;
REAL t2_h_1=t0_h*t0_h;
t2_h=t2_h_1+t2_e;
t2_l_n=t2_l_n_1+t2_e;
}
}
else{
if(t0_h<=0){
REAL t2_l_n_1=t0_l_n*t0_h;
REAL t0_h_n=-t0_h;
REAL t2_h_1=t0_h_n*t0_l_n;
t2_h=t2_h_1+t2_e;
t2_l_n=t2_l_n_1+t2_e;
}
else if(t0_l_n>=0){
REAL t2_l_n_1=t0_l_n*t0_h;
REAL t2_h_1=t0_h*t0_h;
t2_h=t2_h_1+t2_e;
t2_l_n=t2_l_n_1+t2_e;
}
else{
REAL t0_l=-t0_l_n;
REAL t2_l_n_1=t0_l*t0_l_n;
REAL t2_h_1=t0_h*t0_h;
t2_h=t2_h_1+t2_e;
t2_l_n=t2_l_n_1+t2_e;
}
}
REAL t3_e_1=fabs(t3);
REAL t3_e=eps*t3_e_1;
REAL  t3_h, t3_l_n;
if(t0_h<=0){
if(t2_h<=0){
REAL t0_h_n=-t0_h;
REAL t3_l_n_1=t0_h_n*t2_h;
REAL t3_h_1=t0_l_n*t2_l_n;
t3_h=t3_h_1+t3_e;
t3_l_n=t3_l_n_1+t3_e;
}
else if(t2_l_n>=0){
REAL t3_l_n_1=t0_l_n*t2_h;
REAL t3_h_1=t0_l_n*t2_l_n;
t3_h=t3_h_1+t3_e;
t3_l_n=t3_l_n_1+t3_e;
}
else{
REAL t3_l_n_1=t0_l_n*t2_h;
REAL t0_h_n=-t0_h;
REAL t3_h_1=t0_h_n*t2_l_n;
t3_h=t3_h_1+t3_e;
t3_l_n=t3_l_n_1+t3_e;
}
}
else if(t0_l_n>=0){
if(t2_h<=0){
REAL t3_l_n_1=t0_h*t2_l_n;
REAL t3_h_1=t0_l_n*t2_l_n;
t3_h=t3_h_1+t3_e;
t3_l_n=t3_l_n_1+t3_e;
}
else if(t2_l_n>=0){
REAL t3_l_n_1_1=t0_l_n*t2_h;
REAL t3_l_n_1_2=t2_l_n*t0_h;
REAL t3_l_n_1=max(t3_l_n_1_1,t3_l_n_1_2);
REAL t3_h_1_1=t0_l_n*t2_l_n;
REAL t3_h_1_2=t0_h*t2_h;
REAL t3_h_1=max(t3_h_1_1,t3_h_1_2);
t3_h=t3_h_1+t3_e;
t3_l_n=t3_l_n_1+t3_e;
}
else{
REAL t3_l_n_1=t0_l_n*t2_h;
REAL t3_h_1=t0_h*t2_h;
t3_h=t3_h_1+t3_e;
t3_l_n=t3_l_n_1+t3_e;
}
}
else{
if(t2_h<=0){
REAL t3_l_n_1=t2_l_n*t0_h;
REAL t2_h_n=-t2_h;
REAL t3_h_1=t2_h_n*t0_l_n;
t3_h=t3_h_1+t3_e;
t3_l_n=t3_l_n_1+t3_e;
}
else if(t2_l_n>=0){
REAL t3_l_n_1=t2_l_n*t0_h;
REAL t3_h_1=t0_h*t2_h;
t3_h=t3_h_1+t3_e;
t3_l_n=t3_l_n_1+t3_e;
}
else{
REAL t0_l=-t0_l_n;
REAL t3_l_n_1=t0_l*t2_l_n;
REAL t3_h_1=t0_h*t2_h;
t3_h=t3_h_1+t3_e;
t3_l_n=t3_l_n_1+t3_e;
}
}
REAL t4_e_1=fabs(t4);
REAL t4_e=eps*t4_e_1;
REAL t4_h_1=0.00138888888888889*t3_h;
REAL t4_h=t4_h_1+t4_e;
REAL t4_l_n_1=0.00138888888888889*t3_l_n;
REAL t4_l_n=t4_l_n_1+t4_e;
REAL t5_e_1=fabs(t5);
REAL t5_e=eps*t5_e_1;
REAL  t5_h, t5_l_n;
if(t2_h<=0){
if(t2_h<=0){
REAL t2_h_n=-t2_h;
REAL t5_l_n_1=t2_h_n*t2_h;
REAL t5_h_1=t2_l_n*t2_l_n;
t5_h=t5_h_1+t5_e;
t5_l_n=t5_l_n_1+t5_e;
}
else if(t2_l_n>=0){
REAL t5_l_n_1=t2_l_n*t2_h;
REAL t5_h_1=t2_l_n*t2_l_n;
t5_h=t5_h_1+t5_e;
t5_l_n=t5_l_n_1+t5_e;
}
else{
REAL t5_l_n_1=t2_l_n*t2_h;
REAL t2_h_n=-t2_h;
REAL t5_h_1=t2_h_n*t2_l_n;
t5_h=t5_h_1+t5_e;
t5_l_n=t5_l_n_1+t5_e;
}
}
else if(t2_l_n>=0){
if(t2_h<=0){
REAL t5_l_n_1=t2_h*t2_l_n;
REAL t5_h_1=t2_l_n*t2_l_n;
t5_h=t5_h_1+t5_e;
t5_l_n=t5_l_n_1+t5_e;
}
else if(t2_l_n>=0){
REAL t5_l_n_1_1=t2_l_n*t2_h;
REAL t5_l_n_1_2=t2_l_n*t2_h;
REAL t5_l_n_1=max(t5_l_n_1_1,t5_l_n_1_2);
REAL t5_h_1_1=t2_l_n*t2_l_n;
REAL t5_h_1_2=t2_h*t2_h;
REAL t5_h_1=max(t5_h_1_1,t5_h_1_2);
t5_h=t5_h_1+t5_e;
t5_l_n=t5_l_n_1+t5_e;
}
else{
REAL t5_l_n_1=t2_l_n*t2_h;
REAL t5_h_1=t2_h*t2_h;
t5_h=t5_h_1+t5_e;
t5_l_n=t5_l_n_1+t5_e;
}
}
else{
if(t2_h<=0){
REAL t5_l_n_1=t2_l_n*t2_h;
REAL t2_h_n=-t2_h;
REAL t5_h_1=t2_h_n*t2_l_n;
t5_h=t5_h_1+t5_e;
t5_l_n=t5_l_n_1+t5_e;
}
else if(t2_l_n>=0){
REAL t5_l_n_1=t2_l_n*t2_h;
REAL t5_h_1=t2_h*t2_h;
t5_h=t5_h_1+t5_e;
t5_l_n=t5_l_n_1+t5_e;
}
else{
REAL t2_l=-t2_l_n;
REAL t5_l_n_1=t2_l*t2_l_n;
REAL t5_h_1=t2_h*t2_h;
t5_h=t5_h_1+t5_e;
t5_l_n=t5_l_n_1+t5_e;
}
}
REAL t6_e_1=fabs(t6);
REAL t6_e=eps*t6_e_1;
REAL x_n=-x;
REAL  t6_h, t6_l_n;
if(t5_h<=0){
if(x<=0){
REAL t5_h_n=-t5_h;
REAL t6_l_n_1=t5_h_n*x;
REAL t6_h_1=t5_l_n*x_n;
t6_h=t6_h_1+t6_e;
t6_l_n=t6_l_n_1+t6_e;
}
else if(x_n>=0){
REAL t6_l_n_1=t5_l_n*x;
REAL t6_h_1=t5_l_n*x_n;
t6_h=t6_h_1+t6_e;
t6_l_n=t6_l_n_1+t6_e;
}
else{
REAL t6_l_n_1=t5_l_n*x;
REAL t5_h_n=-t5_h;
REAL t6_h_1=t5_h_n*x_n;
t6_h=t6_h_1+t6_e;
t6_l_n=t6_l_n_1+t6_e;
}
}
else if(t5_l_n>=0){
if(x<=0){
REAL t6_l_n_1=t5_h*x_n;
REAL t6_h_1=t5_l_n*x_n;
t6_h=t6_h_1+t6_e;
t6_l_n=t6_l_n_1+t6_e;
}
else if(x_n>=0){
REAL t6_l_n_1_1=t5_l_n*x;
REAL t6_l_n_1_2=x_n*t5_h;
REAL t6_l_n_1=max(t6_l_n_1_1,t6_l_n_1_2);
REAL t6_h_1_1=t5_l_n*x_n;
REAL t6_h_1_2=t5_h*x;
REAL t6_h_1=max(t6_h_1_1,t6_h_1_2);
t6_h=t6_h_1+t6_e;
t6_l_n=t6_l_n_1+t6_e;
}
else{
REAL t6_l_n_1=t5_l_n*x;
REAL t6_h_1=t5_h*x;
t6_h=t6_h_1+t6_e;
t6_l_n=t6_l_n_1+t6_e;
}
}
else{
if(x<=0){
REAL t6_l_n_1=x_n*t5_h;
REAL x_n=-x;
REAL t6_h_1=x_n*t5_l_n;
t6_h=t6_h_1+t6_e;
t6_l_n=t6_l_n_1+t6_e;
}
else if(x_n>=0){
REAL t6_l_n_1=x_n*t5_h;
REAL t6_h_1=t5_h*x;
t6_h=t6_h_1+t6_e;
t6_l_n=t6_l_n_1+t6_e;
}
else{
REAL t5_l=-t5_l_n;
REAL t6_l_n_1=t5_l*x_n;
REAL t6_h_1=t5_h*x;
t6_h=t6_h_1+t6_e;
t6_l_n=t6_l_n_1+t6_e;
}
}
REAL t7_e_1=fabs(t7);
REAL t7_e=eps*t7_e_1;
REAL t7_h_1=2.75573192239859e-6*t6_h;
REAL t7_h=t7_h_1+t7_e;
REAL t7_l_n_1=2.75573192239859e-6*t6_l_n;
REAL t7_l_n=t7_l_n_1+t7_e;
REAL t8_e_1=fabs(t8);
REAL t8_e=eps*t8_e_1;
REAL  t8_h, t8_l_n;
if(t2_h<=0){
if(x<=0){
REAL t2_h_n=-t2_h;
REAL t8_l_n_1=t2_h_n*x;
REAL t8_h_1=t2_l_n*x_n;
t8_h=t8_h_1+t8_e;
t8_l_n=t8_l_n_1+t8_e;
}
else if(x_n>=0){
REAL t8_l_n_1=t2_l_n*x;
REAL t8_h_1=t2_l_n*x_n;
t8_h=t8_h_1+t8_e;
t8_l_n=t8_l_n_1+t8_e;
}
else{
REAL t8_l_n_1=t2_l_n*x;
REAL t2_h_n=-t2_h;
REAL t8_h_1=t2_h_n*x_n;
t8_h=t8_h_1+t8_e;
t8_l_n=t8_l_n_1+t8_e;
}
}
else if(t2_l_n>=0){
if(x<=0){
REAL t8_l_n_1=t2_h*x_n;
REAL t8_h_1=t2_l_n*x_n;
t8_h=t8_h_1+t8_e;
t8_l_n=t8_l_n_1+t8_e;
}
else if(x_n>=0){
REAL t8_l_n_1_1=t2_l_n*x;
REAL t8_l_n_1_2=x_n*t2_h;
REAL t8_l_n_1=max(t8_l_n_1_1,t8_l_n_1_2);
REAL t8_h_1_1=t2_l_n*x_n;
REAL t8_h_1_2=t2_h*x;
REAL t8_h_1=max(t8_h_1_1,t8_h_1_2);
t8_h=t8_h_1+t8_e;
t8_l_n=t8_l_n_1+t8_e;
}
else{
REAL t8_l_n_1=t2_l_n*x;
REAL t8_h_1=t2_h*x;
t8_h=t8_h_1+t8_e;
t8_l_n=t8_l_n_1+t8_e;
}
}
else{
if(x<=0){
REAL t8_l_n_1=x_n*t2_h;
REAL x_n=-x;
REAL t8_h_1=x_n*t2_l_n;
t8_h=t8_h_1+t8_e;
t8_l_n=t8_l_n_1+t8_e;
}
else if(x_n>=0){
REAL t8_l_n_1=x_n*t2_h;
REAL t8_h_1=t2_h*x;
t8_h=t8_h_1+t8_e;
t8_l_n=t8_l_n_1+t8_e;
}
else{
REAL t2_l=-t2_l_n;
REAL t8_l_n_1=t2_l*x_n;
REAL t8_h_1=t2_h*x;
t8_h=t8_h_1+t8_e;
t8_l_n=t8_l_n_1+t8_e;
}
}
REAL t9_e_1=fabs(t9);
REAL t9_e=eps*t9_e_1;
REAL t9_h_1=0.00833333333333333*t8_h;
REAL t9_h=t9_h_1+t9_e;
REAL t9_l_n_1=0.00833333333333333*t8_l_n;
REAL t9_l_n=t9_l_n_1+t9_e;
REAL t10_e_1=fabs(t10);
REAL t10_e=eps*t10_e_1;
REAL  t10_h, t10_l_n;
if(t0_h<=0){
if(t5_h<=0){
REAL t0_h_n=-t0_h;
REAL t10_l_n_1=t0_h_n*t5_h;
REAL t10_h_1=t0_l_n*t5_l_n;
t10_h=t10_h_1+t10_e;
t10_l_n=t10_l_n_1+t10_e;
}
else if(t5_l_n>=0){
REAL t10_l_n_1=t0_l_n*t5_h;
REAL t10_h_1=t0_l_n*t5_l_n;
t10_h=t10_h_1+t10_e;
t10_l_n=t10_l_n_1+t10_e;
}
else{
REAL t10_l_n_1=t0_l_n*t5_h;
REAL t0_h_n=-t0_h;
REAL t10_h_1=t0_h_n*t5_l_n;
t10_h=t10_h_1+t10_e;
t10_l_n=t10_l_n_1+t10_e;
}
}
else if(t0_l_n>=0){
if(t5_h<=0){
REAL t10_l_n_1=t0_h*t5_l_n;
REAL t10_h_1=t0_l_n*t5_l_n;
t10_h=t10_h_1+t10_e;
t10_l_n=t10_l_n_1+t10_e;
}
else if(t5_l_n>=0){
REAL t10_l_n_1_1=t0_l_n*t5_h;
REAL t10_l_n_1_2=t5_l_n*t0_h;
REAL t10_l_n_1=max(t10_l_n_1_1,t10_l_n_1_2);
REAL t10_h_1_1=t0_l_n*t5_l_n;
REAL t10_h_1_2=t0_h*t5_h;
REAL t10_h_1=max(t10_h_1_1,t10_h_1_2);
t10_h=t10_h_1+t10_e;
t10_l_n=t10_l_n_1+t10_e;
}
else{
REAL t10_l_n_1=t0_l_n*t5_h;
REAL t10_h_1=t0_h*t5_h;
t10_h=t10_h_1+t10_e;
t10_l_n=t10_l_n_1+t10_e;
}
}
else{
if(t5_h<=0){
REAL t10_l_n_1=t5_l_n*t0_h;
REAL t5_h_n=-t5_h;
REAL t10_h_1=t5_h_n*t0_l_n;
t10_h=t10_h_1+t10_e;
t10_l_n=t10_l_n_1+t10_e;
}
else if(t5_l_n>=0){
REAL t10_l_n_1=t5_l_n*t0_h;
REAL t10_h_1=t0_h*t5_h;
t10_h=t10_h_1+t10_e;
t10_l_n=t10_l_n_1+t10_e;
}
else{
REAL t0_l=-t0_l_n;
REAL t10_l_n_1=t0_l*t5_l_n;
REAL t10_h_1=t0_h*t5_h;
t10_h=t10_h_1+t10_e;
t10_l_n=t10_l_n_1+t10_e;
}
}
REAL t11_e_1=fabs(t11);
REAL t11_e=eps*t11_e_1;
REAL t11_h_1=2.75573192239859e-7*t10_h;
REAL t11_h=t11_h_1+t11_e;
REAL t11_l_n_1=2.75573192239859e-7*t10_l_n;
REAL t11_l_n=t11_l_n_1+t11_e;
REAL t12_e_1=fabs(t12);
REAL t12_e=eps*t12_e_1;
REAL t12_h_1=0.0416666666666667*t2_h;
REAL t12_h=t12_h_1+t12_e;
REAL t12_l_n_1=0.0416666666666667*t2_l_n;
REAL t12_l_n=t12_l_n_1+t12_e;
REAL t13_e_1=fabs(t13);
REAL t13_e=eps*t13_e_1;
REAL  t13_h, t13_l_n;
if(t0_h<=0){
if(x<=0){
REAL t0_h_n=-t0_h;
REAL t13_l_n_1=t0_h_n*x;
REAL t13_h_1=t0_l_n*x_n;
t13_h=t13_h_1+t13_e;
t13_l_n=t13_l_n_1+t13_e;
}
else if(x_n>=0){
REAL t13_l_n_1=t0_l_n*x;
REAL t13_h_1=t0_l_n*x_n;
t13_h=t13_h_1+t13_e;
t13_l_n=t13_l_n_1+t13_e;
}
else{
REAL t13_l_n_1=t0_l_n*x;
REAL t0_h_n=-t0_h;
REAL t13_h_1=t0_h_n*x_n;
t13_h=t13_h_1+t13_e;
t13_l_n=t13_l_n_1+t13_e;
}
}
else if(t0_l_n>=0){
if(x<=0){
REAL t13_l_n_1=t0_h*x_n;
REAL t13_h_1=t0_l_n*x_n;
t13_h=t13_h_1+t13_e;
t13_l_n=t13_l_n_1+t13_e;
}
else if(x_n>=0){
REAL t13_l_n_1_1=t0_l_n*x;
REAL t13_l_n_1_2=x_n*t0_h;
REAL t13_l_n_1=max(t13_l_n_1_1,t13_l_n_1_2);
REAL t13_h_1_1=t0_l_n*x_n;
REAL t13_h_1_2=t0_h*x;
REAL t13_h_1=max(t13_h_1_1,t13_h_1_2);
t13_h=t13_h_1+t13_e;
t13_l_n=t13_l_n_1+t13_e;
}
else{
REAL t13_l_n_1=t0_l_n*x;
REAL t13_h_1=t0_h*x;
t13_h=t13_h_1+t13_e;
t13_l_n=t13_l_n_1+t13_e;
}
}
else{
if(x<=0){
REAL t13_l_n_1=x_n*t0_h;
REAL x_n=-x;
REAL t13_h_1=x_n*t0_l_n;
t13_h=t13_h_1+t13_e;
t13_l_n=t13_l_n_1+t13_e;
}
else if(x_n>=0){
REAL t13_l_n_1=x_n*t0_h;
REAL t13_h_1=t0_h*x;
t13_h=t13_h_1+t13_e;
t13_l_n=t13_l_n_1+t13_e;
}
else{
REAL t0_l=-t0_l_n;
REAL t13_l_n_1=t0_l*x_n;
REAL t13_h_1=t0_h*x;
t13_h=t13_h_1+t13_e;
t13_l_n=t13_l_n_1+t13_e;
}
}
REAL t14_e_1=fabs(t14);
REAL t14_e=eps*t14_e_1;
REAL t14_h_1=0.166666666666667*t13_h;
REAL t14_h=t14_h_1+t14_e;
REAL t14_l_n_1=0.166666666666667*t13_l_n;
REAL t14_l_n=t14_l_n_1+t14_e;
REAL t15_e_1=fabs(t15);
REAL t15_e=eps*t15_e_1;
REAL t15_h_1=2.48015873015873e-5*t5_h;
REAL t15_h=t15_h_1+t15_e;
REAL t15_l_n_1=2.48015873015873e-5*t5_l_n;
REAL t15_l_n=t15_l_n_1+t15_e;
REAL t16_e_1=fabs(t16);
REAL t16_e=eps*t16_e_1;
REAL  t16_h, t16_l_n;
if(t3_h<=0){
if(x<=0){
REAL t3_h_n=-t3_h;
REAL t16_l_n_1=t3_h_n*x;
REAL t16_h_1=t3_l_n*x_n;
t16_h=t16_h_1+t16_e;
t16_l_n=t16_l_n_1+t16_e;
}
else if(x_n>=0){
REAL t16_l_n_1=t3_l_n*x;
REAL t16_h_1=t3_l_n*x_n;
t16_h=t16_h_1+t16_e;
t16_l_n=t16_l_n_1+t16_e;
}
else{
REAL t16_l_n_1=t3_l_n*x;
REAL t3_h_n=-t3_h;
REAL t16_h_1=t3_h_n*x_n;
t16_h=t16_h_1+t16_e;
t16_l_n=t16_l_n_1+t16_e;
}
}
else if(t3_l_n>=0){
if(x<=0){
REAL t16_l_n_1=t3_h*x_n;
REAL t16_h_1=t3_l_n*x_n;
t16_h=t16_h_1+t16_e;
t16_l_n=t16_l_n_1+t16_e;
}
else if(x_n>=0){
REAL t16_l_n_1_1=t3_l_n*x;
REAL t16_l_n_1_2=x_n*t3_h;
REAL t16_l_n_1=max(t16_l_n_1_1,t16_l_n_1_2);
REAL t16_h_1_1=t3_l_n*x_n;
REAL t16_h_1_2=t3_h*x;
REAL t16_h_1=max(t16_h_1_1,t16_h_1_2);
t16_h=t16_h_1+t16_e;
t16_l_n=t16_l_n_1+t16_e;
}
else{
REAL t16_l_n_1=t3_l_n*x;
REAL t16_h_1=t3_h*x;
t16_h=t16_h_1+t16_e;
t16_l_n=t16_l_n_1+t16_e;
}
}
else{
if(x<=0){
REAL t16_l_n_1=x_n*t3_h;
REAL x_n=-x;
REAL t16_h_1=x_n*t3_l_n;
t16_h=t16_h_1+t16_e;
t16_l_n=t16_l_n_1+t16_e;
}
else if(x_n>=0){
REAL t16_l_n_1=x_n*t3_h;
REAL t16_h_1=t3_h*x;
t16_h=t16_h_1+t16_e;
t16_l_n=t16_l_n_1+t16_e;
}
else{
REAL t3_l=-t3_l_n;
REAL t16_l_n_1=t3_l*x_n;
REAL t16_h_1=t3_h*x;
t16_h=t16_h_1+t16_e;
t16_l_n=t16_l_n_1+t16_e;
}
}
REAL t17_e_1=fabs(t17);
REAL t17_e=eps*t17_e_1;
REAL t17_h_1=0.000198412698412698*t16_h;
REAL t17_h=t17_h_1+t17_e;
REAL t17_l_n_1=0.000198412698412698*t16_l_n;
REAL t17_l_n=t17_l_n_1+t17_e;
REAL t18_e_1=fabs(t18);
REAL t18_e=eps*t18_e_1;
REAL t18_h_1=t1_h+1;
REAL t18_h=t18_h_1+t18_e;
REAL t18_l_n_1=t1_l_n+-1;
REAL t18_l_n=t18_l_n_1+t18_e;
REAL t19_e_1=fabs(t19);
REAL t19_e=eps*t19_e_1;
REAL t19_h_1=t11_h+t18_h;
REAL t19_h=t19_h_1+t19_e;
REAL t19_l_n_1=t11_l_n+t18_l_n;
REAL t19_l_n=t19_l_n_1+t19_e;
REAL t20_e_1=fabs(t20);
REAL t20_e=eps*t20_e_1;
REAL t20_h_1=t12_h+t19_h;
REAL t20_h=t20_h_1+t20_e;
REAL t20_l_n_1=t12_l_n+t19_l_n;
REAL t20_l_n=t20_l_n_1+t20_e;
REAL t21_e_1=fabs(t21);
REAL t21_e=eps*t21_e_1;
REAL t21_h_1=t14_h+t20_h;
REAL t21_h=t21_h_1+t21_e;
REAL t21_l_n_1=t14_l_n+t20_l_n;
REAL t21_l_n=t21_l_n_1+t21_e;
REAL t22_e_1=fabs(t22);
REAL t22_e=eps*t22_e_1;
REAL t22_h_1=t15_h+t21_h;
REAL t22_h=t22_h_1+t22_e;
REAL t22_l_n_1=t15_l_n+t21_l_n;
REAL t22_l_n=t22_l_n_1+t22_e;
REAL t23_e_1=fabs(t23);
REAL t23_e=eps*t23_e_1;
REAL t23_h_1=t17_h+t22_h;
REAL t23_h=t23_h_1+t23_e;
REAL t23_l_n_1=t17_l_n+t22_l_n;
REAL t23_l_n=t23_l_n_1+t23_e;
REAL t24_e_1=fabs(t24);
REAL t24_e=eps*t24_e_1;
REAL t24_h_1=t23_h+t4_h;
REAL t24_h=t24_h_1+t24_e;
REAL t24_l_n_1=t23_l_n+t4_l_n;
REAL t24_l_n=t24_l_n_1+t24_e;
REAL t25_e_1=fabs(t25);
REAL t25_e=eps*t25_e_1;
REAL t25_h_1=t24_h+t7_h;
REAL t25_h=t25_h_1+t25_e;
REAL t25_l_n_1=t24_l_n+t7_l_n;
REAL t25_l_n=t25_l_n_1+t25_e;
REAL t26_e_1=fabs(t26);
REAL t26_e=eps*t26_e_1;
REAL t26_h_1=t25_h+t9_h;
REAL t26_h=t26_h_1+t26_e;
REAL t26_l_n_1=t25_l_n+t9_l_n;
REAL t26_l_n=t26_l_n_1+t26_e;
REAL t27_e_1=fabs(t27);
REAL t27_e=eps*t27_e_1;
REAL t27_h_1=t26_h+x;
REAL t27_h=t27_h_1+t27_e;
REAL t27_l_n_1=t26_l_n+x_n;
REAL t27_l_n=t27_l_n_1+t27_e;
fesetround(roundmode);
*result=t27;
if (inferr==0) {
*hbound=t27_h;
*lbound=-t27_l_n;
}
else{
*hbound=1/0;
*lbound=1/0;
}
}
//Main test Function

int main(){
int range = 126;
int sample=100000;
REAL varlist;
REAL res, err, herr,lerr;
FILE *fsave;
fsave=fopen("ProbSave.txt","w");
FILE *freport;
freport=fopen("ProbReport.txt","w");
int cIA=0;
int cFPS=0;
int cM=0;
int cBoth=0;
int cIAnFPS=0;
int cFPSnIA=0;
int bIA, bFPS;
REAL x;
int i;
for(i=0; i<sample; i++)
{

x=randfl(range);


FPSyn_Polynomial(x,&res,&err);
IA_Polynomial(x,&res,&herr,&lerr);

printf("%e      %e         %e         %f\n ", res, herr-lerr, err, (herr-lerr)/2/err);
fprintf(fsave, "%e      %e      %e      %e      %e         %f\n ", res, err, herr, lerr,herr-lerr, (herr-lerr)/2/err);

if ((0<herr)&&(0>lerr)){bIA=1;}
else { bIA=0;}

if ((res<err)&&(res>-err)){bFPS=1;}
else { bFPS=0;}

cBoth = cBoth+bIA*bFPS;
cIA = cIA + bIA;
cFPS = cFPS+ bFPS;
if ((bIA==1)&&(bFPS==0)) {cIAnFPS=cIAnFPS+1;}
if ((bIA==0)&&(bFPS==1)) {cFPSnIA=cFPSnIA+1;}

}

fclose(fsave);

printf("In error IA: %d\n", cIA);
printf("In error FPSyn: %d\n", cFPS);
printf("In error Both: %d\n", cBoth);

fprintf(freport,"In error IA: %d\n", cIA);
fprintf(freport,"In error FPSyn: %d\n", cFPS);
fprintf(freport,"In error Both: %d\n", cBoth);
fprintf(freport,"In error FPSyn not IA: %d\n", cFPSnIA);
fprintf(freport,"In error IA not FPSyn: %d\n", cIAnFPS);
fclose(freport);

}