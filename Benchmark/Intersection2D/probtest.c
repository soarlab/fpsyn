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
void FPSyn_Intersection2D(REAL x1,REAL x2,REAL x3,REAL x4,REAL y1,REAL y2,REAL y3,REAL y4,REAL* result, REAL* error_bound){
//Calculation steps for input: ;
REAL t0=-x4;
REAL t1=t0+x3;
REAL t2=-x2;
REAL t3=t2+x1;
REAL t4=x1*y2;
REAL t5=t2*y1;
REAL t6=t4+t5;
REAL t7=t1*t6;
REAL t8=x3*y4;
REAL t9=t0*y3;
REAL t10=t8+t9;
REAL t11=t10*t3;
REAL t13=t7-t11;
REAL t15=y3-y4;
REAL t16=t15*t3;
REAL t18=y1-y2;
REAL t19=t1*t18;
REAL t21=t16-t19;
REAL t22=t13/t21;
//Calculation steps for Error-bound: ;
t1=fabs(t1);
t9=fabs(t9);
t11=fabs(t11);
t8=fabs(t8);
t3=fabs(t3);
t4=fabs(t4);
t5=fabs(t5);
t7=fabs(t7);
t16=fabs(t16);
t19=fabs(t19);
t13=fabs(t13);
t21=fabs(t21);
REAL temp0=t4+t5;
REAL temp1=t1*temp0;
REAL temp2=t8+t9;
REAL temp3=t3*temp2;
REAL temp4=temp1+temp3;
REAL temp5=5.55111512312578e-16*temp4;
REAL temp6=t16+t19;
REAL temp7=4.44089209850063e-16*temp6;
*result=t22;
t22=fabs(t22);
REAL A=(t13*temp7+t21*temp5)*1.0000000001;
REAL B=(t21-1.0000000001*temp7)*0.9999999999999998;
REAL E=(A/t21*1.0000000000000002/B*1.0000000000000002+t22*1.1102230246251565e-16*1.0000000000000002)*1.0000000000000002;
*error_bound=t21<=temp7?1/0:E;
}
//IA Function
void IA_Intersection2D(REAL x1,REAL x2,REAL x3,REAL x4,REAL y1,REAL y2,REAL y3,REAL y4,REAL* result, REAL* hbound, REAL* lbound){
REAL t0=-x4;
REAL t1=t0+x3;
REAL t2=-x2;
REAL t3=t2+x1;
REAL t4=x1*y2;
REAL t5=t2*y1;
REAL t6=t4+t5;
REAL t7=t1*t6;
REAL t8=x3*y4;
REAL t9=t0*y3;
REAL t10=t8+t9;
REAL t11=t10*t3;
REAL t13=t7-t11;
REAL t15=y3-y4;
REAL t16=t15*t3;
REAL t18=y1-y2;
REAL t19=t1*t18;
REAL t21=t16-t19;
REAL t22=t13/t21;
//Calculate bound
REAL eps=pow(2,-53);
int inferr=0;
int roundmode = fegetround();
fesetround(FE_UPWARD);
REAL x4_n=-x4;
REAL t0_h=x4_n;
REAL t0_l_n=x4;
REAL t1_e_1=fabs(t1);
REAL t1_e=eps*t1_e_1;
REAL x3_n=-x3;
REAL t1_h_1=t0_h+x3;
REAL t1_h=t1_h_1+t1_e;
REAL t1_l_n_1=t0_l_n+x3_n;
REAL t1_l_n=t1_l_n_1+t1_e;
REAL x2_n=-x2;
REAL t2_h=x2_n;
REAL t2_l_n=x2;
REAL t3_e_1=fabs(t3);
REAL t3_e=eps*t3_e_1;
REAL x1_n=-x1;
REAL t3_h_1=t2_h+x1;
REAL t3_h=t3_h_1+t3_e;
REAL t3_l_n_1=t2_l_n+x1_n;
REAL t3_l_n=t3_l_n_1+t3_e;
REAL t4_e_1=fabs(t4);
REAL t4_e=eps*t4_e_1;
REAL t4_h=t4+t4_e;
REAL t4_l_n=t4_e-t4;
REAL t5_e_1=fabs(t5);
REAL t5_e=eps*t5_e_1;
REAL y1_n=-y1;
REAL  t5_h, t5_l_n;
if(t2_h<=0){
if(y1<=0){
REAL t2_h_n=-t2_h;
REAL t5_l_n_1=t2_h_n*y1;
REAL t5_h_1=t2_l_n*y1_n;
t5_h=t5_h_1+t5_e;
t5_l_n=t5_l_n_1+t5_e;
}
else if(y1_n>=0){
REAL t5_l_n_1=t2_l_n*y1;
REAL t5_h_1=t2_l_n*y1_n;
t5_h=t5_h_1+t5_e;
t5_l_n=t5_l_n_1+t5_e;
}
else{
REAL t5_l_n_1=t2_l_n*y1;
REAL t2_h_n=-t2_h;
REAL t5_h_1=t2_h_n*y1_n;
t5_h=t5_h_1+t5_e;
t5_l_n=t5_l_n_1+t5_e;
}
}
else if(t2_l_n>=0){
if(y1<=0){
REAL t5_l_n_1=t2_h*y1_n;
REAL t5_h_1=t2_l_n*y1_n;
t5_h=t5_h_1+t5_e;
t5_l_n=t5_l_n_1+t5_e;
}
else if(y1_n>=0){
REAL t5_l_n_1_1=t2_l_n*y1;
REAL t5_l_n_1_2=y1_n*t2_h;
REAL t5_l_n_1=max(t5_l_n_1_1,t5_l_n_1_2);
REAL t5_h_1_1=t2_l_n*y1_n;
REAL t5_h_1_2=t2_h*y1;
REAL t5_h_1=max(t5_h_1_1,t5_h_1_2);
t5_h=t5_h_1+t5_e;
t5_l_n=t5_l_n_1+t5_e;
}
else{
REAL t5_l_n_1=t2_l_n*y1;
REAL t5_h_1=t2_h*y1;
t5_h=t5_h_1+t5_e;
t5_l_n=t5_l_n_1+t5_e;
}
}
else{
if(y1<=0){
REAL t5_l_n_1=y1_n*t2_h;
REAL y1_n=-y1;
REAL t5_h_1=y1_n*t2_l_n;
t5_h=t5_h_1+t5_e;
t5_l_n=t5_l_n_1+t5_e;
}
else if(y1_n>=0){
REAL t5_l_n_1=y1_n*t2_h;
REAL t5_h_1=t2_h*y1;
t5_h=t5_h_1+t5_e;
t5_l_n=t5_l_n_1+t5_e;
}
else{
REAL t2_l=-t2_l_n;
REAL t5_l_n_1=t2_l*y1_n;
REAL t5_h_1=t2_h*y1;
t5_h=t5_h_1+t5_e;
t5_l_n=t5_l_n_1+t5_e;
}
}
REAL t6_e_1=fabs(t6);
REAL t6_e=eps*t6_e_1;
REAL t6_h_1=t4_h+t5_h;
REAL t6_h=t6_h_1+t6_e;
REAL t6_l_n_1=t4_l_n+t5_l_n;
REAL t6_l_n=t6_l_n_1+t6_e;
REAL t7_e_1=fabs(t7);
REAL t7_e=eps*t7_e_1;
REAL  t7_h, t7_l_n;
if(t1_h<=0){
if(t6_h<=0){
REAL t1_h_n=-t1_h;
REAL t7_l_n_1=t1_h_n*t6_h;
REAL t7_h_1=t1_l_n*t6_l_n;
t7_h=t7_h_1+t7_e;
t7_l_n=t7_l_n_1+t7_e;
}
else if(t6_l_n>=0){
REAL t7_l_n_1=t1_l_n*t6_h;
REAL t7_h_1=t1_l_n*t6_l_n;
t7_h=t7_h_1+t7_e;
t7_l_n=t7_l_n_1+t7_e;
}
else{
REAL t7_l_n_1=t1_l_n*t6_h;
REAL t1_h_n=-t1_h;
REAL t7_h_1=t1_h_n*t6_l_n;
t7_h=t7_h_1+t7_e;
t7_l_n=t7_l_n_1+t7_e;
}
}
else if(t1_l_n>=0){
if(t6_h<=0){
REAL t7_l_n_1=t1_h*t6_l_n;
REAL t7_h_1=t1_l_n*t6_l_n;
t7_h=t7_h_1+t7_e;
t7_l_n=t7_l_n_1+t7_e;
}
else if(t6_l_n>=0){
REAL t7_l_n_1_1=t1_l_n*t6_h;
REAL t7_l_n_1_2=t6_l_n*t1_h;
REAL t7_l_n_1=max(t7_l_n_1_1,t7_l_n_1_2);
REAL t7_h_1_1=t1_l_n*t6_l_n;
REAL t7_h_1_2=t1_h*t6_h;
REAL t7_h_1=max(t7_h_1_1,t7_h_1_2);
t7_h=t7_h_1+t7_e;
t7_l_n=t7_l_n_1+t7_e;
}
else{
REAL t7_l_n_1=t1_l_n*t6_h;
REAL t7_h_1=t1_h*t6_h;
t7_h=t7_h_1+t7_e;
t7_l_n=t7_l_n_1+t7_e;
}
}
else{
if(t6_h<=0){
REAL t7_l_n_1=t6_l_n*t1_h;
REAL t6_h_n=-t6_h;
REAL t7_h_1=t6_h_n*t1_l_n;
t7_h=t7_h_1+t7_e;
t7_l_n=t7_l_n_1+t7_e;
}
else if(t6_l_n>=0){
REAL t7_l_n_1=t6_l_n*t1_h;
REAL t7_h_1=t1_h*t6_h;
t7_h=t7_h_1+t7_e;
t7_l_n=t7_l_n_1+t7_e;
}
else{
REAL t1_l=-t1_l_n;
REAL t7_l_n_1=t1_l*t6_l_n;
REAL t7_h_1=t1_h*t6_h;
t7_h=t7_h_1+t7_e;
t7_l_n=t7_l_n_1+t7_e;
}
}
REAL t8_e_1=fabs(t8);
REAL t8_e=eps*t8_e_1;
REAL t8_h=t8+t8_e;
REAL t8_l_n=t8_e-t8;
REAL t9_e_1=fabs(t9);
REAL t9_e=eps*t9_e_1;
REAL y3_n=-y3;
REAL  t9_h, t9_l_n;
if(t0_h<=0){
if(y3<=0){
REAL t0_h_n=-t0_h;
REAL t9_l_n_1=t0_h_n*y3;
REAL t9_h_1=t0_l_n*y3_n;
t9_h=t9_h_1+t9_e;
t9_l_n=t9_l_n_1+t9_e;
}
else if(y3_n>=0){
REAL t9_l_n_1=t0_l_n*y3;
REAL t9_h_1=t0_l_n*y3_n;
t9_h=t9_h_1+t9_e;
t9_l_n=t9_l_n_1+t9_e;
}
else{
REAL t9_l_n_1=t0_l_n*y3;
REAL t0_h_n=-t0_h;
REAL t9_h_1=t0_h_n*y3_n;
t9_h=t9_h_1+t9_e;
t9_l_n=t9_l_n_1+t9_e;
}
}
else if(t0_l_n>=0){
if(y3<=0){
REAL t9_l_n_1=t0_h*y3_n;
REAL t9_h_1=t0_l_n*y3_n;
t9_h=t9_h_1+t9_e;
t9_l_n=t9_l_n_1+t9_e;
}
else if(y3_n>=0){
REAL t9_l_n_1_1=t0_l_n*y3;
REAL t9_l_n_1_2=y3_n*t0_h;
REAL t9_l_n_1=max(t9_l_n_1_1,t9_l_n_1_2);
REAL t9_h_1_1=t0_l_n*y3_n;
REAL t9_h_1_2=t0_h*y3;
REAL t9_h_1=max(t9_h_1_1,t9_h_1_2);
t9_h=t9_h_1+t9_e;
t9_l_n=t9_l_n_1+t9_e;
}
else{
REAL t9_l_n_1=t0_l_n*y3;
REAL t9_h_1=t0_h*y3;
t9_h=t9_h_1+t9_e;
t9_l_n=t9_l_n_1+t9_e;
}
}
else{
if(y3<=0){
REAL t9_l_n_1=y3_n*t0_h;
REAL y3_n=-y3;
REAL t9_h_1=y3_n*t0_l_n;
t9_h=t9_h_1+t9_e;
t9_l_n=t9_l_n_1+t9_e;
}
else if(y3_n>=0){
REAL t9_l_n_1=y3_n*t0_h;
REAL t9_h_1=t0_h*y3;
t9_h=t9_h_1+t9_e;
t9_l_n=t9_l_n_1+t9_e;
}
else{
REAL t0_l=-t0_l_n;
REAL t9_l_n_1=t0_l*y3_n;
REAL t9_h_1=t0_h*y3;
t9_h=t9_h_1+t9_e;
t9_l_n=t9_l_n_1+t9_e;
}
}
REAL t10_e_1=fabs(t10);
REAL t10_e=eps*t10_e_1;
REAL t10_h_1=t8_h+t9_h;
REAL t10_h=t10_h_1+t10_e;
REAL t10_l_n_1=t8_l_n+t9_l_n;
REAL t10_l_n=t10_l_n_1+t10_e;
REAL t11_e_1=fabs(t11);
REAL t11_e=eps*t11_e_1;
REAL  t11_h, t11_l_n;
if(t10_h<=0){
if(t3_h<=0){
REAL t10_h_n=-t10_h;
REAL t11_l_n_1=t10_h_n*t3_h;
REAL t11_h_1=t10_l_n*t3_l_n;
t11_h=t11_h_1+t11_e;
t11_l_n=t11_l_n_1+t11_e;
}
else if(t3_l_n>=0){
REAL t11_l_n_1=t10_l_n*t3_h;
REAL t11_h_1=t10_l_n*t3_l_n;
t11_h=t11_h_1+t11_e;
t11_l_n=t11_l_n_1+t11_e;
}
else{
REAL t11_l_n_1=t10_l_n*t3_h;
REAL t10_h_n=-t10_h;
REAL t11_h_1=t10_h_n*t3_l_n;
t11_h=t11_h_1+t11_e;
t11_l_n=t11_l_n_1+t11_e;
}
}
else if(t10_l_n>=0){
if(t3_h<=0){
REAL t11_l_n_1=t10_h*t3_l_n;
REAL t11_h_1=t10_l_n*t3_l_n;
t11_h=t11_h_1+t11_e;
t11_l_n=t11_l_n_1+t11_e;
}
else if(t3_l_n>=0){
REAL t11_l_n_1_1=t10_l_n*t3_h;
REAL t11_l_n_1_2=t3_l_n*t10_h;
REAL t11_l_n_1=max(t11_l_n_1_1,t11_l_n_1_2);
REAL t11_h_1_1=t10_l_n*t3_l_n;
REAL t11_h_1_2=t10_h*t3_h;
REAL t11_h_1=max(t11_h_1_1,t11_h_1_2);
t11_h=t11_h_1+t11_e;
t11_l_n=t11_l_n_1+t11_e;
}
else{
REAL t11_l_n_1=t10_l_n*t3_h;
REAL t11_h_1=t10_h*t3_h;
t11_h=t11_h_1+t11_e;
t11_l_n=t11_l_n_1+t11_e;
}
}
else{
if(t3_h<=0){
REAL t11_l_n_1=t3_l_n*t10_h;
REAL t3_h_n=-t3_h;
REAL t11_h_1=t3_h_n*t10_l_n;
t11_h=t11_h_1+t11_e;
t11_l_n=t11_l_n_1+t11_e;
}
else if(t3_l_n>=0){
REAL t11_l_n_1=t3_l_n*t10_h;
REAL t11_h_1=t10_h*t3_h;
t11_h=t11_h_1+t11_e;
t11_l_n=t11_l_n_1+t11_e;
}
else{
REAL t10_l=-t10_l_n;
REAL t11_l_n_1=t10_l*t3_l_n;
REAL t11_h_1=t10_h*t3_h;
t11_h=t11_h_1+t11_e;
t11_l_n=t11_l_n_1+t11_e;
}
}
REAL t13_e_1=fabs(t13);
REAL t13_e=eps*t13_e_1;
REAL t13_h_1=t7_h+t11_l_n;
REAL t13_h=t13_h_1+t13_e;
REAL t13_l_n_1=t11_h+t7_l_n;
REAL t13_l_n=t13_l_n_1+t13_e;
REAL t15_e_1=fabs(t15);
REAL t15_e=eps*t15_e_1;
REAL t15_h=t15+t15_e;
REAL t15_l_n=t15_e-t15;
REAL t16_e_1=fabs(t16);
REAL t16_e=eps*t16_e_1;
REAL  t16_h, t16_l_n;
if(t15_h<=0){
if(t3_h<=0){
REAL t15_h_n=-t15_h;
REAL t16_l_n_1=t15_h_n*t3_h;
REAL t16_h_1=t15_l_n*t3_l_n;
t16_h=t16_h_1+t16_e;
t16_l_n=t16_l_n_1+t16_e;
}
else if(t3_l_n>=0){
REAL t16_l_n_1=t15_l_n*t3_h;
REAL t16_h_1=t15_l_n*t3_l_n;
t16_h=t16_h_1+t16_e;
t16_l_n=t16_l_n_1+t16_e;
}
else{
REAL t16_l_n_1=t15_l_n*t3_h;
REAL t15_h_n=-t15_h;
REAL t16_h_1=t15_h_n*t3_l_n;
t16_h=t16_h_1+t16_e;
t16_l_n=t16_l_n_1+t16_e;
}
}
else if(t15_l_n>=0){
if(t3_h<=0){
REAL t16_l_n_1=t15_h*t3_l_n;
REAL t16_h_1=t15_l_n*t3_l_n;
t16_h=t16_h_1+t16_e;
t16_l_n=t16_l_n_1+t16_e;
}
else if(t3_l_n>=0){
REAL t16_l_n_1_1=t15_l_n*t3_h;
REAL t16_l_n_1_2=t3_l_n*t15_h;
REAL t16_l_n_1=max(t16_l_n_1_1,t16_l_n_1_2);
REAL t16_h_1_1=t15_l_n*t3_l_n;
REAL t16_h_1_2=t15_h*t3_h;
REAL t16_h_1=max(t16_h_1_1,t16_h_1_2);
t16_h=t16_h_1+t16_e;
t16_l_n=t16_l_n_1+t16_e;
}
else{
REAL t16_l_n_1=t15_l_n*t3_h;
REAL t16_h_1=t15_h*t3_h;
t16_h=t16_h_1+t16_e;
t16_l_n=t16_l_n_1+t16_e;
}
}
else{
if(t3_h<=0){
REAL t16_l_n_1=t3_l_n*t15_h;
REAL t3_h_n=-t3_h;
REAL t16_h_1=t3_h_n*t15_l_n;
t16_h=t16_h_1+t16_e;
t16_l_n=t16_l_n_1+t16_e;
}
else if(t3_l_n>=0){
REAL t16_l_n_1=t3_l_n*t15_h;
REAL t16_h_1=t15_h*t3_h;
t16_h=t16_h_1+t16_e;
t16_l_n=t16_l_n_1+t16_e;
}
else{
REAL t15_l=-t15_l_n;
REAL t16_l_n_1=t15_l*t3_l_n;
REAL t16_h_1=t15_h*t3_h;
t16_h=t16_h_1+t16_e;
t16_l_n=t16_l_n_1+t16_e;
}
}
REAL t18_e_1=fabs(t18);
REAL t18_e=eps*t18_e_1;
REAL t18_h=t18+t18_e;
REAL t18_l_n=t18_e-t18;
REAL t19_e_1=fabs(t19);
REAL t19_e=eps*t19_e_1;
REAL  t19_h, t19_l_n;
if(t1_h<=0){
if(t18_h<=0){
REAL t1_h_n=-t1_h;
REAL t19_l_n_1=t1_h_n*t18_h;
REAL t19_h_1=t1_l_n*t18_l_n;
t19_h=t19_h_1+t19_e;
t19_l_n=t19_l_n_1+t19_e;
}
else if(t18_l_n>=0){
REAL t19_l_n_1=t1_l_n*t18_h;
REAL t19_h_1=t1_l_n*t18_l_n;
t19_h=t19_h_1+t19_e;
t19_l_n=t19_l_n_1+t19_e;
}
else{
REAL t19_l_n_1=t1_l_n*t18_h;
REAL t1_h_n=-t1_h;
REAL t19_h_1=t1_h_n*t18_l_n;
t19_h=t19_h_1+t19_e;
t19_l_n=t19_l_n_1+t19_e;
}
}
else if(t1_l_n>=0){
if(t18_h<=0){
REAL t19_l_n_1=t1_h*t18_l_n;
REAL t19_h_1=t1_l_n*t18_l_n;
t19_h=t19_h_1+t19_e;
t19_l_n=t19_l_n_1+t19_e;
}
else if(t18_l_n>=0){
REAL t19_l_n_1_1=t1_l_n*t18_h;
REAL t19_l_n_1_2=t18_l_n*t1_h;
REAL t19_l_n_1=max(t19_l_n_1_1,t19_l_n_1_2);
REAL t19_h_1_1=t1_l_n*t18_l_n;
REAL t19_h_1_2=t1_h*t18_h;
REAL t19_h_1=max(t19_h_1_1,t19_h_1_2);
t19_h=t19_h_1+t19_e;
t19_l_n=t19_l_n_1+t19_e;
}
else{
REAL t19_l_n_1=t1_l_n*t18_h;
REAL t19_h_1=t1_h*t18_h;
t19_h=t19_h_1+t19_e;
t19_l_n=t19_l_n_1+t19_e;
}
}
else{
if(t18_h<=0){
REAL t19_l_n_1=t18_l_n*t1_h;
REAL t18_h_n=-t18_h;
REAL t19_h_1=t18_h_n*t1_l_n;
t19_h=t19_h_1+t19_e;
t19_l_n=t19_l_n_1+t19_e;
}
else if(t18_l_n>=0){
REAL t19_l_n_1=t18_l_n*t1_h;
REAL t19_h_1=t1_h*t18_h;
t19_h=t19_h_1+t19_e;
t19_l_n=t19_l_n_1+t19_e;
}
else{
REAL t1_l=-t1_l_n;
REAL t19_l_n_1=t1_l*t18_l_n;
REAL t19_h_1=t1_h*t18_h;
t19_h=t19_h_1+t19_e;
t19_l_n=t19_l_n_1+t19_e;
}
}
REAL t21_e_1=fabs(t21);
REAL t21_e=eps*t21_e_1;
REAL t21_h_1=t16_h+t19_l_n;
REAL t21_h=t21_h_1+t21_e;
REAL t21_l_n_1=t19_h+t16_l_n;
REAL t21_l_n=t21_l_n_1+t21_e;
REAL t22_e_1=fabs(t22);
REAL t22_e=eps*t22_e_1;
if(t21_l_n*t21_h>=0) inferr=1;
REAL t21_l_n_inv=1/t21_l_n;
REAL t21_h_inv=1/t21_h;
REAL  t22_h, t22_l_n;
if(t13_h<=0){
if(t21_l_n_inv<=0){
REAL t13_h_n=-t13_h;
REAL t22_l_n_1=t13_h_n*t21_l_n_inv;
REAL t22_h_1=t13_l_n*t21_h_inv;
t22_h=t22_h_1+t22_e;
t22_l_n=t22_l_n_1+t22_e;
}
else if(t21_h_inv>=0){
REAL t22_l_n_1=t13_l_n*t21_l_n_inv;
REAL t22_h_1=t13_l_n*t21_h_inv;
t22_h=t22_h_1+t22_e;
t22_l_n=t22_l_n_1+t22_e;
}
else{
REAL t22_l_n_1=t13_l_n*t21_l_n_inv;
REAL t13_h_n=-t13_h;
REAL t22_h_1=t13_h_n*t21_h_inv;
t22_h=t22_h_1+t22_e;
t22_l_n=t22_l_n_1+t22_e;
}
}
else if(t13_l_n>=0){
if(t21_l_n_inv<=0){
REAL t22_l_n_1=t13_h*t21_h_inv;
REAL t22_h_1=t13_l_n*t21_h_inv;
t22_h=t22_h_1+t22_e;
t22_l_n=t22_l_n_1+t22_e;
}
else if(t21_h_inv>=0){
REAL t22_l_n_1_1=t13_l_n*t21_l_n_inv;
REAL t22_l_n_1_2=t21_h_inv*t13_h;
REAL t22_l_n_1=max(t22_l_n_1_1,t22_l_n_1_2);
REAL t22_h_1_1=t13_l_n*t21_h_inv;
REAL t22_h_1_2=t13_h*t21_l_n_inv;
REAL t22_h_1=max(t22_h_1_1,t22_h_1_2);
t22_h=t22_h_1+t22_e;
t22_l_n=t22_l_n_1+t22_e;
}
else{
REAL t22_l_n_1=t13_l_n*t21_l_n_inv;
REAL t22_h_1=t13_h*t21_l_n_inv;
t22_h=t22_h_1+t22_e;
t22_l_n=t22_l_n_1+t22_e;
}
}
else{
if(t21_l_n_inv<=0){
REAL t22_l_n_1=t21_h_inv*t13_h;
REAL t21_l_n_inv_n=-t21_l_n_inv;
REAL t22_h_1=t21_l_n_inv_n*t13_l_n;
t22_h=t22_h_1+t22_e;
t22_l_n=t22_l_n_1+t22_e;
}
else if(t21_h_inv>=0){
REAL t22_l_n_1=t21_h_inv*t13_h;
REAL t22_h_1=t13_h*t21_l_n_inv;
t22_h=t22_h_1+t22_e;
t22_l_n=t22_l_n_1+t22_e;
}
else{
REAL t13_l=-t13_l_n;
REAL t22_l_n_1=t13_l*t21_h_inv;
REAL t22_h_1=t13_h*t21_l_n_inv;
t22_h=t22_h_1+t22_e;
t22_l_n=t22_l_n_1+t22_e;
}
}
fesetround(roundmode);
*result=t22;
if (inferr==0) {
*hbound=t22_h;
*lbound=-t22_l_n;
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
REAL x1,x2,x3,x4,y1,y2,y3,y4;
int i;
for(i=0; i<sample; i++)
{

x1=randfl(range);
x2=randfl(range);
x3=randfl(range);
x4=randfl(range);
y1=randfl(range);
y2=randfl(range);
y3=randfl(range);
y4=randfl(range);


FPSyn_Intersection2D(x1,x2,x3,x4,y1,y2,y3,y4,&res,&err);
IA_Intersection2D(x1,x2,x3,x4,y1,y2,y3,y4,&res,&herr,&lerr);

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