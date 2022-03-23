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
void FPSyn_ConvexHullArea(REAL x0,REAL x1,REAL x2,REAL x3,REAL x4,REAL x5,REAL x6,REAL x7,REAL x8,REAL x9,REAL y0,REAL y1,REAL y2,REAL y3,REAL y4,REAL y5,REAL y6,REAL y7,REAL y8,REAL y9,REAL* result, REAL* error_bound){
REAL t2 = y1-y9;
REAL t3 = x1-x9;
REAL t4 = y2-y9;
REAL t5 = x2-x9;
REAL t6 = y3-y9;
REAL t7 = x3-x9;
REAL t8 = y4-y9;
REAL t9 = x4-x9;
REAL t10 = y5-y9;
REAL t11 = x5-x9;
REAL t12 = y6-y9;
REAL t13 = x6-x9;
REAL t14 = y7-y9;
REAL t15 = x7-x9;
REAL t17 = 0.5*(t10*t9);
REAL t19 = 0.5*(t11*t12);
REAL t21 = 0.5*(t13*t14);
REAL t24 = 0.5*(t15*(y8-y9));
REAL t27 = 0.5*(t2*(x0-x9));
REAL t29 = 0.5*(t3*t4);
REAL t31 = 0.5*(t5*t6);
REAL t33 = 0.5*(t7*t8);
REAL t35 = -0.5*(t10*t13);
REAL t37 = -0.5*(t11*t8);
REAL t39 = -0.5*(t12*t15);
REAL t42 = -0.5*(t14*(x8-x9));
REAL t44 = -0.5*(t2*t5);
REAL t47 = -0.5*(t3*(y0-y9));
REAL t49 = -0.5*(t4*t7);
REAL t51 = -0.5*(t6*t9);
REAL t53 = t21+(t17+t19);
REAL t55 = t27+(t24+t53);
REAL t57 = t31+(t29+t55);
REAL t59 = t35+(t33+t57);
REAL t61 = t39+(t37+t59);
REAL t63 = t44+(t42+t61);
REAL t65 = t49+(t47+t63);
*result = (t51+t65);
*error_bound = 1.0000000001*((5.55111512312578e-16*((Absolute(t63))+((Absolute(t61))+((Absolute(t55))+((Absolute(t49))+((Absolute(t47))+((Absolute(t37))+((Absolute(t35))+((Absolute(t33))+((Absolute(t29))+((Absolute(t27))+((Absolute(t24))+((Absolute(t19))+(Absolute(t21)))))))))))))))+(5.55111512312578e-16*((Absolute(t65))+((Absolute(t59))+((Absolute(t57))+((Absolute(t53))+((Absolute(t51))+((Absolute(t44))+((Absolute(t42))+((Absolute(t39))+((Absolute(t17))+(Absolute(t31)))))))))))));
}//IA Function
void IA_ConvexHullArea(REAL x0,REAL x1,REAL x2,REAL x3,REAL x4,REAL x5,REAL x6,REAL x7,REAL x8,REAL x9,REAL y0,REAL y1,REAL y2,REAL y3,REAL y4,REAL y5,REAL y6,REAL y7,REAL y8,REAL y9,REAL* result, REAL* hbound, REAL* lbound){
REAL t2=y1-y9;
REAL t3=x1-x9;
REAL t4=y2-y9;
REAL t5=x2-x9;
REAL t6=y3-y9;
REAL t7=x3-x9;
REAL t8=y4-y9;
REAL t9=x4-x9;
REAL t10=y5-y9;
REAL t11=x5-x9;
REAL t12=y6-y9;
REAL t13=x6-x9;
REAL t14=y7-y9;
REAL t15=x7-x9;
REAL t16=t10*t9;
REAL t17=0.5*t16;
REAL t18=t11*t12;
REAL t19=0.5*t18;
REAL t20=t13*t14;
REAL t21=0.5*t20;
REAL t22=y8-y9;
REAL t23=t15*t22;
REAL t24=0.5*t23;
REAL t25=x0-x9;
REAL t26=t2*t25;
REAL t27=0.5*t26;
REAL t28=t3*t4;
REAL t29=0.5*t28;
REAL t30=t5*t6;
REAL t31=0.5*t30;
REAL t32=t7*t8;
REAL t33=0.5*t32;
REAL t34=t10*t13;
REAL t35=-0.5*t34;
REAL t36=t11*t8;
REAL t37=-0.5*t36;
REAL t38=t12*t15;
REAL t39=-0.5*t38;
REAL t40=x8-x9;
REAL t41=t14*t40;
REAL t42=-0.5*t41;
REAL t43=t2*t5;
REAL t44=-0.5*t43;
REAL t45=y0-y9;
REAL t46=t3*t45;
REAL t47=-0.5*t46;
REAL t48=t4*t7;
REAL t49=-0.5*t48;
REAL t50=t6*t9;
REAL t51=-0.5*t50;
REAL t52=t17+t19;
REAL t53=t21+t52;
REAL t54=t24+t53;
REAL t55=t27+t54;
REAL t56=t29+t55;
REAL t57=t31+t56;
REAL t58=t33+t57;
REAL t59=t35+t58;
REAL t60=t37+t59;
REAL t61=t39+t60;
REAL t62=t42+t61;
REAL t63=t44+t62;
REAL t64=t47+t63;
REAL t65=t49+t64;
REAL t66=t51+t65;
//Calculate bound
REAL eps=pow(2,-53);
int inferr=0;
int roundmode = fegetround();
fesetround(FE_UPWARD);
REAL t2_e_1=fabs(t2);
REAL t2_e=eps*t2_e_1;
REAL t2_h=t2+t2_e;
REAL t2_l_n=t2_e-t2;
REAL t3_e_1=fabs(t3);
REAL t3_e=eps*t3_e_1;
REAL t3_h=t3+t3_e;
REAL t3_l_n=t3_e-t3;
REAL t4_e_1=fabs(t4);
REAL t4_e=eps*t4_e_1;
REAL t4_h=t4+t4_e;
REAL t4_l_n=t4_e-t4;
REAL t5_e_1=fabs(t5);
REAL t5_e=eps*t5_e_1;
REAL t5_h=t5+t5_e;
REAL t5_l_n=t5_e-t5;
REAL t6_e_1=fabs(t6);
REAL t6_e=eps*t6_e_1;
REAL t6_h=t6+t6_e;
REAL t6_l_n=t6_e-t6;
REAL t7_e_1=fabs(t7);
REAL t7_e=eps*t7_e_1;
REAL t7_h=t7+t7_e;
REAL t7_l_n=t7_e-t7;
REAL t8_e_1=fabs(t8);
REAL t8_e=eps*t8_e_1;
REAL t8_h=t8+t8_e;
REAL t8_l_n=t8_e-t8;
REAL t9_e_1=fabs(t9);
REAL t9_e=eps*t9_e_1;
REAL t9_h=t9+t9_e;
REAL t9_l_n=t9_e-t9;
REAL t10_e_1=fabs(t10);
REAL t10_e=eps*t10_e_1;
REAL t10_h=t10+t10_e;
REAL t10_l_n=t10_e-t10;
REAL t11_e_1=fabs(t11);
REAL t11_e=eps*t11_e_1;
REAL t11_h=t11+t11_e;
REAL t11_l_n=t11_e-t11;
REAL t12_e_1=fabs(t12);
REAL t12_e=eps*t12_e_1;
REAL t12_h=t12+t12_e;
REAL t12_l_n=t12_e-t12;
REAL t13_e_1=fabs(t13);
REAL t13_e=eps*t13_e_1;
REAL t13_h=t13+t13_e;
REAL t13_l_n=t13_e-t13;
REAL t14_e_1=fabs(t14);
REAL t14_e=eps*t14_e_1;
REAL t14_h=t14+t14_e;
REAL t14_l_n=t14_e-t14;
REAL t15_e_1=fabs(t15);
REAL t15_e=eps*t15_e_1;
REAL t15_h=t15+t15_e;
REAL t15_l_n=t15_e-t15;
REAL t16_e_1=fabs(t16);
REAL t16_e=eps*t16_e_1;
REAL  t16_h, t16_l_n;
if(t10_h<=0){
if(t9_h<=0){
REAL t10_h_n=-t10_h;
REAL t16_l_n_1=t10_h_n*t9_h;
REAL t16_h_1=t10_l_n*t9_l_n;
t16_h=t16_h_1+t16_e;
t16_l_n=t16_l_n_1+t16_e;
}
else if(t9_l_n>=0){
REAL t16_l_n_1=t10_l_n*t9_h;
REAL t16_h_1=t10_l_n*t9_l_n;
t16_h=t16_h_1+t16_e;
t16_l_n=t16_l_n_1+t16_e;
}
else{
REAL t16_l_n_1=t10_l_n*t9_h;
REAL t10_h_n=-t10_h;
REAL t16_h_1=t10_h_n*t9_l_n;
t16_h=t16_h_1+t16_e;
t16_l_n=t16_l_n_1+t16_e;
}
}
else if(t10_l_n>=0){
if(t9_h<=0){
REAL t16_l_n_1=t10_h*t9_l_n;
REAL t16_h_1=t10_l_n*t9_l_n;
t16_h=t16_h_1+t16_e;
t16_l_n=t16_l_n_1+t16_e;
}
else if(t9_l_n>=0){
REAL t16_l_n_1_1=t10_l_n*t9_h;
REAL t16_l_n_1_2=t9_l_n*t10_h;
REAL t16_l_n_1=max(t16_l_n_1_1,t16_l_n_1_2);
REAL t16_h_1_1=t10_l_n*t9_l_n;
REAL t16_h_1_2=t10_h*t9_h;
REAL t16_h_1=max(t16_h_1_1,t16_h_1_2);
t16_h=t16_h_1+t16_e;
t16_l_n=t16_l_n_1+t16_e;
}
else{
REAL t16_l_n_1=t10_l_n*t9_h;
REAL t16_h_1=t10_h*t9_h;
t16_h=t16_h_1+t16_e;
t16_l_n=t16_l_n_1+t16_e;
}
}
else{
if(t9_h<=0){
REAL t16_l_n_1=t9_l_n*t10_h;
REAL t9_h_n=-t9_h;
REAL t16_h_1=t9_h_n*t10_l_n;
t16_h=t16_h_1+t16_e;
t16_l_n=t16_l_n_1+t16_e;
}
else if(t9_l_n>=0){
REAL t16_l_n_1=t9_l_n*t10_h;
REAL t16_h_1=t10_h*t9_h;
t16_h=t16_h_1+t16_e;
t16_l_n=t16_l_n_1+t16_e;
}
else{
REAL t10_l=-t10_l_n;
REAL t16_l_n_1=t10_l*t9_l_n;
REAL t16_h_1=t10_h*t9_h;
t16_h=t16_h_1+t16_e;
t16_l_n=t16_l_n_1+t16_e;
}
}
REAL t17_e_1=fabs(t17);
REAL t17_e=eps*t17_e_1;
REAL t17_h_1=0.5*t16_h;
REAL t17_h=t17_h_1+t17_e;
REAL t17_l_n_1=0.5*t16_l_n;
REAL t17_l_n=t17_l_n_1+t17_e;
REAL t18_e_1=fabs(t18);
REAL t18_e=eps*t18_e_1;
REAL  t18_h, t18_l_n;
if(t11_h<=0){
if(t12_h<=0){
REAL t11_h_n=-t11_h;
REAL t18_l_n_1=t11_h_n*t12_h;
REAL t18_h_1=t11_l_n*t12_l_n;
t18_h=t18_h_1+t18_e;
t18_l_n=t18_l_n_1+t18_e;
}
else if(t12_l_n>=0){
REAL t18_l_n_1=t11_l_n*t12_h;
REAL t18_h_1=t11_l_n*t12_l_n;
t18_h=t18_h_1+t18_e;
t18_l_n=t18_l_n_1+t18_e;
}
else{
REAL t18_l_n_1=t11_l_n*t12_h;
REAL t11_h_n=-t11_h;
REAL t18_h_1=t11_h_n*t12_l_n;
t18_h=t18_h_1+t18_e;
t18_l_n=t18_l_n_1+t18_e;
}
}
else if(t11_l_n>=0){
if(t12_h<=0){
REAL t18_l_n_1=t11_h*t12_l_n;
REAL t18_h_1=t11_l_n*t12_l_n;
t18_h=t18_h_1+t18_e;
t18_l_n=t18_l_n_1+t18_e;
}
else if(t12_l_n>=0){
REAL t18_l_n_1_1=t11_l_n*t12_h;
REAL t18_l_n_1_2=t12_l_n*t11_h;
REAL t18_l_n_1=max(t18_l_n_1_1,t18_l_n_1_2);
REAL t18_h_1_1=t11_l_n*t12_l_n;
REAL t18_h_1_2=t11_h*t12_h;
REAL t18_h_1=max(t18_h_1_1,t18_h_1_2);
t18_h=t18_h_1+t18_e;
t18_l_n=t18_l_n_1+t18_e;
}
else{
REAL t18_l_n_1=t11_l_n*t12_h;
REAL t18_h_1=t11_h*t12_h;
t18_h=t18_h_1+t18_e;
t18_l_n=t18_l_n_1+t18_e;
}
}
else{
if(t12_h<=0){
REAL t18_l_n_1=t12_l_n*t11_h;
REAL t12_h_n=-t12_h;
REAL t18_h_1=t12_h_n*t11_l_n;
t18_h=t18_h_1+t18_e;
t18_l_n=t18_l_n_1+t18_e;
}
else if(t12_l_n>=0){
REAL t18_l_n_1=t12_l_n*t11_h;
REAL t18_h_1=t11_h*t12_h;
t18_h=t18_h_1+t18_e;
t18_l_n=t18_l_n_1+t18_e;
}
else{
REAL t11_l=-t11_l_n;
REAL t18_l_n_1=t11_l*t12_l_n;
REAL t18_h_1=t11_h*t12_h;
t18_h=t18_h_1+t18_e;
t18_l_n=t18_l_n_1+t18_e;
}
}
REAL t19_e_1=fabs(t19);
REAL t19_e=eps*t19_e_1;
REAL t19_h_1=0.5*t18_h;
REAL t19_h=t19_h_1+t19_e;
REAL t19_l_n_1=0.5*t18_l_n;
REAL t19_l_n=t19_l_n_1+t19_e;
REAL t20_e_1=fabs(t20);
REAL t20_e=eps*t20_e_1;
REAL  t20_h, t20_l_n;
if(t13_h<=0){
if(t14_h<=0){
REAL t13_h_n=-t13_h;
REAL t20_l_n_1=t13_h_n*t14_h;
REAL t20_h_1=t13_l_n*t14_l_n;
t20_h=t20_h_1+t20_e;
t20_l_n=t20_l_n_1+t20_e;
}
else if(t14_l_n>=0){
REAL t20_l_n_1=t13_l_n*t14_h;
REAL t20_h_1=t13_l_n*t14_l_n;
t20_h=t20_h_1+t20_e;
t20_l_n=t20_l_n_1+t20_e;
}
else{
REAL t20_l_n_1=t13_l_n*t14_h;
REAL t13_h_n=-t13_h;
REAL t20_h_1=t13_h_n*t14_l_n;
t20_h=t20_h_1+t20_e;
t20_l_n=t20_l_n_1+t20_e;
}
}
else if(t13_l_n>=0){
if(t14_h<=0){
REAL t20_l_n_1=t13_h*t14_l_n;
REAL t20_h_1=t13_l_n*t14_l_n;
t20_h=t20_h_1+t20_e;
t20_l_n=t20_l_n_1+t20_e;
}
else if(t14_l_n>=0){
REAL t20_l_n_1_1=t13_l_n*t14_h;
REAL t20_l_n_1_2=t14_l_n*t13_h;
REAL t20_l_n_1=max(t20_l_n_1_1,t20_l_n_1_2);
REAL t20_h_1_1=t13_l_n*t14_l_n;
REAL t20_h_1_2=t13_h*t14_h;
REAL t20_h_1=max(t20_h_1_1,t20_h_1_2);
t20_h=t20_h_1+t20_e;
t20_l_n=t20_l_n_1+t20_e;
}
else{
REAL t20_l_n_1=t13_l_n*t14_h;
REAL t20_h_1=t13_h*t14_h;
t20_h=t20_h_1+t20_e;
t20_l_n=t20_l_n_1+t20_e;
}
}
else{
if(t14_h<=0){
REAL t20_l_n_1=t14_l_n*t13_h;
REAL t14_h_n=-t14_h;
REAL t20_h_1=t14_h_n*t13_l_n;
t20_h=t20_h_1+t20_e;
t20_l_n=t20_l_n_1+t20_e;
}
else if(t14_l_n>=0){
REAL t20_l_n_1=t14_l_n*t13_h;
REAL t20_h_1=t13_h*t14_h;
t20_h=t20_h_1+t20_e;
t20_l_n=t20_l_n_1+t20_e;
}
else{
REAL t13_l=-t13_l_n;
REAL t20_l_n_1=t13_l*t14_l_n;
REAL t20_h_1=t13_h*t14_h;
t20_h=t20_h_1+t20_e;
t20_l_n=t20_l_n_1+t20_e;
}
}
REAL t21_e_1=fabs(t21);
REAL t21_e=eps*t21_e_1;
REAL t21_h_1=0.5*t20_h;
REAL t21_h=t21_h_1+t21_e;
REAL t21_l_n_1=0.5*t20_l_n;
REAL t21_l_n=t21_l_n_1+t21_e;
REAL t22_e_1=fabs(t22);
REAL t22_e=eps*t22_e_1;
REAL t22_h=t22+t22_e;
REAL t22_l_n=t22_e-t22;
REAL t23_e_1=fabs(t23);
REAL t23_e=eps*t23_e_1;
REAL  t23_h, t23_l_n;
if(t15_h<=0){
if(t22_h<=0){
REAL t15_h_n=-t15_h;
REAL t23_l_n_1=t15_h_n*t22_h;
REAL t23_h_1=t15_l_n*t22_l_n;
t23_h=t23_h_1+t23_e;
t23_l_n=t23_l_n_1+t23_e;
}
else if(t22_l_n>=0){
REAL t23_l_n_1=t15_l_n*t22_h;
REAL t23_h_1=t15_l_n*t22_l_n;
t23_h=t23_h_1+t23_e;
t23_l_n=t23_l_n_1+t23_e;
}
else{
REAL t23_l_n_1=t15_l_n*t22_h;
REAL t15_h_n=-t15_h;
REAL t23_h_1=t15_h_n*t22_l_n;
t23_h=t23_h_1+t23_e;
t23_l_n=t23_l_n_1+t23_e;
}
}
else if(t15_l_n>=0){
if(t22_h<=0){
REAL t23_l_n_1=t15_h*t22_l_n;
REAL t23_h_1=t15_l_n*t22_l_n;
t23_h=t23_h_1+t23_e;
t23_l_n=t23_l_n_1+t23_e;
}
else if(t22_l_n>=0){
REAL t23_l_n_1_1=t15_l_n*t22_h;
REAL t23_l_n_1_2=t22_l_n*t15_h;
REAL t23_l_n_1=max(t23_l_n_1_1,t23_l_n_1_2);
REAL t23_h_1_1=t15_l_n*t22_l_n;
REAL t23_h_1_2=t15_h*t22_h;
REAL t23_h_1=max(t23_h_1_1,t23_h_1_2);
t23_h=t23_h_1+t23_e;
t23_l_n=t23_l_n_1+t23_e;
}
else{
REAL t23_l_n_1=t15_l_n*t22_h;
REAL t23_h_1=t15_h*t22_h;
t23_h=t23_h_1+t23_e;
t23_l_n=t23_l_n_1+t23_e;
}
}
else{
if(t22_h<=0){
REAL t23_l_n_1=t22_l_n*t15_h;
REAL t22_h_n=-t22_h;
REAL t23_h_1=t22_h_n*t15_l_n;
t23_h=t23_h_1+t23_e;
t23_l_n=t23_l_n_1+t23_e;
}
else if(t22_l_n>=0){
REAL t23_l_n_1=t22_l_n*t15_h;
REAL t23_h_1=t15_h*t22_h;
t23_h=t23_h_1+t23_e;
t23_l_n=t23_l_n_1+t23_e;
}
else{
REAL t15_l=-t15_l_n;
REAL t23_l_n_1=t15_l*t22_l_n;
REAL t23_h_1=t15_h*t22_h;
t23_h=t23_h_1+t23_e;
t23_l_n=t23_l_n_1+t23_e;
}
}
REAL t24_e_1=fabs(t24);
REAL t24_e=eps*t24_e_1;
REAL t24_h_1=0.5*t23_h;
REAL t24_h=t24_h_1+t24_e;
REAL t24_l_n_1=0.5*t23_l_n;
REAL t24_l_n=t24_l_n_1+t24_e;
REAL t25_e_1=fabs(t25);
REAL t25_e=eps*t25_e_1;
REAL t25_h=t25+t25_e;
REAL t25_l_n=t25_e-t25;
REAL t26_e_1=fabs(t26);
REAL t26_e=eps*t26_e_1;
REAL  t26_h, t26_l_n;
if(t2_h<=0){
if(t25_h<=0){
REAL t2_h_n=-t2_h;
REAL t26_l_n_1=t2_h_n*t25_h;
REAL t26_h_1=t2_l_n*t25_l_n;
t26_h=t26_h_1+t26_e;
t26_l_n=t26_l_n_1+t26_e;
}
else if(t25_l_n>=0){
REAL t26_l_n_1=t2_l_n*t25_h;
REAL t26_h_1=t2_l_n*t25_l_n;
t26_h=t26_h_1+t26_e;
t26_l_n=t26_l_n_1+t26_e;
}
else{
REAL t26_l_n_1=t2_l_n*t25_h;
REAL t2_h_n=-t2_h;
REAL t26_h_1=t2_h_n*t25_l_n;
t26_h=t26_h_1+t26_e;
t26_l_n=t26_l_n_1+t26_e;
}
}
else if(t2_l_n>=0){
if(t25_h<=0){
REAL t26_l_n_1=t2_h*t25_l_n;
REAL t26_h_1=t2_l_n*t25_l_n;
t26_h=t26_h_1+t26_e;
t26_l_n=t26_l_n_1+t26_e;
}
else if(t25_l_n>=0){
REAL t26_l_n_1_1=t2_l_n*t25_h;
REAL t26_l_n_1_2=t25_l_n*t2_h;
REAL t26_l_n_1=max(t26_l_n_1_1,t26_l_n_1_2);
REAL t26_h_1_1=t2_l_n*t25_l_n;
REAL t26_h_1_2=t2_h*t25_h;
REAL t26_h_1=max(t26_h_1_1,t26_h_1_2);
t26_h=t26_h_1+t26_e;
t26_l_n=t26_l_n_1+t26_e;
}
else{
REAL t26_l_n_1=t2_l_n*t25_h;
REAL t26_h_1=t2_h*t25_h;
t26_h=t26_h_1+t26_e;
t26_l_n=t26_l_n_1+t26_e;
}
}
else{
if(t25_h<=0){
REAL t26_l_n_1=t25_l_n*t2_h;
REAL t25_h_n=-t25_h;
REAL t26_h_1=t25_h_n*t2_l_n;
t26_h=t26_h_1+t26_e;
t26_l_n=t26_l_n_1+t26_e;
}
else if(t25_l_n>=0){
REAL t26_l_n_1=t25_l_n*t2_h;
REAL t26_h_1=t2_h*t25_h;
t26_h=t26_h_1+t26_e;
t26_l_n=t26_l_n_1+t26_e;
}
else{
REAL t2_l=-t2_l_n;
REAL t26_l_n_1=t2_l*t25_l_n;
REAL t26_h_1=t2_h*t25_h;
t26_h=t26_h_1+t26_e;
t26_l_n=t26_l_n_1+t26_e;
}
}
REAL t27_e_1=fabs(t27);
REAL t27_e=eps*t27_e_1;
REAL t27_h_1=0.5*t26_h;
REAL t27_h=t27_h_1+t27_e;
REAL t27_l_n_1=0.5*t26_l_n;
REAL t27_l_n=t27_l_n_1+t27_e;
REAL t28_e_1=fabs(t28);
REAL t28_e=eps*t28_e_1;
REAL  t28_h, t28_l_n;
if(t3_h<=0){
if(t4_h<=0){
REAL t3_h_n=-t3_h;
REAL t28_l_n_1=t3_h_n*t4_h;
REAL t28_h_1=t3_l_n*t4_l_n;
t28_h=t28_h_1+t28_e;
t28_l_n=t28_l_n_1+t28_e;
}
else if(t4_l_n>=0){
REAL t28_l_n_1=t3_l_n*t4_h;
REAL t28_h_1=t3_l_n*t4_l_n;
t28_h=t28_h_1+t28_e;
t28_l_n=t28_l_n_1+t28_e;
}
else{
REAL t28_l_n_1=t3_l_n*t4_h;
REAL t3_h_n=-t3_h;
REAL t28_h_1=t3_h_n*t4_l_n;
t28_h=t28_h_1+t28_e;
t28_l_n=t28_l_n_1+t28_e;
}
}
else if(t3_l_n>=0){
if(t4_h<=0){
REAL t28_l_n_1=t3_h*t4_l_n;
REAL t28_h_1=t3_l_n*t4_l_n;
t28_h=t28_h_1+t28_e;
t28_l_n=t28_l_n_1+t28_e;
}
else if(t4_l_n>=0){
REAL t28_l_n_1_1=t3_l_n*t4_h;
REAL t28_l_n_1_2=t4_l_n*t3_h;
REAL t28_l_n_1=max(t28_l_n_1_1,t28_l_n_1_2);
REAL t28_h_1_1=t3_l_n*t4_l_n;
REAL t28_h_1_2=t3_h*t4_h;
REAL t28_h_1=max(t28_h_1_1,t28_h_1_2);
t28_h=t28_h_1+t28_e;
t28_l_n=t28_l_n_1+t28_e;
}
else{
REAL t28_l_n_1=t3_l_n*t4_h;
REAL t28_h_1=t3_h*t4_h;
t28_h=t28_h_1+t28_e;
t28_l_n=t28_l_n_1+t28_e;
}
}
else{
if(t4_h<=0){
REAL t28_l_n_1=t4_l_n*t3_h;
REAL t4_h_n=-t4_h;
REAL t28_h_1=t4_h_n*t3_l_n;
t28_h=t28_h_1+t28_e;
t28_l_n=t28_l_n_1+t28_e;
}
else if(t4_l_n>=0){
REAL t28_l_n_1=t4_l_n*t3_h;
REAL t28_h_1=t3_h*t4_h;
t28_h=t28_h_1+t28_e;
t28_l_n=t28_l_n_1+t28_e;
}
else{
REAL t3_l=-t3_l_n;
REAL t28_l_n_1=t3_l*t4_l_n;
REAL t28_h_1=t3_h*t4_h;
t28_h=t28_h_1+t28_e;
t28_l_n=t28_l_n_1+t28_e;
}
}
REAL t29_e_1=fabs(t29);
REAL t29_e=eps*t29_e_1;
REAL t29_h_1=0.5*t28_h;
REAL t29_h=t29_h_1+t29_e;
REAL t29_l_n_1=0.5*t28_l_n;
REAL t29_l_n=t29_l_n_1+t29_e;
REAL t30_e_1=fabs(t30);
REAL t30_e=eps*t30_e_1;
REAL  t30_h, t30_l_n;
if(t5_h<=0){
if(t6_h<=0){
REAL t5_h_n=-t5_h;
REAL t30_l_n_1=t5_h_n*t6_h;
REAL t30_h_1=t5_l_n*t6_l_n;
t30_h=t30_h_1+t30_e;
t30_l_n=t30_l_n_1+t30_e;
}
else if(t6_l_n>=0){
REAL t30_l_n_1=t5_l_n*t6_h;
REAL t30_h_1=t5_l_n*t6_l_n;
t30_h=t30_h_1+t30_e;
t30_l_n=t30_l_n_1+t30_e;
}
else{
REAL t30_l_n_1=t5_l_n*t6_h;
REAL t5_h_n=-t5_h;
REAL t30_h_1=t5_h_n*t6_l_n;
t30_h=t30_h_1+t30_e;
t30_l_n=t30_l_n_1+t30_e;
}
}
else if(t5_l_n>=0){
if(t6_h<=0){
REAL t30_l_n_1=t5_h*t6_l_n;
REAL t30_h_1=t5_l_n*t6_l_n;
t30_h=t30_h_1+t30_e;
t30_l_n=t30_l_n_1+t30_e;
}
else if(t6_l_n>=0){
REAL t30_l_n_1_1=t5_l_n*t6_h;
REAL t30_l_n_1_2=t6_l_n*t5_h;
REAL t30_l_n_1=max(t30_l_n_1_1,t30_l_n_1_2);
REAL t30_h_1_1=t5_l_n*t6_l_n;
REAL t30_h_1_2=t5_h*t6_h;
REAL t30_h_1=max(t30_h_1_1,t30_h_1_2);
t30_h=t30_h_1+t30_e;
t30_l_n=t30_l_n_1+t30_e;
}
else{
REAL t30_l_n_1=t5_l_n*t6_h;
REAL t30_h_1=t5_h*t6_h;
t30_h=t30_h_1+t30_e;
t30_l_n=t30_l_n_1+t30_e;
}
}
else{
if(t6_h<=0){
REAL t30_l_n_1=t6_l_n*t5_h;
REAL t6_h_n=-t6_h;
REAL t30_h_1=t6_h_n*t5_l_n;
t30_h=t30_h_1+t30_e;
t30_l_n=t30_l_n_1+t30_e;
}
else if(t6_l_n>=0){
REAL t30_l_n_1=t6_l_n*t5_h;
REAL t30_h_1=t5_h*t6_h;
t30_h=t30_h_1+t30_e;
t30_l_n=t30_l_n_1+t30_e;
}
else{
REAL t5_l=-t5_l_n;
REAL t30_l_n_1=t5_l*t6_l_n;
REAL t30_h_1=t5_h*t6_h;
t30_h=t30_h_1+t30_e;
t30_l_n=t30_l_n_1+t30_e;
}
}
REAL t31_e_1=fabs(t31);
REAL t31_e=eps*t31_e_1;
REAL t31_h_1=0.5*t30_h;
REAL t31_h=t31_h_1+t31_e;
REAL t31_l_n_1=0.5*t30_l_n;
REAL t31_l_n=t31_l_n_1+t31_e;
REAL t32_e_1=fabs(t32);
REAL t32_e=eps*t32_e_1;
REAL  t32_h, t32_l_n;
if(t7_h<=0){
if(t8_h<=0){
REAL t7_h_n=-t7_h;
REAL t32_l_n_1=t7_h_n*t8_h;
REAL t32_h_1=t7_l_n*t8_l_n;
t32_h=t32_h_1+t32_e;
t32_l_n=t32_l_n_1+t32_e;
}
else if(t8_l_n>=0){
REAL t32_l_n_1=t7_l_n*t8_h;
REAL t32_h_1=t7_l_n*t8_l_n;
t32_h=t32_h_1+t32_e;
t32_l_n=t32_l_n_1+t32_e;
}
else{
REAL t32_l_n_1=t7_l_n*t8_h;
REAL t7_h_n=-t7_h;
REAL t32_h_1=t7_h_n*t8_l_n;
t32_h=t32_h_1+t32_e;
t32_l_n=t32_l_n_1+t32_e;
}
}
else if(t7_l_n>=0){
if(t8_h<=0){
REAL t32_l_n_1=t7_h*t8_l_n;
REAL t32_h_1=t7_l_n*t8_l_n;
t32_h=t32_h_1+t32_e;
t32_l_n=t32_l_n_1+t32_e;
}
else if(t8_l_n>=0){
REAL t32_l_n_1_1=t7_l_n*t8_h;
REAL t32_l_n_1_2=t8_l_n*t7_h;
REAL t32_l_n_1=max(t32_l_n_1_1,t32_l_n_1_2);
REAL t32_h_1_1=t7_l_n*t8_l_n;
REAL t32_h_1_2=t7_h*t8_h;
REAL t32_h_1=max(t32_h_1_1,t32_h_1_2);
t32_h=t32_h_1+t32_e;
t32_l_n=t32_l_n_1+t32_e;
}
else{
REAL t32_l_n_1=t7_l_n*t8_h;
REAL t32_h_1=t7_h*t8_h;
t32_h=t32_h_1+t32_e;
t32_l_n=t32_l_n_1+t32_e;
}
}
else{
if(t8_h<=0){
REAL t32_l_n_1=t8_l_n*t7_h;
REAL t8_h_n=-t8_h;
REAL t32_h_1=t8_h_n*t7_l_n;
t32_h=t32_h_1+t32_e;
t32_l_n=t32_l_n_1+t32_e;
}
else if(t8_l_n>=0){
REAL t32_l_n_1=t8_l_n*t7_h;
REAL t32_h_1=t7_h*t8_h;
t32_h=t32_h_1+t32_e;
t32_l_n=t32_l_n_1+t32_e;
}
else{
REAL t7_l=-t7_l_n;
REAL t32_l_n_1=t7_l*t8_l_n;
REAL t32_h_1=t7_h*t8_h;
t32_h=t32_h_1+t32_e;
t32_l_n=t32_l_n_1+t32_e;
}
}
REAL t33_e_1=fabs(t33);
REAL t33_e=eps*t33_e_1;
REAL t33_h_1=0.5*t32_h;
REAL t33_h=t33_h_1+t33_e;
REAL t33_l_n_1=0.5*t32_l_n;
REAL t33_l_n=t33_l_n_1+t33_e;
REAL t34_e_1=fabs(t34);
REAL t34_e=eps*t34_e_1;
REAL  t34_h, t34_l_n;
if(t10_h<=0){
if(t13_h<=0){
REAL t10_h_n=-t10_h;
REAL t34_l_n_1=t10_h_n*t13_h;
REAL t34_h_1=t10_l_n*t13_l_n;
t34_h=t34_h_1+t34_e;
t34_l_n=t34_l_n_1+t34_e;
}
else if(t13_l_n>=0){
REAL t34_l_n_1=t10_l_n*t13_h;
REAL t34_h_1=t10_l_n*t13_l_n;
t34_h=t34_h_1+t34_e;
t34_l_n=t34_l_n_1+t34_e;
}
else{
REAL t34_l_n_1=t10_l_n*t13_h;
REAL t10_h_n=-t10_h;
REAL t34_h_1=t10_h_n*t13_l_n;
t34_h=t34_h_1+t34_e;
t34_l_n=t34_l_n_1+t34_e;
}
}
else if(t10_l_n>=0){
if(t13_h<=0){
REAL t34_l_n_1=t10_h*t13_l_n;
REAL t34_h_1=t10_l_n*t13_l_n;
t34_h=t34_h_1+t34_e;
t34_l_n=t34_l_n_1+t34_e;
}
else if(t13_l_n>=0){
REAL t34_l_n_1_1=t10_l_n*t13_h;
REAL t34_l_n_1_2=t13_l_n*t10_h;
REAL t34_l_n_1=max(t34_l_n_1_1,t34_l_n_1_2);
REAL t34_h_1_1=t10_l_n*t13_l_n;
REAL t34_h_1_2=t10_h*t13_h;
REAL t34_h_1=max(t34_h_1_1,t34_h_1_2);
t34_h=t34_h_1+t34_e;
t34_l_n=t34_l_n_1+t34_e;
}
else{
REAL t34_l_n_1=t10_l_n*t13_h;
REAL t34_h_1=t10_h*t13_h;
t34_h=t34_h_1+t34_e;
t34_l_n=t34_l_n_1+t34_e;
}
}
else{
if(t13_h<=0){
REAL t34_l_n_1=t13_l_n*t10_h;
REAL t13_h_n=-t13_h;
REAL t34_h_1=t13_h_n*t10_l_n;
t34_h=t34_h_1+t34_e;
t34_l_n=t34_l_n_1+t34_e;
}
else if(t13_l_n>=0){
REAL t34_l_n_1=t13_l_n*t10_h;
REAL t34_h_1=t10_h*t13_h;
t34_h=t34_h_1+t34_e;
t34_l_n=t34_l_n_1+t34_e;
}
else{
REAL t10_l=-t10_l_n;
REAL t34_l_n_1=t10_l*t13_l_n;
REAL t34_h_1=t10_h*t13_h;
t34_h=t34_h_1+t34_e;
t34_l_n=t34_l_n_1+t34_e;
}
}
REAL t35_e_1=fabs(t35);
REAL t35_e=eps*t35_e_1;
REAL t35_h_1=0.5*t34_l_n;
REAL t35_h=t35_h_1+t35_e;
REAL t35_l_n_1=0.5*t34_h;
REAL t35_l_n=t35_l_n_1+t35_e;
REAL t36_e_1=fabs(t36);
REAL t36_e=eps*t36_e_1;
REAL  t36_h, t36_l_n;
if(t11_h<=0){
if(t8_h<=0){
REAL t11_h_n=-t11_h;
REAL t36_l_n_1=t11_h_n*t8_h;
REAL t36_h_1=t11_l_n*t8_l_n;
t36_h=t36_h_1+t36_e;
t36_l_n=t36_l_n_1+t36_e;
}
else if(t8_l_n>=0){
REAL t36_l_n_1=t11_l_n*t8_h;
REAL t36_h_1=t11_l_n*t8_l_n;
t36_h=t36_h_1+t36_e;
t36_l_n=t36_l_n_1+t36_e;
}
else{
REAL t36_l_n_1=t11_l_n*t8_h;
REAL t11_h_n=-t11_h;
REAL t36_h_1=t11_h_n*t8_l_n;
t36_h=t36_h_1+t36_e;
t36_l_n=t36_l_n_1+t36_e;
}
}
else if(t11_l_n>=0){
if(t8_h<=0){
REAL t36_l_n_1=t11_h*t8_l_n;
REAL t36_h_1=t11_l_n*t8_l_n;
t36_h=t36_h_1+t36_e;
t36_l_n=t36_l_n_1+t36_e;
}
else if(t8_l_n>=0){
REAL t36_l_n_1_1=t11_l_n*t8_h;
REAL t36_l_n_1_2=t8_l_n*t11_h;
REAL t36_l_n_1=max(t36_l_n_1_1,t36_l_n_1_2);
REAL t36_h_1_1=t11_l_n*t8_l_n;
REAL t36_h_1_2=t11_h*t8_h;
REAL t36_h_1=max(t36_h_1_1,t36_h_1_2);
t36_h=t36_h_1+t36_e;
t36_l_n=t36_l_n_1+t36_e;
}
else{
REAL t36_l_n_1=t11_l_n*t8_h;
REAL t36_h_1=t11_h*t8_h;
t36_h=t36_h_1+t36_e;
t36_l_n=t36_l_n_1+t36_e;
}
}
else{
if(t8_h<=0){
REAL t36_l_n_1=t8_l_n*t11_h;
REAL t8_h_n=-t8_h;
REAL t36_h_1=t8_h_n*t11_l_n;
t36_h=t36_h_1+t36_e;
t36_l_n=t36_l_n_1+t36_e;
}
else if(t8_l_n>=0){
REAL t36_l_n_1=t8_l_n*t11_h;
REAL t36_h_1=t11_h*t8_h;
t36_h=t36_h_1+t36_e;
t36_l_n=t36_l_n_1+t36_e;
}
else{
REAL t11_l=-t11_l_n;
REAL t36_l_n_1=t11_l*t8_l_n;
REAL t36_h_1=t11_h*t8_h;
t36_h=t36_h_1+t36_e;
t36_l_n=t36_l_n_1+t36_e;
}
}
REAL t37_e_1=fabs(t37);
REAL t37_e=eps*t37_e_1;
REAL t37_h_1=0.5*t36_l_n;
REAL t37_h=t37_h_1+t37_e;
REAL t37_l_n_1=0.5*t36_h;
REAL t37_l_n=t37_l_n_1+t37_e;
REAL t38_e_1=fabs(t38);
REAL t38_e=eps*t38_e_1;
REAL  t38_h, t38_l_n;
if(t12_h<=0){
if(t15_h<=0){
REAL t12_h_n=-t12_h;
REAL t38_l_n_1=t12_h_n*t15_h;
REAL t38_h_1=t12_l_n*t15_l_n;
t38_h=t38_h_1+t38_e;
t38_l_n=t38_l_n_1+t38_e;
}
else if(t15_l_n>=0){
REAL t38_l_n_1=t12_l_n*t15_h;
REAL t38_h_1=t12_l_n*t15_l_n;
t38_h=t38_h_1+t38_e;
t38_l_n=t38_l_n_1+t38_e;
}
else{
REAL t38_l_n_1=t12_l_n*t15_h;
REAL t12_h_n=-t12_h;
REAL t38_h_1=t12_h_n*t15_l_n;
t38_h=t38_h_1+t38_e;
t38_l_n=t38_l_n_1+t38_e;
}
}
else if(t12_l_n>=0){
if(t15_h<=0){
REAL t38_l_n_1=t12_h*t15_l_n;
REAL t38_h_1=t12_l_n*t15_l_n;
t38_h=t38_h_1+t38_e;
t38_l_n=t38_l_n_1+t38_e;
}
else if(t15_l_n>=0){
REAL t38_l_n_1_1=t12_l_n*t15_h;
REAL t38_l_n_1_2=t15_l_n*t12_h;
REAL t38_l_n_1=max(t38_l_n_1_1,t38_l_n_1_2);
REAL t38_h_1_1=t12_l_n*t15_l_n;
REAL t38_h_1_2=t12_h*t15_h;
REAL t38_h_1=max(t38_h_1_1,t38_h_1_2);
t38_h=t38_h_1+t38_e;
t38_l_n=t38_l_n_1+t38_e;
}
else{
REAL t38_l_n_1=t12_l_n*t15_h;
REAL t38_h_1=t12_h*t15_h;
t38_h=t38_h_1+t38_e;
t38_l_n=t38_l_n_1+t38_e;
}
}
else{
if(t15_h<=0){
REAL t38_l_n_1=t15_l_n*t12_h;
REAL t15_h_n=-t15_h;
REAL t38_h_1=t15_h_n*t12_l_n;
t38_h=t38_h_1+t38_e;
t38_l_n=t38_l_n_1+t38_e;
}
else if(t15_l_n>=0){
REAL t38_l_n_1=t15_l_n*t12_h;
REAL t38_h_1=t12_h*t15_h;
t38_h=t38_h_1+t38_e;
t38_l_n=t38_l_n_1+t38_e;
}
else{
REAL t12_l=-t12_l_n;
REAL t38_l_n_1=t12_l*t15_l_n;
REAL t38_h_1=t12_h*t15_h;
t38_h=t38_h_1+t38_e;
t38_l_n=t38_l_n_1+t38_e;
}
}
REAL t39_e_1=fabs(t39);
REAL t39_e=eps*t39_e_1;
REAL t39_h_1=0.5*t38_l_n;
REAL t39_h=t39_h_1+t39_e;
REAL t39_l_n_1=0.5*t38_h;
REAL t39_l_n=t39_l_n_1+t39_e;
REAL t40_e_1=fabs(t40);
REAL t40_e=eps*t40_e_1;
REAL t40_h=t40+t40_e;
REAL t40_l_n=t40_e-t40;
REAL t41_e_1=fabs(t41);
REAL t41_e=eps*t41_e_1;
REAL  t41_h, t41_l_n;
if(t14_h<=0){
if(t40_h<=0){
REAL t14_h_n=-t14_h;
REAL t41_l_n_1=t14_h_n*t40_h;
REAL t41_h_1=t14_l_n*t40_l_n;
t41_h=t41_h_1+t41_e;
t41_l_n=t41_l_n_1+t41_e;
}
else if(t40_l_n>=0){
REAL t41_l_n_1=t14_l_n*t40_h;
REAL t41_h_1=t14_l_n*t40_l_n;
t41_h=t41_h_1+t41_e;
t41_l_n=t41_l_n_1+t41_e;
}
else{
REAL t41_l_n_1=t14_l_n*t40_h;
REAL t14_h_n=-t14_h;
REAL t41_h_1=t14_h_n*t40_l_n;
t41_h=t41_h_1+t41_e;
t41_l_n=t41_l_n_1+t41_e;
}
}
else if(t14_l_n>=0){
if(t40_h<=0){
REAL t41_l_n_1=t14_h*t40_l_n;
REAL t41_h_1=t14_l_n*t40_l_n;
t41_h=t41_h_1+t41_e;
t41_l_n=t41_l_n_1+t41_e;
}
else if(t40_l_n>=0){
REAL t41_l_n_1_1=t14_l_n*t40_h;
REAL t41_l_n_1_2=t40_l_n*t14_h;
REAL t41_l_n_1=max(t41_l_n_1_1,t41_l_n_1_2);
REAL t41_h_1_1=t14_l_n*t40_l_n;
REAL t41_h_1_2=t14_h*t40_h;
REAL t41_h_1=max(t41_h_1_1,t41_h_1_2);
t41_h=t41_h_1+t41_e;
t41_l_n=t41_l_n_1+t41_e;
}
else{
REAL t41_l_n_1=t14_l_n*t40_h;
REAL t41_h_1=t14_h*t40_h;
t41_h=t41_h_1+t41_e;
t41_l_n=t41_l_n_1+t41_e;
}
}
else{
if(t40_h<=0){
REAL t41_l_n_1=t40_l_n*t14_h;
REAL t40_h_n=-t40_h;
REAL t41_h_1=t40_h_n*t14_l_n;
t41_h=t41_h_1+t41_e;
t41_l_n=t41_l_n_1+t41_e;
}
else if(t40_l_n>=0){
REAL t41_l_n_1=t40_l_n*t14_h;
REAL t41_h_1=t14_h*t40_h;
t41_h=t41_h_1+t41_e;
t41_l_n=t41_l_n_1+t41_e;
}
else{
REAL t14_l=-t14_l_n;
REAL t41_l_n_1=t14_l*t40_l_n;
REAL t41_h_1=t14_h*t40_h;
t41_h=t41_h_1+t41_e;
t41_l_n=t41_l_n_1+t41_e;
}
}
REAL t42_e_1=fabs(t42);
REAL t42_e=eps*t42_e_1;
REAL t42_h_1=0.5*t41_l_n;
REAL t42_h=t42_h_1+t42_e;
REAL t42_l_n_1=0.5*t41_h;
REAL t42_l_n=t42_l_n_1+t42_e;
REAL t43_e_1=fabs(t43);
REAL t43_e=eps*t43_e_1;
REAL  t43_h, t43_l_n;
if(t2_h<=0){
if(t5_h<=0){
REAL t2_h_n=-t2_h;
REAL t43_l_n_1=t2_h_n*t5_h;
REAL t43_h_1=t2_l_n*t5_l_n;
t43_h=t43_h_1+t43_e;
t43_l_n=t43_l_n_1+t43_e;
}
else if(t5_l_n>=0){
REAL t43_l_n_1=t2_l_n*t5_h;
REAL t43_h_1=t2_l_n*t5_l_n;
t43_h=t43_h_1+t43_e;
t43_l_n=t43_l_n_1+t43_e;
}
else{
REAL t43_l_n_1=t2_l_n*t5_h;
REAL t2_h_n=-t2_h;
REAL t43_h_1=t2_h_n*t5_l_n;
t43_h=t43_h_1+t43_e;
t43_l_n=t43_l_n_1+t43_e;
}
}
else if(t2_l_n>=0){
if(t5_h<=0){
REAL t43_l_n_1=t2_h*t5_l_n;
REAL t43_h_1=t2_l_n*t5_l_n;
t43_h=t43_h_1+t43_e;
t43_l_n=t43_l_n_1+t43_e;
}
else if(t5_l_n>=0){
REAL t43_l_n_1_1=t2_l_n*t5_h;
REAL t43_l_n_1_2=t5_l_n*t2_h;
REAL t43_l_n_1=max(t43_l_n_1_1,t43_l_n_1_2);
REAL t43_h_1_1=t2_l_n*t5_l_n;
REAL t43_h_1_2=t2_h*t5_h;
REAL t43_h_1=max(t43_h_1_1,t43_h_1_2);
t43_h=t43_h_1+t43_e;
t43_l_n=t43_l_n_1+t43_e;
}
else{
REAL t43_l_n_1=t2_l_n*t5_h;
REAL t43_h_1=t2_h*t5_h;
t43_h=t43_h_1+t43_e;
t43_l_n=t43_l_n_1+t43_e;
}
}
else{
if(t5_h<=0){
REAL t43_l_n_1=t5_l_n*t2_h;
REAL t5_h_n=-t5_h;
REAL t43_h_1=t5_h_n*t2_l_n;
t43_h=t43_h_1+t43_e;
t43_l_n=t43_l_n_1+t43_e;
}
else if(t5_l_n>=0){
REAL t43_l_n_1=t5_l_n*t2_h;
REAL t43_h_1=t2_h*t5_h;
t43_h=t43_h_1+t43_e;
t43_l_n=t43_l_n_1+t43_e;
}
else{
REAL t2_l=-t2_l_n;
REAL t43_l_n_1=t2_l*t5_l_n;
REAL t43_h_1=t2_h*t5_h;
t43_h=t43_h_1+t43_e;
t43_l_n=t43_l_n_1+t43_e;
}
}
REAL t44_e_1=fabs(t44);
REAL t44_e=eps*t44_e_1;
REAL t44_h_1=0.5*t43_l_n;
REAL t44_h=t44_h_1+t44_e;
REAL t44_l_n_1=0.5*t43_h;
REAL t44_l_n=t44_l_n_1+t44_e;
REAL t45_e_1=fabs(t45);
REAL t45_e=eps*t45_e_1;
REAL t45_h=t45+t45_e;
REAL t45_l_n=t45_e-t45;
REAL t46_e_1=fabs(t46);
REAL t46_e=eps*t46_e_1;
REAL  t46_h, t46_l_n;
if(t3_h<=0){
if(t45_h<=0){
REAL t3_h_n=-t3_h;
REAL t46_l_n_1=t3_h_n*t45_h;
REAL t46_h_1=t3_l_n*t45_l_n;
t46_h=t46_h_1+t46_e;
t46_l_n=t46_l_n_1+t46_e;
}
else if(t45_l_n>=0){
REAL t46_l_n_1=t3_l_n*t45_h;
REAL t46_h_1=t3_l_n*t45_l_n;
t46_h=t46_h_1+t46_e;
t46_l_n=t46_l_n_1+t46_e;
}
else{
REAL t46_l_n_1=t3_l_n*t45_h;
REAL t3_h_n=-t3_h;
REAL t46_h_1=t3_h_n*t45_l_n;
t46_h=t46_h_1+t46_e;
t46_l_n=t46_l_n_1+t46_e;
}
}
else if(t3_l_n>=0){
if(t45_h<=0){
REAL t46_l_n_1=t3_h*t45_l_n;
REAL t46_h_1=t3_l_n*t45_l_n;
t46_h=t46_h_1+t46_e;
t46_l_n=t46_l_n_1+t46_e;
}
else if(t45_l_n>=0){
REAL t46_l_n_1_1=t3_l_n*t45_h;
REAL t46_l_n_1_2=t45_l_n*t3_h;
REAL t46_l_n_1=max(t46_l_n_1_1,t46_l_n_1_2);
REAL t46_h_1_1=t3_l_n*t45_l_n;
REAL t46_h_1_2=t3_h*t45_h;
REAL t46_h_1=max(t46_h_1_1,t46_h_1_2);
t46_h=t46_h_1+t46_e;
t46_l_n=t46_l_n_1+t46_e;
}
else{
REAL t46_l_n_1=t3_l_n*t45_h;
REAL t46_h_1=t3_h*t45_h;
t46_h=t46_h_1+t46_e;
t46_l_n=t46_l_n_1+t46_e;
}
}
else{
if(t45_h<=0){
REAL t46_l_n_1=t45_l_n*t3_h;
REAL t45_h_n=-t45_h;
REAL t46_h_1=t45_h_n*t3_l_n;
t46_h=t46_h_1+t46_e;
t46_l_n=t46_l_n_1+t46_e;
}
else if(t45_l_n>=0){
REAL t46_l_n_1=t45_l_n*t3_h;
REAL t46_h_1=t3_h*t45_h;
t46_h=t46_h_1+t46_e;
t46_l_n=t46_l_n_1+t46_e;
}
else{
REAL t3_l=-t3_l_n;
REAL t46_l_n_1=t3_l*t45_l_n;
REAL t46_h_1=t3_h*t45_h;
t46_h=t46_h_1+t46_e;
t46_l_n=t46_l_n_1+t46_e;
}
}
REAL t47_e_1=fabs(t47);
REAL t47_e=eps*t47_e_1;
REAL t47_h_1=0.5*t46_l_n;
REAL t47_h=t47_h_1+t47_e;
REAL t47_l_n_1=0.5*t46_h;
REAL t47_l_n=t47_l_n_1+t47_e;
REAL t48_e_1=fabs(t48);
REAL t48_e=eps*t48_e_1;
REAL  t48_h, t48_l_n;
if(t4_h<=0){
if(t7_h<=0){
REAL t4_h_n=-t4_h;
REAL t48_l_n_1=t4_h_n*t7_h;
REAL t48_h_1=t4_l_n*t7_l_n;
t48_h=t48_h_1+t48_e;
t48_l_n=t48_l_n_1+t48_e;
}
else if(t7_l_n>=0){
REAL t48_l_n_1=t4_l_n*t7_h;
REAL t48_h_1=t4_l_n*t7_l_n;
t48_h=t48_h_1+t48_e;
t48_l_n=t48_l_n_1+t48_e;
}
else{
REAL t48_l_n_1=t4_l_n*t7_h;
REAL t4_h_n=-t4_h;
REAL t48_h_1=t4_h_n*t7_l_n;
t48_h=t48_h_1+t48_e;
t48_l_n=t48_l_n_1+t48_e;
}
}
else if(t4_l_n>=0){
if(t7_h<=0){
REAL t48_l_n_1=t4_h*t7_l_n;
REAL t48_h_1=t4_l_n*t7_l_n;
t48_h=t48_h_1+t48_e;
t48_l_n=t48_l_n_1+t48_e;
}
else if(t7_l_n>=0){
REAL t48_l_n_1_1=t4_l_n*t7_h;
REAL t48_l_n_1_2=t7_l_n*t4_h;
REAL t48_l_n_1=max(t48_l_n_1_1,t48_l_n_1_2);
REAL t48_h_1_1=t4_l_n*t7_l_n;
REAL t48_h_1_2=t4_h*t7_h;
REAL t48_h_1=max(t48_h_1_1,t48_h_1_2);
t48_h=t48_h_1+t48_e;
t48_l_n=t48_l_n_1+t48_e;
}
else{
REAL t48_l_n_1=t4_l_n*t7_h;
REAL t48_h_1=t4_h*t7_h;
t48_h=t48_h_1+t48_e;
t48_l_n=t48_l_n_1+t48_e;
}
}
else{
if(t7_h<=0){
REAL t48_l_n_1=t7_l_n*t4_h;
REAL t7_h_n=-t7_h;
REAL t48_h_1=t7_h_n*t4_l_n;
t48_h=t48_h_1+t48_e;
t48_l_n=t48_l_n_1+t48_e;
}
else if(t7_l_n>=0){
REAL t48_l_n_1=t7_l_n*t4_h;
REAL t48_h_1=t4_h*t7_h;
t48_h=t48_h_1+t48_e;
t48_l_n=t48_l_n_1+t48_e;
}
else{
REAL t4_l=-t4_l_n;
REAL t48_l_n_1=t4_l*t7_l_n;
REAL t48_h_1=t4_h*t7_h;
t48_h=t48_h_1+t48_e;
t48_l_n=t48_l_n_1+t48_e;
}
}
REAL t49_e_1=fabs(t49);
REAL t49_e=eps*t49_e_1;
REAL t49_h_1=0.5*t48_l_n;
REAL t49_h=t49_h_1+t49_e;
REAL t49_l_n_1=0.5*t48_h;
REAL t49_l_n=t49_l_n_1+t49_e;
REAL t50_e_1=fabs(t50);
REAL t50_e=eps*t50_e_1;
REAL  t50_h, t50_l_n;
if(t6_h<=0){
if(t9_h<=0){
REAL t6_h_n=-t6_h;
REAL t50_l_n_1=t6_h_n*t9_h;
REAL t50_h_1=t6_l_n*t9_l_n;
t50_h=t50_h_1+t50_e;
t50_l_n=t50_l_n_1+t50_e;
}
else if(t9_l_n>=0){
REAL t50_l_n_1=t6_l_n*t9_h;
REAL t50_h_1=t6_l_n*t9_l_n;
t50_h=t50_h_1+t50_e;
t50_l_n=t50_l_n_1+t50_e;
}
else{
REAL t50_l_n_1=t6_l_n*t9_h;
REAL t6_h_n=-t6_h;
REAL t50_h_1=t6_h_n*t9_l_n;
t50_h=t50_h_1+t50_e;
t50_l_n=t50_l_n_1+t50_e;
}
}
else if(t6_l_n>=0){
if(t9_h<=0){
REAL t50_l_n_1=t6_h*t9_l_n;
REAL t50_h_1=t6_l_n*t9_l_n;
t50_h=t50_h_1+t50_e;
t50_l_n=t50_l_n_1+t50_e;
}
else if(t9_l_n>=0){
REAL t50_l_n_1_1=t6_l_n*t9_h;
REAL t50_l_n_1_2=t9_l_n*t6_h;
REAL t50_l_n_1=max(t50_l_n_1_1,t50_l_n_1_2);
REAL t50_h_1_1=t6_l_n*t9_l_n;
REAL t50_h_1_2=t6_h*t9_h;
REAL t50_h_1=max(t50_h_1_1,t50_h_1_2);
t50_h=t50_h_1+t50_e;
t50_l_n=t50_l_n_1+t50_e;
}
else{
REAL t50_l_n_1=t6_l_n*t9_h;
REAL t50_h_1=t6_h*t9_h;
t50_h=t50_h_1+t50_e;
t50_l_n=t50_l_n_1+t50_e;
}
}
else{
if(t9_h<=0){
REAL t50_l_n_1=t9_l_n*t6_h;
REAL t9_h_n=-t9_h;
REAL t50_h_1=t9_h_n*t6_l_n;
t50_h=t50_h_1+t50_e;
t50_l_n=t50_l_n_1+t50_e;
}
else if(t9_l_n>=0){
REAL t50_l_n_1=t9_l_n*t6_h;
REAL t50_h_1=t6_h*t9_h;
t50_h=t50_h_1+t50_e;
t50_l_n=t50_l_n_1+t50_e;
}
else{
REAL t6_l=-t6_l_n;
REAL t50_l_n_1=t6_l*t9_l_n;
REAL t50_h_1=t6_h*t9_h;
t50_h=t50_h_1+t50_e;
t50_l_n=t50_l_n_1+t50_e;
}
}
REAL t51_e_1=fabs(t51);
REAL t51_e=eps*t51_e_1;
REAL t51_h_1=0.5*t50_l_n;
REAL t51_h=t51_h_1+t51_e;
REAL t51_l_n_1=0.5*t50_h;
REAL t51_l_n=t51_l_n_1+t51_e;
REAL t52_e_1=fabs(t52);
REAL t52_e=eps*t52_e_1;
REAL t52_h_1=t17_h+t19_h;
REAL t52_h=t52_h_1+t52_e;
REAL t52_l_n_1=t17_l_n+t19_l_n;
REAL t52_l_n=t52_l_n_1+t52_e;
REAL t53_e_1=fabs(t53);
REAL t53_e=eps*t53_e_1;
REAL t53_h_1=t21_h+t52_h;
REAL t53_h=t53_h_1+t53_e;
REAL t53_l_n_1=t21_l_n+t52_l_n;
REAL t53_l_n=t53_l_n_1+t53_e;
REAL t54_e_1=fabs(t54);
REAL t54_e=eps*t54_e_1;
REAL t54_h_1=t24_h+t53_h;
REAL t54_h=t54_h_1+t54_e;
REAL t54_l_n_1=t24_l_n+t53_l_n;
REAL t54_l_n=t54_l_n_1+t54_e;
REAL t55_e_1=fabs(t55);
REAL t55_e=eps*t55_e_1;
REAL t55_h_1=t27_h+t54_h;
REAL t55_h=t55_h_1+t55_e;
REAL t55_l_n_1=t27_l_n+t54_l_n;
REAL t55_l_n=t55_l_n_1+t55_e;
REAL t56_e_1=fabs(t56);
REAL t56_e=eps*t56_e_1;
REAL t56_h_1=t29_h+t55_h;
REAL t56_h=t56_h_1+t56_e;
REAL t56_l_n_1=t29_l_n+t55_l_n;
REAL t56_l_n=t56_l_n_1+t56_e;
REAL t57_e_1=fabs(t57);
REAL t57_e=eps*t57_e_1;
REAL t57_h_1=t31_h+t56_h;
REAL t57_h=t57_h_1+t57_e;
REAL t57_l_n_1=t31_l_n+t56_l_n;
REAL t57_l_n=t57_l_n_1+t57_e;
REAL t58_e_1=fabs(t58);
REAL t58_e=eps*t58_e_1;
REAL t58_h_1=t33_h+t57_h;
REAL t58_h=t58_h_1+t58_e;
REAL t58_l_n_1=t33_l_n+t57_l_n;
REAL t58_l_n=t58_l_n_1+t58_e;
REAL t59_e_1=fabs(t59);
REAL t59_e=eps*t59_e_1;
REAL t59_h_1=t35_h+t58_h;
REAL t59_h=t59_h_1+t59_e;
REAL t59_l_n_1=t35_l_n+t58_l_n;
REAL t59_l_n=t59_l_n_1+t59_e;
REAL t60_e_1=fabs(t60);
REAL t60_e=eps*t60_e_1;
REAL t60_h_1=t37_h+t59_h;
REAL t60_h=t60_h_1+t60_e;
REAL t60_l_n_1=t37_l_n+t59_l_n;
REAL t60_l_n=t60_l_n_1+t60_e;
REAL t61_e_1=fabs(t61);
REAL t61_e=eps*t61_e_1;
REAL t61_h_1=t39_h+t60_h;
REAL t61_h=t61_h_1+t61_e;
REAL t61_l_n_1=t39_l_n+t60_l_n;
REAL t61_l_n=t61_l_n_1+t61_e;
REAL t62_e_1=fabs(t62);
REAL t62_e=eps*t62_e_1;
REAL t62_h_1=t42_h+t61_h;
REAL t62_h=t62_h_1+t62_e;
REAL t62_l_n_1=t42_l_n+t61_l_n;
REAL t62_l_n=t62_l_n_1+t62_e;
REAL t63_e_1=fabs(t63);
REAL t63_e=eps*t63_e_1;
REAL t63_h_1=t44_h+t62_h;
REAL t63_h=t63_h_1+t63_e;
REAL t63_l_n_1=t44_l_n+t62_l_n;
REAL t63_l_n=t63_l_n_1+t63_e;
REAL t64_e_1=fabs(t64);
REAL t64_e=eps*t64_e_1;
REAL t64_h_1=t47_h+t63_h;
REAL t64_h=t64_h_1+t64_e;
REAL t64_l_n_1=t47_l_n+t63_l_n;
REAL t64_l_n=t64_l_n_1+t64_e;
REAL t65_e_1=fabs(t65);
REAL t65_e=eps*t65_e_1;
REAL t65_h_1=t49_h+t64_h;
REAL t65_h=t65_h_1+t65_e;
REAL t65_l_n_1=t49_l_n+t64_l_n;
REAL t65_l_n=t65_l_n_1+t65_e;
REAL t66_e_1=fabs(t66);
REAL t66_e=eps*t66_e_1;
REAL t66_h_1=t51_h+t65_h;
REAL t66_h=t66_h_1+t66_e;
REAL t66_l_n_1=t51_l_n+t65_l_n;
REAL t66_l_n=t66_l_n_1+t66_e;
fesetround(roundmode);
*result=t66;
if (inferr==0) {
*hbound=t66_h;
*lbound=-t66_l_n;
}
else{
*hbound=1/0;
*lbound=1/0;
}
}
//Exact Function
void Exact_ConvexHullArea(REAL x0,REAL x1,REAL x2,REAL x3,REAL x4,REAL x5,REAL x6,REAL x7,REAL x8,REAL x9,REAL y0,REAL y1,REAL y2,REAL y3,REAL y4,REAL y5,REAL y6,REAL y7,REAL y8,REAL y9){
mpfr_t MPx0;
mpfr_init2(MPx0,1024);
mpfr_t MPx1;
mpfr_init2(MPx1,1024);
mpfr_t MPx2;
mpfr_init2(MPx2,1024);
mpfr_t MPx3;
mpfr_init2(MPx3,1024);
mpfr_t MPx4;
mpfr_init2(MPx4,1024);
mpfr_t MPx5;
mpfr_init2(MPx5,1024);
mpfr_t MPx6;
mpfr_init2(MPx6,1024);
mpfr_t MPx7;
mpfr_init2(MPx7,1024);
mpfr_t MPx8;
mpfr_init2(MPx8,1024);
mpfr_t MPx9;
mpfr_init2(MPx9,1024);
mpfr_t MPy0;
mpfr_init2(MPy0,1024);
mpfr_t MPy1;
mpfr_init2(MPy1,1024);
mpfr_t MPy2;
mpfr_init2(MPy2,1024);
mpfr_t MPy3;
mpfr_init2(MPy3,1024);
mpfr_t MPy4;
mpfr_init2(MPy4,1024);
mpfr_t MPy5;
mpfr_init2(MPy5,1024);
mpfr_t MPy6;
mpfr_init2(MPy6,1024);
mpfr_t MPy7;
mpfr_init2(MPy7,1024);
mpfr_t MPy8;
mpfr_init2(MPy8,1024);
mpfr_t MPy9;
mpfr_init2(MPy9,1024);
mpfr_t MPt2;
mpfr_init2(MPt2,1024);
mpfr_sub(MPt2,MPy1,MPy9,MPFR_RNDN);
mpfr_t MPt3;
mpfr_init2(MPt3,1024);
mpfr_sub(MPt3,MPx1,MPx9,MPFR_RNDN);
mpfr_t MPt4;
mpfr_init2(MPt4,1024);
mpfr_sub(MPt4,MPy2,MPy9,MPFR_RNDN);
mpfr_t MPt5;
mpfr_init2(MPt5,1024);
mpfr_sub(MPt5,MPx2,MPx9,MPFR_RNDN);
mpfr_t MPt6;
mpfr_init2(MPt6,1024);
mpfr_sub(MPt6,MPy3,MPy9,MPFR_RNDN);
mpfr_t MPt7;
mpfr_init2(MPt7,1024);
mpfr_sub(MPt7,MPx3,MPx9,MPFR_RNDN);
mpfr_t MPt8;
mpfr_init2(MPt8,1024);
mpfr_sub(MPt8,MPy4,MPy9,MPFR_RNDN);
mpfr_t MPt9;
mpfr_init2(MPt9,1024);
mpfr_sub(MPt9,MPx4,MPx9,MPFR_RNDN);
mpfr_t MPt10;
mpfr_init2(MPt10,1024);
mpfr_sub(MPt10,MPy5,MPy9,MPFR_RNDN);
mpfr_t MPt11;
mpfr_init2(MPt11,1024);
mpfr_sub(MPt11,MPx5,MPx9,MPFR_RNDN);
mpfr_t MPt12;
mpfr_init2(MPt12,1024);
mpfr_sub(MPt12,MPy6,MPy9,MPFR_RNDN);
mpfr_t MPt13;
mpfr_init2(MPt13,1024);
mpfr_sub(MPt13,MPx6,MPx9,MPFR_RNDN);
mpfr_t MPt14;
mpfr_init2(MPt14,1024);
mpfr_sub(MPt14,MPy7,MPy9,MPFR_RNDN);
mpfr_t MPt15;
mpfr_init2(MPt15,1024);
mpfr_sub(MPt15,MPx7,MPx9,MPFR_RNDN);
mpfr_t MPt16;
mpfr_init2(MPt16,1024);
mpfr_mul(MPt16,MPt10,MPt9,MPFR_RNDN);
mpfr_t MPt17;
mpfr_init2(MPt17,1024);
mpfr_mul_d(MPt17,MPt16,0.5,MPFR_RNDN);
mpfr_t MPt18;
mpfr_init2(MPt18,1024);
mpfr_mul(MPt18,MPt11,MPt12,MPFR_RNDN);
mpfr_t MPt19;
mpfr_init2(MPt19,1024);
mpfr_mul_d(MPt19,MPt18,0.5,MPFR_RNDN);
mpfr_t MPt20;
mpfr_init2(MPt20,1024);
mpfr_mul(MPt20,MPt13,MPt14,MPFR_RNDN);
mpfr_t MPt21;
mpfr_init2(MPt21,1024);
mpfr_mul_d(MPt21,MPt20,0.5,MPFR_RNDN);
mpfr_t MPt22;
mpfr_init2(MPt22,1024);
mpfr_sub(MPt22,MPy8,MPy9,MPFR_RNDN);
mpfr_t MPt23;
mpfr_init2(MPt23,1024);
mpfr_mul(MPt23,MPt15,MPt22,MPFR_RNDN);
mpfr_t MPt24;
mpfr_init2(MPt24,1024);
mpfr_mul_d(MPt24,MPt23,0.5,MPFR_RNDN);
mpfr_t MPt25;
mpfr_init2(MPt25,1024);
mpfr_sub(MPt25,MPx0,MPx9,MPFR_RNDN);
mpfr_t MPt26;
mpfr_init2(MPt26,1024);
mpfr_mul(MPt26,MPt2,MPt25,MPFR_RNDN);
mpfr_t MPt27;
mpfr_init2(MPt27,1024);
mpfr_mul_d(MPt27,MPt26,0.5,MPFR_RNDN);
mpfr_t MPt28;
mpfr_init2(MPt28,1024);
mpfr_mul(MPt28,MPt3,MPt4,MPFR_RNDN);
mpfr_t MPt29;
mpfr_init2(MPt29,1024);
mpfr_mul_d(MPt29,MPt28,0.5,MPFR_RNDN);
mpfr_t MPt30;
mpfr_init2(MPt30,1024);
mpfr_mul(MPt30,MPt5,MPt6,MPFR_RNDN);
mpfr_t MPt31;
mpfr_init2(MPt31,1024);
mpfr_mul_d(MPt31,MPt30,0.5,MPFR_RNDN);
mpfr_t MPt32;
mpfr_init2(MPt32,1024);
mpfr_mul(MPt32,MPt7,MPt8,MPFR_RNDN);
mpfr_t MPt33;
mpfr_init2(MPt33,1024);
mpfr_mul_d(MPt33,MPt32,0.5,MPFR_RNDN);
mpfr_t MPt34;
mpfr_init2(MPt34,1024);
mpfr_mul(MPt34,MPt10,MPt13,MPFR_RNDN);
mpfr_t MPt35;
mpfr_init2(MPt35,1024);
mpfr_mul_d(MPt35,MPt34,-0.5,MPFR_RNDN);
mpfr_t MPt36;
mpfr_init2(MPt36,1024);
mpfr_mul(MPt36,MPt11,MPt8,MPFR_RNDN);
mpfr_t MPt37;
mpfr_init2(MPt37,1024);
mpfr_mul_d(MPt37,MPt36,-0.5,MPFR_RNDN);
mpfr_t MPt38;
mpfr_init2(MPt38,1024);
mpfr_mul(MPt38,MPt12,MPt15,MPFR_RNDN);
mpfr_t MPt39;
mpfr_init2(MPt39,1024);
mpfr_mul_d(MPt39,MPt38,-0.5,MPFR_RNDN);
mpfr_t MPt40;
mpfr_init2(MPt40,1024);
mpfr_sub(MPt40,MPx8,MPx9,MPFR_RNDN);
mpfr_t MPt41;
mpfr_init2(MPt41,1024);
mpfr_mul(MPt41,MPt14,MPt40,MPFR_RNDN);
mpfr_t MPt42;
mpfr_init2(MPt42,1024);
mpfr_mul_d(MPt42,MPt41,-0.5,MPFR_RNDN);
mpfr_t MPt43;
mpfr_init2(MPt43,1024);
mpfr_mul(MPt43,MPt2,MPt5,MPFR_RNDN);
mpfr_t MPt44;
mpfr_init2(MPt44,1024);
mpfr_mul_d(MPt44,MPt43,-0.5,MPFR_RNDN);
mpfr_t MPt45;
mpfr_init2(MPt45,1024);
mpfr_sub(MPt45,MPy0,MPy9,MPFR_RNDN);
mpfr_t MPt46;
mpfr_init2(MPt46,1024);
mpfr_mul(MPt46,MPt3,MPt45,MPFR_RNDN);
mpfr_t MPt47;
mpfr_init2(MPt47,1024);
mpfr_mul_d(MPt47,MPt46,-0.5,MPFR_RNDN);
mpfr_t MPt48;
mpfr_init2(MPt48,1024);
mpfr_mul(MPt48,MPt4,MPt7,MPFR_RNDN);
mpfr_t MPt49;
mpfr_init2(MPt49,1024);
mpfr_mul_d(MPt49,MPt48,-0.5,MPFR_RNDN);
mpfr_t MPt50;
mpfr_init2(MPt50,1024);
mpfr_mul(MPt50,MPt6,MPt9,MPFR_RNDN);
mpfr_t MPt51;
mpfr_init2(MPt51,1024);
mpfr_mul_d(MPt51,MPt50,-0.5,MPFR_RNDN);
mpfr_t MPt52;
mpfr_init2(MPt52,1024);
mpfr_add(MPt52,MPt17,MPt19,MPFR_RNDN);
mpfr_t MPt53;
mpfr_init2(MPt53,1024);
mpfr_add(MPt53,MPt21,MPt52,MPFR_RNDN);
mpfr_t MPt54;
mpfr_init2(MPt54,1024);
mpfr_add(MPt54,MPt24,MPt53,MPFR_RNDN);
mpfr_t MPt55;
mpfr_init2(MPt55,1024);
mpfr_add(MPt55,MPt27,MPt54,MPFR_RNDN);
mpfr_t MPt56;
mpfr_init2(MPt56,1024);
mpfr_add(MPt56,MPt29,MPt55,MPFR_RNDN);
mpfr_t MPt57;
mpfr_init2(MPt57,1024);
mpfr_add(MPt57,MPt31,MPt56,MPFR_RNDN);
mpfr_t MPt58;
mpfr_init2(MPt58,1024);
mpfr_add(MPt58,MPt33,MPt57,MPFR_RNDN);
mpfr_t MPt59;
mpfr_init2(MPt59,1024);
mpfr_add(MPt59,MPt35,MPt58,MPFR_RNDN);
mpfr_t MPt60;
mpfr_init2(MPt60,1024);
mpfr_add(MPt60,MPt37,MPt59,MPFR_RNDN);
mpfr_t MPt61;
mpfr_init2(MPt61,1024);
mpfr_add(MPt61,MPt39,MPt60,MPFR_RNDN);
mpfr_t MPt62;
mpfr_init2(MPt62,1024);
mpfr_add(MPt62,MPt42,MPt61,MPFR_RNDN);
mpfr_t MPt63;
mpfr_init2(MPt63,1024);
mpfr_add(MPt63,MPt44,MPt62,MPFR_RNDN);
mpfr_t MPt64;
mpfr_init2(MPt64,1024);
mpfr_add(MPt64,MPt47,MPt63,MPFR_RNDN);
mpfr_t MPt65;
mpfr_init2(MPt65,1024);
mpfr_add(MPt65,MPt49,MPt64,MPFR_RNDN);
mpfr_t MPt66;
mpfr_init2(MPt66,1024);
mpfr_add(MPt66,MPt51,MPt65,MPFR_RNDN);
mpfr_clear(MPx0);
mpfr_clear(MPx1);
mpfr_clear(MPx2);
mpfr_clear(MPx3);
mpfr_clear(MPx4);
mpfr_clear(MPx5);
mpfr_clear(MPx6);
mpfr_clear(MPx7);
mpfr_clear(MPx8);
mpfr_clear(MPx9);
mpfr_clear(MPy0);
mpfr_clear(MPy1);
mpfr_clear(MPy2);
mpfr_clear(MPy3);
mpfr_clear(MPy4);
mpfr_clear(MPy5);
mpfr_clear(MPy6);
mpfr_clear(MPy7);
mpfr_clear(MPy8);
mpfr_clear(MPy9);
mpfr_clear(MPt2);
mpfr_clear(MPt3);
mpfr_clear(MPt4);
mpfr_clear(MPt5);
mpfr_clear(MPt6);
mpfr_clear(MPt7);
mpfr_clear(MPt8);
mpfr_clear(MPt9);
mpfr_clear(MPt10);
mpfr_clear(MPt11);
mpfr_clear(MPt12);
mpfr_clear(MPt13);
mpfr_clear(MPt14);
mpfr_clear(MPt15);
mpfr_clear(MPt16);
mpfr_clear(MPt17);
mpfr_clear(MPt18);
mpfr_clear(MPt19);
mpfr_clear(MPt20);
mpfr_clear(MPt21);
mpfr_clear(MPt22);
mpfr_clear(MPt23);
mpfr_clear(MPt24);
mpfr_clear(MPt25);
mpfr_clear(MPt26);
mpfr_clear(MPt27);
mpfr_clear(MPt28);
mpfr_clear(MPt29);
mpfr_clear(MPt30);
mpfr_clear(MPt31);
mpfr_clear(MPt32);
mpfr_clear(MPt33);
mpfr_clear(MPt34);
mpfr_clear(MPt35);
mpfr_clear(MPt36);
mpfr_clear(MPt37);
mpfr_clear(MPt38);
mpfr_clear(MPt39);
mpfr_clear(MPt40);
mpfr_clear(MPt41);
mpfr_clear(MPt42);
mpfr_clear(MPt43);
mpfr_clear(MPt44);
mpfr_clear(MPt45);
mpfr_clear(MPt46);
mpfr_clear(MPt47);
mpfr_clear(MPt48);
mpfr_clear(MPt49);
mpfr_clear(MPt50);
mpfr_clear(MPt51);
mpfr_clear(MPt52);
mpfr_clear(MPt53);
mpfr_clear(MPt54);
mpfr_clear(MPt55);
mpfr_clear(MPt56);
mpfr_clear(MPt57);
mpfr_clear(MPt58);
mpfr_clear(MPt59);
mpfr_clear(MPt60);
mpfr_clear(MPt61);
mpfr_clear(MPt62);
mpfr_clear(MPt63);
mpfr_clear(MPt64);
mpfr_clear(MPt65);
mpfr_clear(MPt66);
}
//FPSyn Adapt Function
void FPSyn_Adapt_ConvexHullArea(REAL x0,REAL x1,REAL x2,REAL x3,REAL x4,REAL x5,REAL x6,REAL x7,REAL x8,REAL x9,REAL y0,REAL y1,REAL y2,REAL y3,REAL y4,REAL y5,REAL y6,REAL y7,REAL y8,REAL y9){
REAL res, err;
FPSyn_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9,&res,&err);
if ((res<err)&&(res>-err)) Exact_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
}
//IA Adapt Function
void IA_Adapt_ConvexHullArea(REAL x0,REAL x1,REAL x2,REAL x3,REAL x4,REAL x5,REAL x6,REAL x7,REAL x8,REAL x9,REAL y0,REAL y1,REAL y2,REAL y3,REAL y4,REAL y5,REAL y6,REAL y7,REAL y8,REAL y9){
REAL res, lerr, herr;
IA_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9,&res,&herr,&lerr);
if ((0<herr)&&(0>lerr)) Exact_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
}
//Main test Function
int main(){
int sample=10000;
int numtest=100;
FILE *save;
save=fopen("RuntestSave.txt","w");
FILE *report;
report=fopen("RuntestReport.txt","w");
int range = 126;
REAL varlist;
REAL res, err,herr,lerr;
REAL FPStotaltime=0; 
REAL IAtotaltime=0;
REAL FPStime=0; 
REAL IAtime=0;
clock_t before, after;
REAL x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9;
int i,j;
for(j=0; j<sample; j++){
FPStime=0;
IAtime=0;


x0=randfl(range);
x1=randfl(range);
x2=randfl(range);
x3=randfl(range);
x4=randfl(range);
x5=randfl(range);
x6=randfl(range);
x7=randfl(range);
x8=randfl(range);
x9=randfl(range);
y0=randfl(range);
y1=randfl(range);
y2=randfl(range);
y3=randfl(range);
y4=randfl(range);
y5=randfl(range);
y6=randfl(range);
y7=randfl(range);
y8=randfl(range);
y9=randfl(range);

before = clock();
FPSyn_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
FPSyn_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
FPSyn_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
FPSyn_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
FPSyn_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
FPSyn_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
FPSyn_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
FPSyn_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
FPSyn_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
FPSyn_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
FPSyn_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
FPSyn_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
FPSyn_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
FPSyn_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
FPSyn_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
FPSyn_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
FPSyn_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
FPSyn_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
FPSyn_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
FPSyn_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
FPSyn_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
FPSyn_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
FPSyn_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
FPSyn_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
FPSyn_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
FPSyn_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
FPSyn_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
FPSyn_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
FPSyn_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
FPSyn_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
FPSyn_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
FPSyn_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
FPSyn_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
FPSyn_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
FPSyn_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
FPSyn_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
FPSyn_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
FPSyn_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
FPSyn_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
FPSyn_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
FPSyn_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
FPSyn_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
FPSyn_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
FPSyn_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
FPSyn_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
FPSyn_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
FPSyn_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
FPSyn_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
FPSyn_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
FPSyn_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
FPSyn_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
FPSyn_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
FPSyn_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
FPSyn_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
FPSyn_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
FPSyn_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
FPSyn_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
FPSyn_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
FPSyn_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
FPSyn_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
FPSyn_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
FPSyn_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
FPSyn_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
FPSyn_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
FPSyn_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
FPSyn_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
FPSyn_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
FPSyn_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
FPSyn_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
FPSyn_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
FPSyn_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
FPSyn_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
FPSyn_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
FPSyn_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
FPSyn_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
FPSyn_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
FPSyn_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
FPSyn_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
FPSyn_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
FPSyn_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
FPSyn_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
FPSyn_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
FPSyn_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
FPSyn_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
FPSyn_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
FPSyn_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
FPSyn_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
FPSyn_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
FPSyn_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
FPSyn_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
FPSyn_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
FPSyn_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
FPSyn_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
FPSyn_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
FPSyn_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
FPSyn_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
FPSyn_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
FPSyn_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
FPSyn_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
FPSyn_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
FPSyn_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
FPSyn_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
FPSyn_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
FPSyn_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
FPSyn_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
FPSyn_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
FPSyn_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
FPSyn_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
FPSyn_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
FPSyn_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
FPSyn_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
FPSyn_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
FPSyn_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
FPSyn_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
FPSyn_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
FPSyn_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
FPSyn_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
FPSyn_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
FPSyn_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
FPSyn_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
FPSyn_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
FPSyn_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
FPSyn_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
FPSyn_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
FPSyn_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
FPSyn_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
FPSyn_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
FPSyn_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
FPSyn_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
FPSyn_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
FPSyn_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
FPSyn_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
FPSyn_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
FPSyn_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
FPSyn_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
FPSyn_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
FPSyn_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
FPSyn_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
FPSyn_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
FPSyn_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
FPSyn_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
FPSyn_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
FPSyn_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
FPSyn_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
FPSyn_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
FPSyn_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
FPSyn_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
FPSyn_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
FPSyn_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
FPSyn_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
FPSyn_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
FPSyn_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
FPSyn_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
FPSyn_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
FPSyn_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
FPSyn_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
FPSyn_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
FPSyn_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
FPSyn_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
FPSyn_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
FPSyn_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
FPSyn_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
FPSyn_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
FPSyn_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
FPSyn_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
FPSyn_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
FPSyn_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
FPSyn_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
FPSyn_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
FPSyn_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
FPSyn_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
FPSyn_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
FPSyn_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
FPSyn_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
FPSyn_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
FPSyn_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
FPSyn_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
FPSyn_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
FPSyn_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
FPSyn_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
FPSyn_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
FPSyn_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
FPSyn_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
FPSyn_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
FPSyn_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
FPSyn_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
FPSyn_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
FPSyn_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
FPSyn_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
FPSyn_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
FPSyn_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
FPSyn_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
FPSyn_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
FPSyn_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
FPSyn_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
FPSyn_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
FPSyn_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
FPSyn_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
FPSyn_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
FPSyn_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
FPSyn_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
FPSyn_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
FPSyn_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
FPSyn_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
FPSyn_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
FPSyn_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
FPSyn_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
FPSyn_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
FPSyn_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
FPSyn_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
FPSyn_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
FPSyn_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
FPSyn_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
FPSyn_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
FPSyn_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
FPSyn_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
FPSyn_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
FPSyn_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
FPSyn_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
FPSyn_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
FPSyn_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
FPSyn_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
FPSyn_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
FPSyn_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
FPSyn_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
FPSyn_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
FPSyn_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
FPSyn_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
FPSyn_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
FPSyn_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
FPSyn_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
FPSyn_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
FPSyn_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
FPSyn_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
FPSyn_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
FPSyn_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
FPSyn_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
FPSyn_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
FPSyn_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
FPSyn_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
FPSyn_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
FPSyn_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
FPSyn_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
FPSyn_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
FPSyn_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
FPSyn_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
FPSyn_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
FPSyn_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
FPSyn_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
FPSyn_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
FPSyn_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
FPSyn_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
FPSyn_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
FPSyn_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
FPSyn_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
FPSyn_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
FPSyn_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
FPSyn_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
FPSyn_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
FPSyn_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
FPSyn_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
FPSyn_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
FPSyn_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
FPSyn_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
FPSyn_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
FPSyn_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
FPSyn_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
FPSyn_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
FPSyn_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
FPSyn_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
FPSyn_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
FPSyn_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
FPSyn_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
FPSyn_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
FPSyn_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
FPSyn_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
FPSyn_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
FPSyn_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
FPSyn_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
FPSyn_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
FPSyn_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
FPSyn_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
FPSyn_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
FPSyn_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
FPSyn_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
FPSyn_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
FPSyn_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
FPSyn_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
FPSyn_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
FPSyn_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
FPSyn_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
FPSyn_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
FPSyn_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
FPSyn_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
FPSyn_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
FPSyn_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
FPSyn_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
FPSyn_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
FPSyn_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
FPSyn_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
FPSyn_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
FPSyn_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
FPSyn_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
FPSyn_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
FPSyn_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
FPSyn_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
FPSyn_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
FPSyn_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
FPSyn_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
FPSyn_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
FPSyn_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
FPSyn_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
FPSyn_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
FPSyn_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
FPSyn_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
FPSyn_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
FPSyn_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
FPSyn_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
FPSyn_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
FPSyn_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
FPSyn_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
FPSyn_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
FPSyn_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
FPSyn_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
FPSyn_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
FPSyn_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
FPSyn_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
FPSyn_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
FPSyn_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
FPSyn_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
FPSyn_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
FPSyn_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
FPSyn_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
FPSyn_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
FPSyn_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
FPSyn_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
FPSyn_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
FPSyn_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
FPSyn_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
FPSyn_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
FPSyn_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
FPSyn_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
FPSyn_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
FPSyn_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
FPSyn_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
FPSyn_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
FPSyn_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
FPSyn_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
FPSyn_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
FPSyn_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
FPSyn_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
FPSyn_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
FPSyn_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
FPSyn_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
FPSyn_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
FPSyn_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
FPSyn_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
FPSyn_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
FPSyn_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
FPSyn_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
FPSyn_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
FPSyn_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
FPSyn_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
FPSyn_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
FPSyn_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
FPSyn_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
FPSyn_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
FPSyn_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
FPSyn_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
FPSyn_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
FPSyn_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
FPSyn_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
FPSyn_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
FPSyn_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
FPSyn_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
FPSyn_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
FPSyn_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
FPSyn_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
FPSyn_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
FPSyn_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
FPSyn_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
FPSyn_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
FPSyn_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
FPSyn_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
FPSyn_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
FPSyn_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
FPSyn_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
FPSyn_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
FPSyn_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
FPSyn_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
FPSyn_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
FPSyn_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
FPSyn_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
FPSyn_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
FPSyn_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
FPSyn_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
FPSyn_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
FPSyn_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
FPSyn_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
FPSyn_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
FPSyn_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
FPSyn_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
FPSyn_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
FPSyn_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
FPSyn_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
FPSyn_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
FPSyn_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
FPSyn_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
FPSyn_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
FPSyn_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
FPSyn_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
FPSyn_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
FPSyn_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
FPSyn_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
FPSyn_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
FPSyn_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
FPSyn_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
FPSyn_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
FPSyn_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
FPSyn_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
FPSyn_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
FPSyn_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
FPSyn_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
FPSyn_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
FPSyn_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
FPSyn_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
FPSyn_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
FPSyn_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
FPSyn_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
FPSyn_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
FPSyn_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
FPSyn_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
FPSyn_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
FPSyn_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
FPSyn_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
FPSyn_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
FPSyn_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
FPSyn_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
FPSyn_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
FPSyn_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
FPSyn_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
FPSyn_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
FPSyn_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
FPSyn_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
FPSyn_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
FPSyn_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
FPSyn_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
FPSyn_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
FPSyn_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
FPSyn_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
FPSyn_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
FPSyn_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
FPSyn_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
FPSyn_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
FPSyn_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
FPSyn_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
FPSyn_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
FPSyn_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
FPSyn_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
FPSyn_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
FPSyn_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
FPSyn_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
FPSyn_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
FPSyn_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
FPSyn_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
FPSyn_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
FPSyn_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
FPSyn_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
FPSyn_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
FPSyn_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
FPSyn_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
FPSyn_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
FPSyn_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
FPSyn_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
FPSyn_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
FPSyn_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
FPSyn_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
FPSyn_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
FPSyn_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
FPSyn_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
FPSyn_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
FPSyn_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
FPSyn_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
FPSyn_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
FPSyn_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
FPSyn_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
FPSyn_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
FPSyn_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
FPSyn_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
FPSyn_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
FPSyn_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
FPSyn_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
FPSyn_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
FPSyn_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
FPSyn_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
FPSyn_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
FPSyn_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
FPSyn_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
FPSyn_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
FPSyn_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
FPSyn_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
FPSyn_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
FPSyn_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
FPSyn_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
FPSyn_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
FPSyn_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
FPSyn_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
FPSyn_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
FPSyn_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
FPSyn_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
FPSyn_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
FPSyn_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
FPSyn_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
FPSyn_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
FPSyn_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
FPSyn_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
FPSyn_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
FPSyn_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
FPSyn_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
FPSyn_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
FPSyn_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
FPSyn_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
FPSyn_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
FPSyn_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
FPSyn_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
FPSyn_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
FPSyn_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
FPSyn_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
FPSyn_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
FPSyn_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
FPSyn_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
FPSyn_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
FPSyn_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
FPSyn_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
FPSyn_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
FPSyn_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
FPSyn_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
FPSyn_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
FPSyn_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
FPSyn_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
FPSyn_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
FPSyn_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
FPSyn_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
FPSyn_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
FPSyn_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
FPSyn_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
FPSyn_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
FPSyn_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
FPSyn_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
FPSyn_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
FPSyn_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
FPSyn_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
FPSyn_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
FPSyn_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
FPSyn_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
FPSyn_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
FPSyn_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
FPSyn_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
FPSyn_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
FPSyn_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
FPSyn_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
FPSyn_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
FPSyn_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
FPSyn_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
FPSyn_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
FPSyn_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
FPSyn_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
FPSyn_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
FPSyn_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
FPSyn_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
FPSyn_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
FPSyn_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
FPSyn_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
FPSyn_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
FPSyn_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
FPSyn_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
FPSyn_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
FPSyn_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
FPSyn_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
FPSyn_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
FPSyn_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
FPSyn_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
FPSyn_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
FPSyn_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
FPSyn_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
FPSyn_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
FPSyn_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
FPSyn_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
FPSyn_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
FPSyn_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
FPSyn_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
FPSyn_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
FPSyn_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
FPSyn_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
FPSyn_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
FPSyn_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
FPSyn_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
FPSyn_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
FPSyn_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
FPSyn_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
FPSyn_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
FPSyn_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
FPSyn_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
FPSyn_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
FPSyn_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
FPSyn_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
FPSyn_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
FPSyn_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
FPSyn_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
FPSyn_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
FPSyn_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
FPSyn_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
FPSyn_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
FPSyn_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
FPSyn_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
FPSyn_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
FPSyn_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
FPSyn_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
FPSyn_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
FPSyn_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
FPSyn_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
FPSyn_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
FPSyn_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
FPSyn_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
FPSyn_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
FPSyn_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
FPSyn_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
FPSyn_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
FPSyn_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
FPSyn_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
FPSyn_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
FPSyn_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
FPSyn_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
FPSyn_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
FPSyn_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
FPSyn_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
FPSyn_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
FPSyn_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
FPSyn_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
FPSyn_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
FPSyn_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
FPSyn_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
FPSyn_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
FPSyn_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
FPSyn_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
FPSyn_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
FPSyn_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
FPSyn_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
FPSyn_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
FPSyn_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
FPSyn_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
FPSyn_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
FPSyn_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
FPSyn_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
FPSyn_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
FPSyn_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
FPSyn_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
FPSyn_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
FPSyn_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
FPSyn_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
FPSyn_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
FPSyn_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
FPSyn_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
FPSyn_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
FPSyn_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
FPSyn_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
FPSyn_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
FPSyn_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
FPSyn_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
FPSyn_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
FPSyn_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
FPSyn_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
FPSyn_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
FPSyn_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
FPSyn_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
FPSyn_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
FPSyn_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
FPSyn_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
FPSyn_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
FPSyn_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
FPSyn_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
FPSyn_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
FPSyn_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
FPSyn_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
FPSyn_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
FPSyn_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
FPSyn_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
FPSyn_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
FPSyn_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
FPSyn_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
FPSyn_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
FPSyn_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
FPSyn_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
FPSyn_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
FPSyn_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
FPSyn_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
FPSyn_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
FPSyn_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
FPSyn_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
FPSyn_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
FPSyn_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
FPSyn_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
FPSyn_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
FPSyn_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
FPSyn_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
FPSyn_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
FPSyn_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
FPSyn_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
FPSyn_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
FPSyn_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
FPSyn_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
FPSyn_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
FPSyn_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
FPSyn_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
FPSyn_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
FPSyn_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
FPSyn_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
FPSyn_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
FPSyn_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
FPSyn_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
FPSyn_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
FPSyn_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
FPSyn_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
FPSyn_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
FPSyn_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
FPSyn_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
FPSyn_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
FPSyn_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
FPSyn_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
FPSyn_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
FPSyn_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
FPSyn_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
FPSyn_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
FPSyn_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
FPSyn_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
FPSyn_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
FPSyn_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
FPSyn_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
FPSyn_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
FPSyn_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
FPSyn_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
FPSyn_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
FPSyn_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
FPSyn_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
FPSyn_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
FPSyn_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
FPSyn_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
FPSyn_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
FPSyn_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
FPSyn_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
FPSyn_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
FPSyn_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
FPSyn_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
FPSyn_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
FPSyn_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
FPSyn_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
FPSyn_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
FPSyn_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
FPSyn_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
FPSyn_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
FPSyn_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
FPSyn_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
FPSyn_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
FPSyn_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
FPSyn_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
FPSyn_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
FPSyn_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
FPSyn_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
FPSyn_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
FPSyn_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
FPSyn_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
FPSyn_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
FPSyn_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
FPSyn_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
FPSyn_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
FPSyn_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
FPSyn_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
FPSyn_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
FPSyn_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
FPSyn_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
FPSyn_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
FPSyn_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
FPSyn_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
FPSyn_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
FPSyn_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
FPSyn_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
FPSyn_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
FPSyn_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
FPSyn_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
FPSyn_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
FPSyn_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
FPSyn_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
FPSyn_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
FPSyn_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
FPSyn_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
FPSyn_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
FPSyn_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
FPSyn_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
FPSyn_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
FPSyn_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
FPSyn_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
FPSyn_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
FPSyn_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
FPSyn_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
FPSyn_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
FPSyn_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
FPSyn_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
FPSyn_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
FPSyn_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
FPSyn_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
FPSyn_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
FPSyn_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
FPSyn_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
FPSyn_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
FPSyn_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
FPSyn_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
FPSyn_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
FPSyn_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
FPSyn_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
FPSyn_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
FPSyn_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
FPSyn_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
FPSyn_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
FPSyn_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
FPSyn_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
FPSyn_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
FPSyn_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
FPSyn_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
FPSyn_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
FPSyn_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
FPSyn_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
FPSyn_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
FPSyn_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
FPSyn_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
FPSyn_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
FPSyn_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
FPSyn_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
FPSyn_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
FPSyn_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
FPSyn_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
FPSyn_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
FPSyn_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
FPSyn_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
FPSyn_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
FPSyn_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
FPSyn_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
FPSyn_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
FPSyn_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
FPSyn_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
FPSyn_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
FPSyn_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
FPSyn_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
FPSyn_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
FPSyn_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
FPSyn_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
FPSyn_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
FPSyn_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
FPSyn_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
FPSyn_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
FPSyn_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
FPSyn_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
FPSyn_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
FPSyn_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
FPSyn_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
FPSyn_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
FPSyn_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
FPSyn_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
FPSyn_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
FPSyn_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
FPSyn_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
FPSyn_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
FPSyn_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
FPSyn_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
FPSyn_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
FPSyn_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
FPSyn_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
FPSyn_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
FPSyn_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
FPSyn_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
FPSyn_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
FPSyn_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
FPSyn_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
FPSyn_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
FPSyn_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
FPSyn_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
FPSyn_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
FPSyn_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
FPSyn_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
FPSyn_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
FPSyn_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
FPSyn_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
FPSyn_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
FPSyn_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
FPSyn_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
FPSyn_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
FPSyn_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
FPSyn_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
FPSyn_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
FPSyn_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
FPSyn_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
FPSyn_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
FPSyn_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
FPSyn_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
FPSyn_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
FPSyn_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
FPSyn_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
FPSyn_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
FPSyn_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
FPSyn_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
FPSyn_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
FPSyn_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
FPSyn_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
FPSyn_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
FPSyn_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
FPSyn_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
FPSyn_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
FPSyn_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
FPSyn_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
FPSyn_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
FPSyn_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
FPSyn_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
FPSyn_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
FPSyn_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
FPSyn_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
FPSyn_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
FPSyn_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
FPSyn_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
FPSyn_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
FPSyn_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
FPSyn_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
FPSyn_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
FPSyn_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
FPSyn_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
FPSyn_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
FPSyn_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
FPSyn_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
FPSyn_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
FPSyn_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
FPSyn_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
FPSyn_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
FPSyn_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
FPSyn_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
FPSyn_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
FPSyn_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
FPSyn_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
FPSyn_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
FPSyn_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
FPSyn_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
FPSyn_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
FPSyn_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
FPSyn_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
FPSyn_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
FPSyn_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
FPSyn_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
FPSyn_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
FPSyn_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
FPSyn_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
FPSyn_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
FPSyn_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
FPSyn_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
FPSyn_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
FPSyn_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
FPSyn_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
FPSyn_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
FPSyn_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
FPSyn_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
FPSyn_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
FPSyn_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
FPSyn_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
FPSyn_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
FPSyn_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
FPSyn_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
FPSyn_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
FPSyn_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
FPSyn_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
FPSyn_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
FPSyn_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
FPSyn_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
FPSyn_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
FPSyn_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
FPSyn_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
FPSyn_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
FPSyn_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
FPSyn_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
FPSyn_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
FPSyn_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
FPSyn_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
FPSyn_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
FPSyn_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
FPSyn_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
FPSyn_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
FPSyn_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
FPSyn_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
FPSyn_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
FPSyn_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
FPSyn_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
FPSyn_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
FPSyn_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
FPSyn_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
FPSyn_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
after = clock();
FPStime = ((double)(after - before))/CLOCKS_PER_SEC;
FPStotaltime = FPStotaltime + FPStime;


before = clock();
IA_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
IA_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
IA_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
IA_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
IA_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
IA_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
IA_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
IA_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
IA_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
IA_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
IA_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
IA_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
IA_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
IA_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
IA_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
IA_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
IA_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
IA_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
IA_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
IA_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
IA_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
IA_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
IA_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
IA_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
IA_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
IA_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
IA_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
IA_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
IA_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
IA_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
IA_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
IA_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
IA_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
IA_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
IA_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
IA_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
IA_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
IA_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
IA_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
IA_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
IA_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
IA_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
IA_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
IA_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
IA_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
IA_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
IA_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
IA_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
IA_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
IA_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
IA_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
IA_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
IA_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
IA_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
IA_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
IA_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
IA_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
IA_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
IA_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
IA_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
IA_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
IA_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
IA_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
IA_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
IA_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
IA_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
IA_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
IA_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
IA_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
IA_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
IA_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
IA_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
IA_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
IA_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
IA_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
IA_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
IA_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
IA_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
IA_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
IA_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
IA_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
IA_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
IA_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
IA_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
IA_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
IA_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
IA_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
IA_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
IA_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
IA_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
IA_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
IA_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
IA_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
IA_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
IA_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
IA_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
IA_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
IA_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
IA_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
IA_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
IA_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
IA_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
IA_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
IA_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
IA_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
IA_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
IA_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
IA_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
IA_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
IA_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
IA_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
IA_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
IA_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
IA_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
IA_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
IA_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
IA_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
IA_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
IA_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
IA_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
IA_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
IA_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
IA_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
IA_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
IA_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
IA_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
IA_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
IA_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
IA_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
IA_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
IA_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
IA_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
IA_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
IA_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
IA_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
IA_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
IA_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
IA_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
IA_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
IA_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
IA_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
IA_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
IA_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
IA_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
IA_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
IA_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
IA_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
IA_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
IA_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
IA_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
IA_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
IA_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
IA_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
IA_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
IA_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
IA_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
IA_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
IA_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
IA_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
IA_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
IA_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
IA_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
IA_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
IA_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
IA_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
IA_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
IA_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
IA_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
IA_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
IA_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
IA_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
IA_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
IA_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
IA_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
IA_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
IA_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
IA_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
IA_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
IA_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
IA_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
IA_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
IA_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
IA_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
IA_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
IA_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
IA_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
IA_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
IA_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
IA_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
IA_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
IA_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
IA_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
IA_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
IA_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
IA_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
IA_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
IA_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
IA_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
IA_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
IA_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
IA_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
IA_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
IA_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
IA_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
IA_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
IA_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
IA_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
IA_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
IA_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
IA_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
IA_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
IA_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
IA_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
IA_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
IA_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
IA_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
IA_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
IA_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
IA_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
IA_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
IA_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
IA_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
IA_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
IA_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
IA_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
IA_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
IA_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
IA_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
IA_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
IA_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
IA_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
IA_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
IA_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
IA_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
IA_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
IA_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
IA_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
IA_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
IA_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
IA_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
IA_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
IA_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
IA_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
IA_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
IA_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
IA_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
IA_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
IA_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
IA_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
IA_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
IA_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
IA_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
IA_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
IA_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
IA_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
IA_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
IA_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
IA_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
IA_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
IA_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
IA_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
IA_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
IA_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
IA_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
IA_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
IA_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
IA_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
IA_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
IA_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
IA_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
IA_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
IA_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
IA_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
IA_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
IA_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
IA_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
IA_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
IA_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
IA_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
IA_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
IA_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
IA_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
IA_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
IA_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
IA_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
IA_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
IA_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
IA_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
IA_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
IA_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
IA_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
IA_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
IA_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
IA_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
IA_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
IA_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
IA_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
IA_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
IA_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
IA_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
IA_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
IA_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
IA_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
IA_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
IA_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
IA_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
IA_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
IA_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
IA_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
IA_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
IA_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
IA_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
IA_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
IA_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
IA_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
IA_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
IA_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
IA_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
IA_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
IA_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
IA_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
IA_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
IA_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
IA_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
IA_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
IA_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
IA_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
IA_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
IA_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
IA_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
IA_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
IA_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
IA_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
IA_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
IA_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
IA_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
IA_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
IA_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
IA_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
IA_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
IA_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
IA_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
IA_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
IA_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
IA_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
IA_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
IA_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
IA_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
IA_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
IA_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
IA_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
IA_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
IA_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
IA_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
IA_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
IA_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
IA_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
IA_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
IA_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
IA_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
IA_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
IA_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
IA_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
IA_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
IA_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
IA_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
IA_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
IA_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
IA_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
IA_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
IA_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
IA_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
IA_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
IA_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
IA_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
IA_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
IA_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
IA_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
IA_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
IA_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
IA_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
IA_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
IA_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
IA_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
IA_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
IA_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
IA_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
IA_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
IA_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
IA_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
IA_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
IA_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
IA_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
IA_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
IA_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
IA_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
IA_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
IA_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
IA_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
IA_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
IA_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
IA_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
IA_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
IA_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
IA_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
IA_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
IA_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
IA_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
IA_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
IA_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
IA_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
IA_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
IA_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
IA_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
IA_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
IA_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
IA_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
IA_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
IA_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
IA_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
IA_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
IA_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
IA_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
IA_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
IA_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
IA_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
IA_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
IA_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
IA_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
IA_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
IA_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
IA_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
IA_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
IA_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
IA_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
IA_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
IA_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
IA_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
IA_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
IA_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
IA_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
IA_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
IA_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
IA_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
IA_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
IA_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
IA_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
IA_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
IA_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
IA_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
IA_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
IA_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
IA_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
IA_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
IA_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
IA_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
IA_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
IA_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
IA_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
IA_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
IA_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
IA_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
IA_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
IA_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
IA_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
IA_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
IA_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
IA_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
IA_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
IA_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
IA_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
IA_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
IA_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
IA_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
IA_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
IA_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
IA_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
IA_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
IA_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
IA_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
IA_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
IA_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
IA_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
IA_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
IA_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
IA_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
IA_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
IA_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
IA_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
IA_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
IA_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
IA_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
IA_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
IA_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
IA_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
IA_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
IA_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
IA_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
IA_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
IA_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
IA_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
IA_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
IA_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
IA_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
IA_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
IA_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
IA_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
IA_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
IA_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
IA_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
IA_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
IA_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
IA_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
IA_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
IA_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
IA_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
IA_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
IA_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
IA_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
IA_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
IA_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
IA_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
IA_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
IA_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
IA_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
IA_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
IA_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
IA_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
IA_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
IA_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
IA_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
IA_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
IA_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
IA_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
IA_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
IA_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
IA_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
IA_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
IA_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
IA_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
IA_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
IA_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
IA_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
IA_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
IA_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
IA_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
IA_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
IA_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
IA_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
IA_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
IA_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
IA_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
IA_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
IA_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
IA_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
IA_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
IA_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
IA_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
IA_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
IA_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
IA_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
IA_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
IA_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
IA_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
IA_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
IA_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
IA_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
IA_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
IA_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
IA_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
IA_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
IA_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
IA_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
IA_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
IA_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
IA_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
IA_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
IA_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
IA_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
IA_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
IA_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
IA_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
IA_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
IA_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
IA_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
IA_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
IA_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
IA_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
IA_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
IA_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
IA_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
IA_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
IA_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
IA_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
IA_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
IA_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
IA_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
IA_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
IA_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
IA_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
IA_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
IA_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
IA_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
IA_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
IA_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
IA_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
IA_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
IA_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
IA_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
IA_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
IA_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
IA_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
IA_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
IA_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
IA_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
IA_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
IA_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
IA_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
IA_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
IA_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
IA_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
IA_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
IA_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
IA_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
IA_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
IA_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
IA_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
IA_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
IA_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
IA_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
IA_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
IA_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
IA_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
IA_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
IA_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
IA_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
IA_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
IA_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
IA_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
IA_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
IA_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
IA_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
IA_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
IA_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
IA_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
IA_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
IA_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
IA_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
IA_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
IA_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
IA_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
IA_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
IA_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
IA_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
IA_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
IA_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
IA_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
IA_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
IA_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
IA_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
IA_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
IA_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
IA_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
IA_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
IA_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
IA_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
IA_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
IA_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
IA_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
IA_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
IA_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
IA_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
IA_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
IA_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
IA_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
IA_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
IA_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
IA_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
IA_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
IA_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
IA_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
IA_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
IA_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
IA_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
IA_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
IA_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
IA_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
IA_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
IA_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
IA_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
IA_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
IA_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
IA_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
IA_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
IA_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
IA_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
IA_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
IA_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
IA_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
IA_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
IA_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
IA_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
IA_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
IA_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
IA_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
IA_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
IA_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
IA_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
IA_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
IA_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
IA_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
IA_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
IA_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
IA_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
IA_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
IA_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
IA_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
IA_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
IA_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
IA_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
IA_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
IA_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
IA_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
IA_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
IA_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
IA_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
IA_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
IA_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
IA_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
IA_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
IA_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
IA_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
IA_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
IA_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
IA_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
IA_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
IA_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
IA_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
IA_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
IA_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
IA_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
IA_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
IA_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
IA_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
IA_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
IA_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
IA_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
IA_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
IA_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
IA_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
IA_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
IA_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
IA_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
IA_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
IA_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
IA_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
IA_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
IA_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
IA_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
IA_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
IA_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
IA_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
IA_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
IA_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
IA_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
IA_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
IA_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
IA_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
IA_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
IA_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
IA_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
IA_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
IA_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
IA_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
IA_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
IA_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
IA_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
IA_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
IA_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
IA_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
IA_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
IA_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
IA_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
IA_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
IA_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
IA_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
IA_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
IA_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
IA_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
IA_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
IA_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
IA_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
IA_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
IA_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
IA_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
IA_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
IA_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
IA_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
IA_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
IA_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
IA_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
IA_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
IA_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
IA_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
IA_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
IA_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
IA_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
IA_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
IA_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
IA_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
IA_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
IA_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
IA_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
IA_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
IA_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
IA_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
IA_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
IA_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
IA_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
IA_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
IA_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
IA_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
IA_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
IA_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
IA_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
IA_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
IA_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
IA_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
IA_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
IA_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
IA_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
IA_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
IA_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
IA_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
IA_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
IA_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
IA_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
IA_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
IA_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
IA_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
IA_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
IA_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
IA_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
IA_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
IA_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
IA_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
IA_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
IA_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
IA_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
IA_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
IA_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
IA_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
IA_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
IA_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
IA_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
IA_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
IA_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
IA_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
IA_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
IA_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
IA_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
IA_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
IA_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
IA_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
IA_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
IA_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
IA_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
IA_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
IA_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
IA_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
IA_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
IA_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
IA_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
IA_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
IA_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
IA_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
IA_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
IA_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
IA_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
IA_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
IA_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
IA_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
IA_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
IA_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
IA_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
IA_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
IA_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
IA_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
IA_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
IA_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
IA_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
IA_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
IA_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
IA_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
IA_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
IA_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
IA_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
IA_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
IA_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
IA_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
IA_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
IA_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
IA_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
IA_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
IA_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
IA_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
IA_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
IA_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
IA_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
IA_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
IA_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
IA_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
IA_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
IA_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
IA_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
IA_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
IA_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
IA_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
IA_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
IA_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
IA_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
IA_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
IA_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
IA_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
IA_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
IA_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
IA_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
IA_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
IA_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
IA_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
IA_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
IA_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
IA_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
IA_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
IA_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
IA_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
IA_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
IA_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
IA_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
IA_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
IA_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
IA_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
IA_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
IA_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
IA_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
IA_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
IA_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
IA_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
IA_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
IA_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
IA_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
IA_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
IA_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
IA_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
IA_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
IA_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
IA_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
IA_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
IA_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
IA_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
IA_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
IA_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
IA_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
IA_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
IA_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
IA_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
IA_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
IA_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
IA_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
IA_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
IA_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
IA_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
IA_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
IA_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
IA_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
IA_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
IA_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
IA_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
IA_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
IA_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
IA_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
IA_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
IA_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
IA_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
IA_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
IA_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
IA_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
IA_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
IA_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
IA_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
IA_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
IA_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
IA_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
IA_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
IA_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
IA_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
IA_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
IA_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
IA_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
IA_Adapt_ConvexHullArea(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9);
after = clock();
IAtime = ((double)(after - before))/CLOCKS_PER_SEC;
IAtotaltime = IAtotaltime + IAtime;


printf("Total time FPS: %e \n", FPStime);
printf("Total time IA: %e \n\n", IAtime);

fprintf(save,"%e\n",FPStime);
fprintf(save,"%e\n",IAtime);
}

printf("Average FPS time : %f\n",FPStotaltime/sample);
printf("Average IA time :%f\n",IAtotaltime/sample);


fprintf(report,"Average FPS time : %f\n",FPStotaltime/sample);
fprintf(report,"Average IA time : %f\n",IAtotaltime/sample);

fclose(save);
fclose(report);
}