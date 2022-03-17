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
void FPSyn_Orient3D(REAL ax,REAL ay,REAL az,REAL bx,REAL by,REAL bz,REAL cx,REAL cy,REAL cz,REAL dx,REAL dy,REAL dz,REAL* result, REAL* error_bound){
//Calculation steps for input: ;
REAL t2=by-dy;
REAL t4=cz-dz;
REAL t5=bz-dz;
REAL t6=cy-dy;
REAL t7=ay-dy;
REAL t8=az-dz;
REAL t9=ax-dx;
REAL t10=t2*t4;
REAL t11=t5*t6;
REAL t13=t10-t11;
REAL t14=t13*t9;
REAL t15=cx-dx;
REAL t16=t5*t7;
REAL t17=t2*t8;
REAL t19=t16-t17;
REAL t20=t15*t19;
REAL t21=bx-dx;
REAL t22=t4*t7;
REAL t23=t6*t8;
REAL t25=t22-t23;
REAL t26=t21*t25;
REAL t28=t14+t20;
REAL t29=t28-t26;
//Calculation steps for Error-bound: ;
t15=fabs(t15);
t22=fabs(t22);
t23=fabs(t23);
t10=fabs(t10);
t11=fabs(t11);
t29=fabs(t29);
t17=fabs(t17);
t16=fabs(t16);
t21=fabs(t21);
t9=fabs(t9);
REAL tem0=t16+t17;
REAL tem1=t15*tem0;
REAL tem2=t22+t23;
REAL tem3=t21*tem2;
REAL tem4=t10+t11;
REAL tem5=t9*tem4;
REAL tem6=tem1+tem3;
REAL tem7=tem5+tem6;
REAL tem8=7.7715611723761e-16*tem7;
REAL tem9=1.11022302462516e-16*t29;
REAL tem10=tem8+tem9;
*result=t29;
*error_bound=1.0000000001*tem10;
}
//Manual Function
void Manual_Orient3D(REAL ax,REAL ay,REAL az,REAL bx,REAL by,REAL bz,REAL cx,REAL cy,REAL cz,REAL dx,REAL dy,REAL dz,REAL* result, REAL* error_bound){

  REAL adx, bdx, cdx, ady, bdy, cdy, adz, bdz, cdz;
  REAL bdxcdy, cdxbdy, cdxady, adxcdy, adxbdy, bdxady;
  REAL det;
  REAL permanent;

  adx = ax - dx;
  bdx = bx - dx;
  cdx = cx - dx;
  ady = ay - dy;
  bdy = by - dy;
  cdy = cy - dy;
  adz = az - dz;
  bdz = bz - dz;
  cdz = cz - dz;

  bdxcdy = bdx * cdy;
  cdxbdy = cdx * bdy;

  cdxady = cdx * ady;
  adxcdy = adx * cdy;

  adxbdy = adx * bdy;
  bdxady = bdx * ady;

  det = adz * (bdxcdy - cdxbdy)
      + bdz * (cdxady - adxcdy)
      + cdz * (adxbdy - bdxady);

  permanent = (Absolute(bdxcdy) + Absolute(cdxbdy)) * Absolute(adz)
            + (Absolute(cdxady) + Absolute(adxcdy)) * Absolute(bdz)
            + (Absolute(adxbdy) + Absolute(bdxady)) * Absolute(cdz);

*result=det;
*error_bound = 7.771561172376103e-16 * permanent;

}
//IA Function
void IA_Orient3D(REAL ax,REAL ay,REAL az,REAL bx,REAL by,REAL bz,REAL cx,REAL cy,REAL cz,REAL dx,REAL dy,REAL dz,REAL* result, REAL* hbound, REAL* lbound){
REAL t2=by-dy;
REAL t4=cz-dz;
REAL t5=bz-dz;
REAL t6=cy-dy;
REAL t7=ay-dy;
REAL t8=az-dz;
REAL t9=ax-dx;
REAL t10=t2*t4;
REAL t11=t5*t6;
REAL t13=t10-t11;
REAL t14=t13*t9;
REAL t15=cx-dx;
REAL t16=t5*t7;
REAL t17=t2*t8;
REAL t19=t16-t17;
REAL t20=t15*t19;
REAL t21=bx-dx;
REAL t22=t4*t7;
REAL t23=t6*t8;
REAL t25=t22-t23;
REAL t26=t21*t25;
REAL t28=t14+t20;
REAL t29=t28-t26;
//Calculate bound
REAL eps=pow(2,-53);
int inferr=0;
int roundmode = fegetround();
fesetround(FE_UPWARD);
REAL t2_e_1=fabs(t2);
REAL t2_e=eps*t2_e_1;
REAL t2_h=t2+t2_e;
REAL t2_l_n=t2_e-t2;
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
REAL  t10_h, t10_l_n;
if(t2_h<=0){
if(t4_h<=0){
REAL t2_h_n=-t2_h;
REAL t10_l_n_1=t2_h_n*t4_h;
REAL t10_h_1=t2_l_n*t4_l_n;
t10_h=t10_h_1+t10_e;
t10_l_n=t10_l_n_1+t10_e;
}
else if(t4_l_n>=0){
REAL t10_l_n_1=t2_l_n*t4_h;
REAL t10_h_1=t2_l_n*t4_l_n;
t10_h=t10_h_1+t10_e;
t10_l_n=t10_l_n_1+t10_e;
}
else{
REAL t10_l_n_1=t2_l_n*t4_h;
REAL t2_h_n=-t2_h;
REAL t10_h_1=t2_h_n*t4_l_n;
t10_h=t10_h_1+t10_e;
t10_l_n=t10_l_n_1+t10_e;
}
}
else if(t2_l_n>=0){
if(t4_h<=0){
REAL t10_l_n_1=t2_h*t4_l_n;
REAL t10_h_1=t2_l_n*t4_l_n;
t10_h=t10_h_1+t10_e;
t10_l_n=t10_l_n_1+t10_e;
}
else if(t4_l_n>=0){
REAL t10_l_n_1_1=t2_l_n*t4_h;
REAL t10_l_n_1_2=t4_l_n*t2_h;
REAL t10_l_n_1=max(t10_l_n_1_1,t10_l_n_1_2);
REAL t10_h_1_1=t2_l_n*t4_l_n;
REAL t10_h_1_2=t2_h*t4_h;
REAL t10_h_1=max(t10_h_1_1,t10_h_1_2);
t10_h=t10_h_1+t10_e;
t10_l_n=t10_l_n_1+t10_e;
}
else{
REAL t10_l_n_1=t2_l_n*t4_h;
REAL t10_h_1=t2_h*t4_h;
t10_h=t10_h_1+t10_e;
t10_l_n=t10_l_n_1+t10_e;
}
}
else{
if(t4_h<=0){
REAL t10_l_n_1=t4_l_n*t2_h;
REAL t4_h_n=-t4_h;
REAL t10_h_1=t4_h_n*t2_l_n;
t10_h=t10_h_1+t10_e;
t10_l_n=t10_l_n_1+t10_e;
}
else if(t4_l_n>=0){
REAL t10_l_n_1=t4_l_n*t2_h;
REAL t10_h_1=t2_h*t4_h;
t10_h=t10_h_1+t10_e;
t10_l_n=t10_l_n_1+t10_e;
}
else{
REAL t2_l=-t2_l_n;
REAL t10_l_n_1=t2_l*t4_l_n;
REAL t10_h_1=t2_h*t4_h;
t10_h=t10_h_1+t10_e;
t10_l_n=t10_l_n_1+t10_e;
}
}
REAL t11_e_1=fabs(t11);
REAL t11_e=eps*t11_e_1;
REAL  t11_h, t11_l_n;
if(t5_h<=0){
if(t6_h<=0){
REAL t5_h_n=-t5_h;
REAL t11_l_n_1=t5_h_n*t6_h;
REAL t11_h_1=t5_l_n*t6_l_n;
t11_h=t11_h_1+t11_e;
t11_l_n=t11_l_n_1+t11_e;
}
else if(t6_l_n>=0){
REAL t11_l_n_1=t5_l_n*t6_h;
REAL t11_h_1=t5_l_n*t6_l_n;
t11_h=t11_h_1+t11_e;
t11_l_n=t11_l_n_1+t11_e;
}
else{
REAL t11_l_n_1=t5_l_n*t6_h;
REAL t5_h_n=-t5_h;
REAL t11_h_1=t5_h_n*t6_l_n;
t11_h=t11_h_1+t11_e;
t11_l_n=t11_l_n_1+t11_e;
}
}
else if(t5_l_n>=0){
if(t6_h<=0){
REAL t11_l_n_1=t5_h*t6_l_n;
REAL t11_h_1=t5_l_n*t6_l_n;
t11_h=t11_h_1+t11_e;
t11_l_n=t11_l_n_1+t11_e;
}
else if(t6_l_n>=0){
REAL t11_l_n_1_1=t5_l_n*t6_h;
REAL t11_l_n_1_2=t6_l_n*t5_h;
REAL t11_l_n_1=max(t11_l_n_1_1,t11_l_n_1_2);
REAL t11_h_1_1=t5_l_n*t6_l_n;
REAL t11_h_1_2=t5_h*t6_h;
REAL t11_h_1=max(t11_h_1_1,t11_h_1_2);
t11_h=t11_h_1+t11_e;
t11_l_n=t11_l_n_1+t11_e;
}
else{
REAL t11_l_n_1=t5_l_n*t6_h;
REAL t11_h_1=t5_h*t6_h;
t11_h=t11_h_1+t11_e;
t11_l_n=t11_l_n_1+t11_e;
}
}
else{
if(t6_h<=0){
REAL t11_l_n_1=t6_l_n*t5_h;
REAL t6_h_n=-t6_h;
REAL t11_h_1=t6_h_n*t5_l_n;
t11_h=t11_h_1+t11_e;
t11_l_n=t11_l_n_1+t11_e;
}
else if(t6_l_n>=0){
REAL t11_l_n_1=t6_l_n*t5_h;
REAL t11_h_1=t5_h*t6_h;
t11_h=t11_h_1+t11_e;
t11_l_n=t11_l_n_1+t11_e;
}
else{
REAL t5_l=-t5_l_n;
REAL t11_l_n_1=t5_l*t6_l_n;
REAL t11_h_1=t5_h*t6_h;
t11_h=t11_h_1+t11_e;
t11_l_n=t11_l_n_1+t11_e;
}
}
REAL t13_e_1=fabs(t13);
REAL t13_e=eps*t13_e_1;
REAL t13_h_1=t10_h+t11_l_n;
REAL t13_h=t13_h_1+t13_e;
REAL t13_l_n_1=t11_h+t10_l_n;
REAL t13_l_n=t13_l_n_1+t13_e;
REAL t14_e_1=fabs(t14);
REAL t14_e=eps*t14_e_1;
REAL  t14_h, t14_l_n;
if(t13_h<=0){
if(t9_h<=0){
REAL t13_h_n=-t13_h;
REAL t14_l_n_1=t13_h_n*t9_h;
REAL t14_h_1=t13_l_n*t9_l_n;
t14_h=t14_h_1+t14_e;
t14_l_n=t14_l_n_1+t14_e;
}
else if(t9_l_n>=0){
REAL t14_l_n_1=t13_l_n*t9_h;
REAL t14_h_1=t13_l_n*t9_l_n;
t14_h=t14_h_1+t14_e;
t14_l_n=t14_l_n_1+t14_e;
}
else{
REAL t14_l_n_1=t13_l_n*t9_h;
REAL t13_h_n=-t13_h;
REAL t14_h_1=t13_h_n*t9_l_n;
t14_h=t14_h_1+t14_e;
t14_l_n=t14_l_n_1+t14_e;
}
}
else if(t13_l_n>=0){
if(t9_h<=0){
REAL t14_l_n_1=t13_h*t9_l_n;
REAL t14_h_1=t13_l_n*t9_l_n;
t14_h=t14_h_1+t14_e;
t14_l_n=t14_l_n_1+t14_e;
}
else if(t9_l_n>=0){
REAL t14_l_n_1_1=t13_l_n*t9_h;
REAL t14_l_n_1_2=t9_l_n*t13_h;
REAL t14_l_n_1=max(t14_l_n_1_1,t14_l_n_1_2);
REAL t14_h_1_1=t13_l_n*t9_l_n;
REAL t14_h_1_2=t13_h*t9_h;
REAL t14_h_1=max(t14_h_1_1,t14_h_1_2);
t14_h=t14_h_1+t14_e;
t14_l_n=t14_l_n_1+t14_e;
}
else{
REAL t14_l_n_1=t13_l_n*t9_h;
REAL t14_h_1=t13_h*t9_h;
t14_h=t14_h_1+t14_e;
t14_l_n=t14_l_n_1+t14_e;
}
}
else{
if(t9_h<=0){
REAL t14_l_n_1=t9_l_n*t13_h;
REAL t9_h_n=-t9_h;
REAL t14_h_1=t9_h_n*t13_l_n;
t14_h=t14_h_1+t14_e;
t14_l_n=t14_l_n_1+t14_e;
}
else if(t9_l_n>=0){
REAL t14_l_n_1=t9_l_n*t13_h;
REAL t14_h_1=t13_h*t9_h;
t14_h=t14_h_1+t14_e;
t14_l_n=t14_l_n_1+t14_e;
}
else{
REAL t13_l=-t13_l_n;
REAL t14_l_n_1=t13_l*t9_l_n;
REAL t14_h_1=t13_h*t9_h;
t14_h=t14_h_1+t14_e;
t14_l_n=t14_l_n_1+t14_e;
}
}
REAL t15_e_1=fabs(t15);
REAL t15_e=eps*t15_e_1;
REAL t15_h=t15+t15_e;
REAL t15_l_n=t15_e-t15;
REAL t16_e_1=fabs(t16);
REAL t16_e=eps*t16_e_1;
REAL  t16_h, t16_l_n;
if(t5_h<=0){
if(t7_h<=0){
REAL t5_h_n=-t5_h;
REAL t16_l_n_1=t5_h_n*t7_h;
REAL t16_h_1=t5_l_n*t7_l_n;
t16_h=t16_h_1+t16_e;
t16_l_n=t16_l_n_1+t16_e;
}
else if(t7_l_n>=0){
REAL t16_l_n_1=t5_l_n*t7_h;
REAL t16_h_1=t5_l_n*t7_l_n;
t16_h=t16_h_1+t16_e;
t16_l_n=t16_l_n_1+t16_e;
}
else{
REAL t16_l_n_1=t5_l_n*t7_h;
REAL t5_h_n=-t5_h;
REAL t16_h_1=t5_h_n*t7_l_n;
t16_h=t16_h_1+t16_e;
t16_l_n=t16_l_n_1+t16_e;
}
}
else if(t5_l_n>=0){
if(t7_h<=0){
REAL t16_l_n_1=t5_h*t7_l_n;
REAL t16_h_1=t5_l_n*t7_l_n;
t16_h=t16_h_1+t16_e;
t16_l_n=t16_l_n_1+t16_e;
}
else if(t7_l_n>=0){
REAL t16_l_n_1_1=t5_l_n*t7_h;
REAL t16_l_n_1_2=t7_l_n*t5_h;
REAL t16_l_n_1=max(t16_l_n_1_1,t16_l_n_1_2);
REAL t16_h_1_1=t5_l_n*t7_l_n;
REAL t16_h_1_2=t5_h*t7_h;
REAL t16_h_1=max(t16_h_1_1,t16_h_1_2);
t16_h=t16_h_1+t16_e;
t16_l_n=t16_l_n_1+t16_e;
}
else{
REAL t16_l_n_1=t5_l_n*t7_h;
REAL t16_h_1=t5_h*t7_h;
t16_h=t16_h_1+t16_e;
t16_l_n=t16_l_n_1+t16_e;
}
}
else{
if(t7_h<=0){
REAL t16_l_n_1=t7_l_n*t5_h;
REAL t7_h_n=-t7_h;
REAL t16_h_1=t7_h_n*t5_l_n;
t16_h=t16_h_1+t16_e;
t16_l_n=t16_l_n_1+t16_e;
}
else if(t7_l_n>=0){
REAL t16_l_n_1=t7_l_n*t5_h;
REAL t16_h_1=t5_h*t7_h;
t16_h=t16_h_1+t16_e;
t16_l_n=t16_l_n_1+t16_e;
}
else{
REAL t5_l=-t5_l_n;
REAL t16_l_n_1=t5_l*t7_l_n;
REAL t16_h_1=t5_h*t7_h;
t16_h=t16_h_1+t16_e;
t16_l_n=t16_l_n_1+t16_e;
}
}
REAL t17_e_1=fabs(t17);
REAL t17_e=eps*t17_e_1;
REAL  t17_h, t17_l_n;
if(t2_h<=0){
if(t8_h<=0){
REAL t2_h_n=-t2_h;
REAL t17_l_n_1=t2_h_n*t8_h;
REAL t17_h_1=t2_l_n*t8_l_n;
t17_h=t17_h_1+t17_e;
t17_l_n=t17_l_n_1+t17_e;
}
else if(t8_l_n>=0){
REAL t17_l_n_1=t2_l_n*t8_h;
REAL t17_h_1=t2_l_n*t8_l_n;
t17_h=t17_h_1+t17_e;
t17_l_n=t17_l_n_1+t17_e;
}
else{
REAL t17_l_n_1=t2_l_n*t8_h;
REAL t2_h_n=-t2_h;
REAL t17_h_1=t2_h_n*t8_l_n;
t17_h=t17_h_1+t17_e;
t17_l_n=t17_l_n_1+t17_e;
}
}
else if(t2_l_n>=0){
if(t8_h<=0){
REAL t17_l_n_1=t2_h*t8_l_n;
REAL t17_h_1=t2_l_n*t8_l_n;
t17_h=t17_h_1+t17_e;
t17_l_n=t17_l_n_1+t17_e;
}
else if(t8_l_n>=0){
REAL t17_l_n_1_1=t2_l_n*t8_h;
REAL t17_l_n_1_2=t8_l_n*t2_h;
REAL t17_l_n_1=max(t17_l_n_1_1,t17_l_n_1_2);
REAL t17_h_1_1=t2_l_n*t8_l_n;
REAL t17_h_1_2=t2_h*t8_h;
REAL t17_h_1=max(t17_h_1_1,t17_h_1_2);
t17_h=t17_h_1+t17_e;
t17_l_n=t17_l_n_1+t17_e;
}
else{
REAL t17_l_n_1=t2_l_n*t8_h;
REAL t17_h_1=t2_h*t8_h;
t17_h=t17_h_1+t17_e;
t17_l_n=t17_l_n_1+t17_e;
}
}
else{
if(t8_h<=0){
REAL t17_l_n_1=t8_l_n*t2_h;
REAL t8_h_n=-t8_h;
REAL t17_h_1=t8_h_n*t2_l_n;
t17_h=t17_h_1+t17_e;
t17_l_n=t17_l_n_1+t17_e;
}
else if(t8_l_n>=0){
REAL t17_l_n_1=t8_l_n*t2_h;
REAL t17_h_1=t2_h*t8_h;
t17_h=t17_h_1+t17_e;
t17_l_n=t17_l_n_1+t17_e;
}
else{
REAL t2_l=-t2_l_n;
REAL t17_l_n_1=t2_l*t8_l_n;
REAL t17_h_1=t2_h*t8_h;
t17_h=t17_h_1+t17_e;
t17_l_n=t17_l_n_1+t17_e;
}
}
REAL t19_e_1=fabs(t19);
REAL t19_e=eps*t19_e_1;
REAL t19_h_1=t16_h+t17_l_n;
REAL t19_h=t19_h_1+t19_e;
REAL t19_l_n_1=t17_h+t16_l_n;
REAL t19_l_n=t19_l_n_1+t19_e;
REAL t20_e_1=fabs(t20);
REAL t20_e=eps*t20_e_1;
REAL  t20_h, t20_l_n;
if(t15_h<=0){
if(t19_h<=0){
REAL t15_h_n=-t15_h;
REAL t20_l_n_1=t15_h_n*t19_h;
REAL t20_h_1=t15_l_n*t19_l_n;
t20_h=t20_h_1+t20_e;
t20_l_n=t20_l_n_1+t20_e;
}
else if(t19_l_n>=0){
REAL t20_l_n_1=t15_l_n*t19_h;
REAL t20_h_1=t15_l_n*t19_l_n;
t20_h=t20_h_1+t20_e;
t20_l_n=t20_l_n_1+t20_e;
}
else{
REAL t20_l_n_1=t15_l_n*t19_h;
REAL t15_h_n=-t15_h;
REAL t20_h_1=t15_h_n*t19_l_n;
t20_h=t20_h_1+t20_e;
t20_l_n=t20_l_n_1+t20_e;
}
}
else if(t15_l_n>=0){
if(t19_h<=0){
REAL t20_l_n_1=t15_h*t19_l_n;
REAL t20_h_1=t15_l_n*t19_l_n;
t20_h=t20_h_1+t20_e;
t20_l_n=t20_l_n_1+t20_e;
}
else if(t19_l_n>=0){
REAL t20_l_n_1_1=t15_l_n*t19_h;
REAL t20_l_n_1_2=t19_l_n*t15_h;
REAL t20_l_n_1=max(t20_l_n_1_1,t20_l_n_1_2);
REAL t20_h_1_1=t15_l_n*t19_l_n;
REAL t20_h_1_2=t15_h*t19_h;
REAL t20_h_1=max(t20_h_1_1,t20_h_1_2);
t20_h=t20_h_1+t20_e;
t20_l_n=t20_l_n_1+t20_e;
}
else{
REAL t20_l_n_1=t15_l_n*t19_h;
REAL t20_h_1=t15_h*t19_h;
t20_h=t20_h_1+t20_e;
t20_l_n=t20_l_n_1+t20_e;
}
}
else{
if(t19_h<=0){
REAL t20_l_n_1=t19_l_n*t15_h;
REAL t19_h_n=-t19_h;
REAL t20_h_1=t19_h_n*t15_l_n;
t20_h=t20_h_1+t20_e;
t20_l_n=t20_l_n_1+t20_e;
}
else if(t19_l_n>=0){
REAL t20_l_n_1=t19_l_n*t15_h;
REAL t20_h_1=t15_h*t19_h;
t20_h=t20_h_1+t20_e;
t20_l_n=t20_l_n_1+t20_e;
}
else{
REAL t15_l=-t15_l_n;
REAL t20_l_n_1=t15_l*t19_l_n;
REAL t20_h_1=t15_h*t19_h;
t20_h=t20_h_1+t20_e;
t20_l_n=t20_l_n_1+t20_e;
}
}
REAL t21_e_1=fabs(t21);
REAL t21_e=eps*t21_e_1;
REAL t21_h=t21+t21_e;
REAL t21_l_n=t21_e-t21;
REAL t22_e_1=fabs(t22);
REAL t22_e=eps*t22_e_1;
REAL  t22_h, t22_l_n;
if(t4_h<=0){
if(t7_h<=0){
REAL t4_h_n=-t4_h;
REAL t22_l_n_1=t4_h_n*t7_h;
REAL t22_h_1=t4_l_n*t7_l_n;
t22_h=t22_h_1+t22_e;
t22_l_n=t22_l_n_1+t22_e;
}
else if(t7_l_n>=0){
REAL t22_l_n_1=t4_l_n*t7_h;
REAL t22_h_1=t4_l_n*t7_l_n;
t22_h=t22_h_1+t22_e;
t22_l_n=t22_l_n_1+t22_e;
}
else{
REAL t22_l_n_1=t4_l_n*t7_h;
REAL t4_h_n=-t4_h;
REAL t22_h_1=t4_h_n*t7_l_n;
t22_h=t22_h_1+t22_e;
t22_l_n=t22_l_n_1+t22_e;
}
}
else if(t4_l_n>=0){
if(t7_h<=0){
REAL t22_l_n_1=t4_h*t7_l_n;
REAL t22_h_1=t4_l_n*t7_l_n;
t22_h=t22_h_1+t22_e;
t22_l_n=t22_l_n_1+t22_e;
}
else if(t7_l_n>=0){
REAL t22_l_n_1_1=t4_l_n*t7_h;
REAL t22_l_n_1_2=t7_l_n*t4_h;
REAL t22_l_n_1=max(t22_l_n_1_1,t22_l_n_1_2);
REAL t22_h_1_1=t4_l_n*t7_l_n;
REAL t22_h_1_2=t4_h*t7_h;
REAL t22_h_1=max(t22_h_1_1,t22_h_1_2);
t22_h=t22_h_1+t22_e;
t22_l_n=t22_l_n_1+t22_e;
}
else{
REAL t22_l_n_1=t4_l_n*t7_h;
REAL t22_h_1=t4_h*t7_h;
t22_h=t22_h_1+t22_e;
t22_l_n=t22_l_n_1+t22_e;
}
}
else{
if(t7_h<=0){
REAL t22_l_n_1=t7_l_n*t4_h;
REAL t7_h_n=-t7_h;
REAL t22_h_1=t7_h_n*t4_l_n;
t22_h=t22_h_1+t22_e;
t22_l_n=t22_l_n_1+t22_e;
}
else if(t7_l_n>=0){
REAL t22_l_n_1=t7_l_n*t4_h;
REAL t22_h_1=t4_h*t7_h;
t22_h=t22_h_1+t22_e;
t22_l_n=t22_l_n_1+t22_e;
}
else{
REAL t4_l=-t4_l_n;
REAL t22_l_n_1=t4_l*t7_l_n;
REAL t22_h_1=t4_h*t7_h;
t22_h=t22_h_1+t22_e;
t22_l_n=t22_l_n_1+t22_e;
}
}
REAL t23_e_1=fabs(t23);
REAL t23_e=eps*t23_e_1;
REAL  t23_h, t23_l_n;
if(t6_h<=0){
if(t8_h<=0){
REAL t6_h_n=-t6_h;
REAL t23_l_n_1=t6_h_n*t8_h;
REAL t23_h_1=t6_l_n*t8_l_n;
t23_h=t23_h_1+t23_e;
t23_l_n=t23_l_n_1+t23_e;
}
else if(t8_l_n>=0){
REAL t23_l_n_1=t6_l_n*t8_h;
REAL t23_h_1=t6_l_n*t8_l_n;
t23_h=t23_h_1+t23_e;
t23_l_n=t23_l_n_1+t23_e;
}
else{
REAL t23_l_n_1=t6_l_n*t8_h;
REAL t6_h_n=-t6_h;
REAL t23_h_1=t6_h_n*t8_l_n;
t23_h=t23_h_1+t23_e;
t23_l_n=t23_l_n_1+t23_e;
}
}
else if(t6_l_n>=0){
if(t8_h<=0){
REAL t23_l_n_1=t6_h*t8_l_n;
REAL t23_h_1=t6_l_n*t8_l_n;
t23_h=t23_h_1+t23_e;
t23_l_n=t23_l_n_1+t23_e;
}
else if(t8_l_n>=0){
REAL t23_l_n_1_1=t6_l_n*t8_h;
REAL t23_l_n_1_2=t8_l_n*t6_h;
REAL t23_l_n_1=max(t23_l_n_1_1,t23_l_n_1_2);
REAL t23_h_1_1=t6_l_n*t8_l_n;
REAL t23_h_1_2=t6_h*t8_h;
REAL t23_h_1=max(t23_h_1_1,t23_h_1_2);
t23_h=t23_h_1+t23_e;
t23_l_n=t23_l_n_1+t23_e;
}
else{
REAL t23_l_n_1=t6_l_n*t8_h;
REAL t23_h_1=t6_h*t8_h;
t23_h=t23_h_1+t23_e;
t23_l_n=t23_l_n_1+t23_e;
}
}
else{
if(t8_h<=0){
REAL t23_l_n_1=t8_l_n*t6_h;
REAL t8_h_n=-t8_h;
REAL t23_h_1=t8_h_n*t6_l_n;
t23_h=t23_h_1+t23_e;
t23_l_n=t23_l_n_1+t23_e;
}
else if(t8_l_n>=0){
REAL t23_l_n_1=t8_l_n*t6_h;
REAL t23_h_1=t6_h*t8_h;
t23_h=t23_h_1+t23_e;
t23_l_n=t23_l_n_1+t23_e;
}
else{
REAL t6_l=-t6_l_n;
REAL t23_l_n_1=t6_l*t8_l_n;
REAL t23_h_1=t6_h*t8_h;
t23_h=t23_h_1+t23_e;
t23_l_n=t23_l_n_1+t23_e;
}
}
REAL t25_e_1=fabs(t25);
REAL t25_e=eps*t25_e_1;
REAL t25_h_1=t22_h+t23_l_n;
REAL t25_h=t25_h_1+t25_e;
REAL t25_l_n_1=t23_h+t22_l_n;
REAL t25_l_n=t25_l_n_1+t25_e;
REAL t26_e_1=fabs(t26);
REAL t26_e=eps*t26_e_1;
REAL  t26_h, t26_l_n;
if(t21_h<=0){
if(t25_h<=0){
REAL t21_h_n=-t21_h;
REAL t26_l_n_1=t21_h_n*t25_h;
REAL t26_h_1=t21_l_n*t25_l_n;
t26_h=t26_h_1+t26_e;
t26_l_n=t26_l_n_1+t26_e;
}
else if(t25_l_n>=0){
REAL t26_l_n_1=t21_l_n*t25_h;
REAL t26_h_1=t21_l_n*t25_l_n;
t26_h=t26_h_1+t26_e;
t26_l_n=t26_l_n_1+t26_e;
}
else{
REAL t26_l_n_1=t21_l_n*t25_h;
REAL t21_h_n=-t21_h;
REAL t26_h_1=t21_h_n*t25_l_n;
t26_h=t26_h_1+t26_e;
t26_l_n=t26_l_n_1+t26_e;
}
}
else if(t21_l_n>=0){
if(t25_h<=0){
REAL t26_l_n_1=t21_h*t25_l_n;
REAL t26_h_1=t21_l_n*t25_l_n;
t26_h=t26_h_1+t26_e;
t26_l_n=t26_l_n_1+t26_e;
}
else if(t25_l_n>=0){
REAL t26_l_n_1_1=t21_l_n*t25_h;
REAL t26_l_n_1_2=t25_l_n*t21_h;
REAL t26_l_n_1=max(t26_l_n_1_1,t26_l_n_1_2);
REAL t26_h_1_1=t21_l_n*t25_l_n;
REAL t26_h_1_2=t21_h*t25_h;
REAL t26_h_1=max(t26_h_1_1,t26_h_1_2);
t26_h=t26_h_1+t26_e;
t26_l_n=t26_l_n_1+t26_e;
}
else{
REAL t26_l_n_1=t21_l_n*t25_h;
REAL t26_h_1=t21_h*t25_h;
t26_h=t26_h_1+t26_e;
t26_l_n=t26_l_n_1+t26_e;
}
}
else{
if(t25_h<=0){
REAL t26_l_n_1=t25_l_n*t21_h;
REAL t25_h_n=-t25_h;
REAL t26_h_1=t25_h_n*t21_l_n;
t26_h=t26_h_1+t26_e;
t26_l_n=t26_l_n_1+t26_e;
}
else if(t25_l_n>=0){
REAL t26_l_n_1=t25_l_n*t21_h;
REAL t26_h_1=t21_h*t25_h;
t26_h=t26_h_1+t26_e;
t26_l_n=t26_l_n_1+t26_e;
}
else{
REAL t21_l=-t21_l_n;
REAL t26_l_n_1=t21_l*t25_l_n;
REAL t26_h_1=t21_h*t25_h;
t26_h=t26_h_1+t26_e;
t26_l_n=t26_l_n_1+t26_e;
}
}
REAL t28_e_1=fabs(t28);
REAL t28_e=eps*t28_e_1;
REAL t28_h_1=t14_h+t20_h;
REAL t28_h=t28_h_1+t28_e;
REAL t28_l_n_1=t14_l_n+t20_l_n;
REAL t28_l_n=t28_l_n_1+t28_e;
REAL t29_e_1=fabs(t29);
REAL t29_e=eps*t29_e_1;
REAL t29_h_1=t28_h+t26_l_n;
REAL t29_h=t29_h_1+t29_e;
REAL t29_l_n_1=t26_h+t28_l_n;
REAL t29_l_n=t29_l_n_1+t29_e;
fesetround(roundmode);
*result=t29;
if (inferr==0) {
*hbound=t29_h;
*lbound=-t29_l_n;
}
else{
*hbound=1/0;
*lbound=1/0;
}
}
//Main test Function

int main(){
int range=100;
int sample=100000;
REAL varlist;
REAL res, err, merr, herr,lerr;
FILE *fsave;
fsave=fopen("ProbSave.txt","w");
FILE *freport;
freport=fopen("ProbReport.txt","w");
int cIA=0;
int cFPS=0;
int cMan=0;
int cIAFPS=0;
int cIAnFPS=0;
int cFPSnIA=0;
int cManFPS=0;
int cMannFPS=0;
int cFPSnMan=0;
int bIA, bFPS, bMan;
REAL ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz;
int i;
for(i=0; i<sample; i++)
{

ax=randfl(range);
ay=randfl(range);
az=randfl(range);
bx=randfl(range);
by=randfl(range);
bz=randfl(range);
cx=randfl(range);
cy=randfl(range);
cz=randfl(range);
dx=randfl(range);
dy=randfl(range);
dz=randfl(range);


FPSyn_Orient3D(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,&res,&err);
Manual_Orient3D(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,&res,&merr);
IA_Orient3D(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,&res,&herr,&lerr);

printf("%e      %e      %e      %e      %e      %e         %f         %f\n ", res, err, merr, herr, lerr,herr-lerr, (herr-lerr)/2/err, merr/err);
fprintf(fsave, "%e      %e      %e      %e      %e      %e         %f         %f\n ", res, err, merr, herr, lerr,herr-lerr, (herr-lerr)/2/err, merr/err);

if ((0<herr)&&(0>lerr)){bIA=1;}
else { bIA=0;}

if ((res<merr)&&(res>-merr)){bMan=1;}
else { bMan=0;}

if ((res<err)&&(res>-err)){bFPS=1;}
else { bFPS=0;}

cIA = cIA + bIA;
cFPS = cFPS+ bFPS;
cMan= cMan + bMan;

cIAFPS= cIAFPS+bIA*bFPS;
if ((bIA==1)&&(bFPS==0)) {cIAnFPS=cIAnFPS+1;}
if ((bIA==0)&&(bFPS==1)) {cFPSnIA=cFPSnIA+1;}

cManFPS= cManFPS+bMan*bFPS;
if ((bMan==1)&&(bFPS==0)) {cMannFPS=cMannFPS+1;}
if ((bMan==0)&&(bFPS==1)) {cFPSnMan=cFPSnMan+1;}

}

fclose(fsave);

printf("In error IA: %d\n", cIA);
printf("In error FPSyn: %d\n", cFPS);
printf("In error Manual: %d\n", cMan);

fprintf(freport,"In error IA: %d\n", cIA);
fprintf(freport,"In error Manual: %d\n", cMan);
fprintf(freport,"In error FPSyn: %d\n", cFPS);

fprintf(freport,"In error IA & FPSyn: %d\n", cIAFPS);
fprintf(freport,"In error FPSyn not IA: %d\n", cFPSnIA);
fprintf(freport,"In error IA not FPSyn: %d\n", cIAnFPS);

fprintf(freport,"In error Manual & FPSyn: %d\n", cManFPS);
fprintf(freport,"In error FPSyn not Manual: %d\n", cFPSnMan);
fprintf(freport,"In error Manual not FPSyn: %d\n", cMannFPS);

fclose(freport);

}