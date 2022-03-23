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
void FPSyn_Orient2D(REAL ax,REAL ay,REAL bx,REAL by,REAL cx,REAL cy,REAL* result, REAL* error_bound){
REAL t4 = (ax-cx)*(by-cy);
REAL t7 = (ay-cy)*(bx-cx);
*result = (t4-t7);
*error_bound = 1.0000000001*(4.44089209850063e-16*((Absolute(t4))+(Absolute(t7))));
}//Manual Function
void Manual_Orient2D(REAL ax,REAL ay,REAL bx,REAL by,REAL cx,REAL cy,REAL* result, REAL* error_bound){
  REAL detleft, detright, det;
  REAL detsum;

  detleft = (ax - cx) * (by - cy);
  detright = (ay - cy) * (bx - cx);
  det = detleft - detright;
  detsum = Absolute(detleft) + Absolute(detright);

*result=det;
*error_bound=3.3306690738754716e-16*detsum;
}
//IA Function
void IA_Orient2D(REAL ax,REAL ay,REAL bx,REAL by,REAL cx,REAL cy,REAL* result, REAL* hbound, REAL* lbound){
REAL t2=ax-cx;
REAL t3=by-cy;
REAL t4=t2*t3;
REAL t5=ay-cy;
REAL t6=bx-cx;
REAL t7=t5*t6;
REAL t9=t4-t7;
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
REAL  t4_h, t4_l_n;
if(t2_h<=0){
if(t3_h<=0){
REAL t2_h_n=-t2_h;
REAL t4_l_n_1=t2_h_n*t3_h;
REAL t4_h_1=t2_l_n*t3_l_n;
t4_h=t4_h_1+t4_e;
t4_l_n=t4_l_n_1+t4_e;
}
else if(t3_l_n>=0){
REAL t4_l_n_1=t2_l_n*t3_h;
REAL t4_h_1=t2_l_n*t3_l_n;
t4_h=t4_h_1+t4_e;
t4_l_n=t4_l_n_1+t4_e;
}
else{
REAL t4_l_n_1=t2_l_n*t3_h;
REAL t2_h_n=-t2_h;
REAL t4_h_1=t2_h_n*t3_l_n;
t4_h=t4_h_1+t4_e;
t4_l_n=t4_l_n_1+t4_e;
}
}
else if(t2_l_n>=0){
if(t3_h<=0){
REAL t4_l_n_1=t2_h*t3_l_n;
REAL t4_h_1=t2_l_n*t3_l_n;
t4_h=t4_h_1+t4_e;
t4_l_n=t4_l_n_1+t4_e;
}
else if(t3_l_n>=0){
REAL t4_l_n_1_1=t2_l_n*t3_h;
REAL t4_l_n_1_2=t3_l_n*t2_h;
REAL t4_l_n_1=max(t4_l_n_1_1,t4_l_n_1_2);
REAL t4_h_1_1=t2_l_n*t3_l_n;
REAL t4_h_1_2=t2_h*t3_h;
REAL t4_h_1=max(t4_h_1_1,t4_h_1_2);
t4_h=t4_h_1+t4_e;
t4_l_n=t4_l_n_1+t4_e;
}
else{
REAL t4_l_n_1=t2_l_n*t3_h;
REAL t4_h_1=t2_h*t3_h;
t4_h=t4_h_1+t4_e;
t4_l_n=t4_l_n_1+t4_e;
}
}
else{
if(t3_h<=0){
REAL t4_l_n_1=t3_l_n*t2_h;
REAL t3_h_n=-t3_h;
REAL t4_h_1=t3_h_n*t2_l_n;
t4_h=t4_h_1+t4_e;
t4_l_n=t4_l_n_1+t4_e;
}
else if(t3_l_n>=0){
REAL t4_l_n_1=t3_l_n*t2_h;
REAL t4_h_1=t2_h*t3_h;
t4_h=t4_h_1+t4_e;
t4_l_n=t4_l_n_1+t4_e;
}
else{
REAL t2_l=-t2_l_n;
REAL t4_l_n_1=t2_l*t3_l_n;
REAL t4_h_1=t2_h*t3_h;
t4_h=t4_h_1+t4_e;
t4_l_n=t4_l_n_1+t4_e;
}
}
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
REAL  t7_h, t7_l_n;
if(t5_h<=0){
if(t6_h<=0){
REAL t5_h_n=-t5_h;
REAL t7_l_n_1=t5_h_n*t6_h;
REAL t7_h_1=t5_l_n*t6_l_n;
t7_h=t7_h_1+t7_e;
t7_l_n=t7_l_n_1+t7_e;
}
else if(t6_l_n>=0){
REAL t7_l_n_1=t5_l_n*t6_h;
REAL t7_h_1=t5_l_n*t6_l_n;
t7_h=t7_h_1+t7_e;
t7_l_n=t7_l_n_1+t7_e;
}
else{
REAL t7_l_n_1=t5_l_n*t6_h;
REAL t5_h_n=-t5_h;
REAL t7_h_1=t5_h_n*t6_l_n;
t7_h=t7_h_1+t7_e;
t7_l_n=t7_l_n_1+t7_e;
}
}
else if(t5_l_n>=0){
if(t6_h<=0){
REAL t7_l_n_1=t5_h*t6_l_n;
REAL t7_h_1=t5_l_n*t6_l_n;
t7_h=t7_h_1+t7_e;
t7_l_n=t7_l_n_1+t7_e;
}
else if(t6_l_n>=0){
REAL t7_l_n_1_1=t5_l_n*t6_h;
REAL t7_l_n_1_2=t6_l_n*t5_h;
REAL t7_l_n_1=max(t7_l_n_1_1,t7_l_n_1_2);
REAL t7_h_1_1=t5_l_n*t6_l_n;
REAL t7_h_1_2=t5_h*t6_h;
REAL t7_h_1=max(t7_h_1_1,t7_h_1_2);
t7_h=t7_h_1+t7_e;
t7_l_n=t7_l_n_1+t7_e;
}
else{
REAL t7_l_n_1=t5_l_n*t6_h;
REAL t7_h_1=t5_h*t6_h;
t7_h=t7_h_1+t7_e;
t7_l_n=t7_l_n_1+t7_e;
}
}
else{
if(t6_h<=0){
REAL t7_l_n_1=t6_l_n*t5_h;
REAL t6_h_n=-t6_h;
REAL t7_h_1=t6_h_n*t5_l_n;
t7_h=t7_h_1+t7_e;
t7_l_n=t7_l_n_1+t7_e;
}
else if(t6_l_n>=0){
REAL t7_l_n_1=t6_l_n*t5_h;
REAL t7_h_1=t5_h*t6_h;
t7_h=t7_h_1+t7_e;
t7_l_n=t7_l_n_1+t7_e;
}
else{
REAL t5_l=-t5_l_n;
REAL t7_l_n_1=t5_l*t6_l_n;
REAL t7_h_1=t5_h*t6_h;
t7_h=t7_h_1+t7_e;
t7_l_n=t7_l_n_1+t7_e;
}
}
REAL t9_e_1=fabs(t9);
REAL t9_e=eps*t9_e_1;
REAL t9_h_1=t4_h+t7_l_n;
REAL t9_h=t9_h_1+t9_e;
REAL t9_l_n_1=t7_h+t4_l_n;
REAL t9_l_n=t9_l_n_1+t9_e;
fesetround(roundmode);
*result=t9;
if (inferr==0) {
*hbound=t9_h;
*lbound=-t9_l_n;
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
REAL ax,ay,bx,by,cx,cy;
int i;
for(i=0; i<sample; i++)
{

ax=randfl(range);
ay=randfl(range);
bx=randfl(range);
by=randfl(range);
cx=randfl(range);
cy=randfl(range);


FPSyn_Orient2D(ax,ay,bx,by,cx,cy,&res,&err);
Manual_Orient2D(ax,ay,bx,by,cx,cy,&res,&merr);
IA_Orient2D(ax,ay,bx,by,cx,cy,&res,&herr,&lerr);

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