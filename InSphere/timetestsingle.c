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
void FPSyn_InSphere(REAL ax,REAL ay,REAL az,REAL bx,REAL by,REAL bz,REAL cx,REAL cy,REAL cz,REAL dx,REAL dy,REAL dz,REAL ex,REAL ey,REAL ez,REAL* result, REAL* error_bound){
//Calculation steps for input: ;
REAL t1=dx-ex;
REAL t3=dy-ey;
REAL t5=dz-ez;
REAL t6=az-ez;
REAL t7=bx-ex;
REAL t8=cy-ey;
REAL t9=by-ey;
REAL t10=cx-ex;
REAL t11=t7*t8;
REAL t12=t10*t9;
REAL t14=t11-t12;
REAL t15=cz-ez;
REAL t16=ax-ex;
REAL t17=ay-ey;
REAL t18=t16*t9;
REAL t19=t17*t7;
REAL t21=t18-t19;
REAL t22=bz-ez;
REAL t23=t16*t8;
REAL t24=t10*t17;
REAL t26=t23-t24;
REAL t27=t10*t3;
REAL t28=t1*t8;
REAL t30=t27-t28;
REAL t31=t16*t3;
REAL t32=t1*t17;
REAL t34=t31-t32;
REAL t35=t3*t7;
REAL t36=t1*t9;
REAL t38=t35-t36;
REAL t39=t1*t1;
REAL t40=t3*t3;
REAL t41=t5*t5;
REAL t42=t39+t40;
REAL t43=t41+t42;
REAL t44=t14*t6;
REAL t45=t15*t21;
REAL t46=t22*t26;
REAL t48=t44+t45;
REAL t49=t48-t46;
REAL t50=t43*t49;
REAL t51=t22*t22;
REAL t52=t7*t7;
REAL t53=t9*t9;
REAL t54=t51+t52;
REAL t55=t53+t54;
REAL t56=t26*t5;
REAL t57=t30*t6;
REAL t58=t15*t34;
REAL t60=t56+t57;
REAL t61=t60-t58;
REAL t62=t55*t61;
REAL t63=t10*t10;
REAL t64=t15*t15;
REAL t65=t8*t8;
REAL t66=t63+t64;
REAL t67=t65+t66;
REAL t68=t21*t5;
REAL t69=t38*t6;
REAL t70=t22*t34;
REAL t72=t68+t69;
REAL t73=t72-t70;
REAL t74=t67*t73;
REAL t76=t16*t16;
REAL t77=t17*t17;
REAL t78=t6*t6;
REAL t79=t76+t77;
REAL t80=t78+t79;
REAL t81=t14*t5;
REAL t82=t22*t30;
REAL t83=t15*t38;
REAL t85=t81+t82;
REAL t86=t85-t83;
REAL t87=t80*t86;
REAL t89=t50+t62;
REAL t90=t89-t74;
REAL t91=t90-t87;
//Calculation steps for Error-bound: ;
t18=fabs(t18);
t31=fabs(t31);
t24=fabs(t24);
t28=fabs(t28);
t87=fabs(t87);
t32=fabs(t32);
t15=fabs(t15);
t22=fabs(t22);
t62=fabs(t62);
t50=fabs(t50);
t5=fabs(t5);
t36=fabs(t36);
t19=fabs(t19);
t23=fabs(t23);
t27=fabs(t27);
t35=fabs(t35);
t11=fabs(t11);
t12=fabs(t12);
t74=fabs(t74);
t6=fabs(t6);
REAL tem0=t50+t62;
REAL tem1=t74+tem0;
REAL tem2=t87+tem1;
REAL tem3=9.99200722162642e-16*tem2;
REAL tem4=t18+t19;
REAL tem5=t23+t24;
REAL tem6=t11+t12;
REAL tem7=t31+t32;
REAL tem8=t27+t28;
REAL tem9=t35+t36;
REAL tem10=t15*tem4;
REAL tem11=t22*tem5;
REAL tem12=t6*tem6;
REAL tem13=tem10+tem11;
REAL tem14=tem12+tem13;
REAL tem15=t43*tem14;
REAL tem16=8.88178419700127e-16*tem15;
REAL tem17=t15*tem7;
REAL tem18=t5*tem5;
REAL tem19=t6*tem8;
REAL tem20=tem17+tem18;
REAL tem21=tem19+tem20;
REAL tem22=t55*tem21;
REAL tem23=8.88178419700127e-16*tem22;
REAL tem24=t22*tem7;
REAL tem25=t5*tem4;
REAL tem26=t6*tem9;
REAL tem27=tem24+tem25;
REAL tem28=tem26+tem27;
REAL tem29=t67*tem28;
REAL tem30=8.88178419700127e-16*tem29;
REAL tem31=t15*tem9;
REAL tem32=t22*tem8;
REAL tem33=t5*tem6;
REAL tem34=tem31+tem32;
REAL tem35=tem33+tem34;
REAL tem36=t80*tem35;
REAL tem37=8.88178419700127e-16*tem36;
REAL tem38=tem16+tem23;
REAL tem39=tem30+tem38;
REAL tem40=tem37+tem39;
REAL tem41=tem3+tem40;
*result=t91;
*error_bound=1.0000000001*tem41;
}
//Manual Function
void Manual_InSphere(REAL ax, REAL ay,REAL az,REAL bx,REAL by,REAL bz,REAL cx,REAL cy,REAL cz,REAL dx,REAL dy,REAL dz,REAL ex,REAL ey,REAL ez, REAL* result, REAL* errbound){

  REAL aex, bex, cex, dex;
  REAL aey, bey, cey, dey;
  REAL aez, bez, cez, dez;
  REAL aexbey, bexaey, bexcey, cexbey, cexdey, dexcey, dexaey, aexdey;
  REAL aexcey, cexaey, bexdey, dexbey;
  REAL alift, blift, clift, dlift;
  REAL ab, bc, cd, da, ac, bd;
  REAL abc, bcd, cda, dab;
  REAL aezplus, bezplus, cezplus, dezplus;
  REAL aexbeyplus, bexaeyplus, bexceyplus, cexbeyplus;
  REAL cexdeyplus, dexceyplus, dexaeyplus, aexdeyplus;
  REAL aexceyplus, cexaeyplus, bexdeyplus, dexbeyplus;
  REAL det;
  REAL permanent;

  aex = ax - ex;
  bex = bx - ex;
  cex = cx - ex;
  dex = dx - ex;
  
  aey = ay - ey;
  bey = by - ey;
  cey = cy - ey;
  dey = dy - ey;
  
  aez = az - ez;
  bez = bz - ez;
  cez = cz - ez;
  dez = dz - ez;


  aexbey = aex * bey;
  bexaey = bex * aey;
  ab = aexbey - bexaey;
  bexcey = bex * cey;
  cexbey = cex * bey;
  bc = bexcey - cexbey;
  cexdey = cex * dey;
  dexcey = dex * cey;
  cd = cexdey - dexcey;
  dexaey = dex * aey;
  aexdey = aex * dey;
  da = dexaey - aexdey;

  aexcey = aex * cey;
  cexaey = cex * aey;
  ac = aexcey - cexaey;
  bexdey = bex * dey;
  dexbey = dex * bey;
  bd = bexdey - dexbey;

  abc = aez * bc - bez * ac + cez * ab;
  bcd = bez * cd - cez * bd + dez * bc;
  cda = cez * da + dez * ac + aez * cd;
  dab = dez * ab + aez * bd + bez * da;

  alift = aex * aex + aey * aey + aez * aez;
  blift = bex * bex + bey * bey + bez * bez;
  clift = cex * cex + cey * cey + cez * cez;
  dlift = dex * dex + dey * dey + dez * dez;

  det = (dlift * abc - clift * dab) + (blift * cda - alift * bcd);

  aezplus = fabs(aez);
  bezplus = fabs(bez);
  cezplus = fabs(cez);
  dezplus = fabs(dez);
  aexbeyplus = fabs(aexbey);
  bexaeyplus = fabs(bexaey);
  bexceyplus = fabs(bexcey);
  cexbeyplus = fabs(cexbey);
  cexdeyplus = fabs(cexdey);
  dexceyplus = fabs(dexcey);
  dexaeyplus = fabs(dexaey);
  aexdeyplus = fabs(aexdey);
  aexceyplus = fabs(aexcey);
  cexaeyplus = fabs(cexaey);
  bexdeyplus = fabs(bexdey);
  dexbeyplus = fabs(dexbey);
  permanent = ((cexdeyplus + dexceyplus) * bezplus
               + (dexbeyplus + bexdeyplus) * cezplus
               + (bexceyplus + cexbeyplus) * dezplus)
            * alift
            + ((dexaeyplus + aexdeyplus) * cezplus
               + (aexceyplus + cexaeyplus) * dezplus
               + (cexdeyplus + dexceyplus) * aezplus)
            * blift
            + ((aexbeyplus + bexaeyplus) * dezplus
               + (bexdeyplus + dexbeyplus) * aezplus
               + (dexaeyplus + aexdeyplus) * bezplus)
            * clift
            + ((bexceyplus + cexbeyplus) * aezplus
               + (cexaeyplus + aexceyplus) * bezplus
               + (aexbeyplus + bexaeyplus) * cezplus)
            * dlift;
  *errbound = 1.7763568394002532e-15 * permanent;
  *result = det	;
}
//IA Function
void IA_InSphere(REAL ax,REAL ay,REAL az,REAL bx,REAL by,REAL bz,REAL cx,REAL cy,REAL cz,REAL dx,REAL dy,REAL dz,REAL ex,REAL ey,REAL ez,REAL* result, REAL* hbound, REAL* lbound){
REAL t1=dx-ex;
REAL t3=dy-ey;
REAL t5=dz-ez;
REAL t6=az-ez;
REAL t7=bx-ex;
REAL t8=cy-ey;
REAL t9=by-ey;
REAL t10=cx-ex;
REAL t11=t7*t8;
REAL t12=t10*t9;
REAL t14=t11-t12;
REAL t15=cz-ez;
REAL t16=ax-ex;
REAL t17=ay-ey;
REAL t18=t16*t9;
REAL t19=t17*t7;
REAL t21=t18-t19;
REAL t22=bz-ez;
REAL t23=t16*t8;
REAL t24=t10*t17;
REAL t26=t23-t24;
REAL t27=t10*t3;
REAL t28=t1*t8;
REAL t30=t27-t28;
REAL t31=t16*t3;
REAL t32=t1*t17;
REAL t34=t31-t32;
REAL t35=t3*t7;
REAL t36=t1*t9;
REAL t38=t35-t36;
REAL t39=t1*t1;
REAL t40=t3*t3;
REAL t41=t5*t5;
REAL t42=t39+t40;
REAL t43=t41+t42;
REAL t44=t14*t6;
REAL t45=t15*t21;
REAL t46=t22*t26;
REAL t48=t44+t45;
REAL t49=t48-t46;
REAL t50=t43*t49;
REAL t51=t22*t22;
REAL t52=t7*t7;
REAL t53=t9*t9;
REAL t54=t51+t52;
REAL t55=t53+t54;
REAL t56=t26*t5;
REAL t57=t30*t6;
REAL t58=t15*t34;
REAL t60=t56+t57;
REAL t61=t60-t58;
REAL t62=t55*t61;
REAL t63=t10*t10;
REAL t64=t15*t15;
REAL t65=t8*t8;
REAL t66=t63+t64;
REAL t67=t65+t66;
REAL t68=t21*t5;
REAL t69=t38*t6;
REAL t70=t22*t34;
REAL t72=t68+t69;
REAL t73=t72-t70;
REAL t74=t67*t73;
REAL t76=t16*t16;
REAL t77=t17*t17;
REAL t78=t6*t6;
REAL t79=t76+t77;
REAL t80=t78+t79;
REAL t81=t14*t5;
REAL t82=t22*t30;
REAL t83=t15*t38;
REAL t85=t81+t82;
REAL t86=t85-t83;
REAL t87=t80*t86;
REAL t89=t50+t62;
REAL t90=t89-t74;
REAL t91=t90-t87;
//Calculate bound
REAL eps=pow(2,-53);
int inferr=0;
int roundmode = fegetround();
fesetround(FE_UPWARD);
REAL t1_e_1=fabs(t1);
REAL t1_e=eps*t1_e_1;
REAL t1_h=t1+t1_e;
REAL t1_l_n=t1_e-t1;
REAL t3_e_1=fabs(t3);
REAL t3_e=eps*t3_e_1;
REAL t3_h=t3+t3_e;
REAL t3_l_n=t3_e-t3;
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
REAL  t11_h, t11_l_n;
if(t7_h<=0){
if(t8_h<=0){
REAL t7_h_n=-t7_h;
REAL t11_l_n_1=t7_h_n*t8_h;
REAL t11_h_1=t7_l_n*t8_l_n;
t11_h=t11_h_1+t11_e;
t11_l_n=t11_l_n_1+t11_e;
}
else if(t8_l_n>=0){
REAL t11_l_n_1=t7_l_n*t8_h;
REAL t11_h_1=t7_l_n*t8_l_n;
t11_h=t11_h_1+t11_e;
t11_l_n=t11_l_n_1+t11_e;
}
else{
REAL t11_l_n_1=t7_l_n*t8_h;
REAL t7_h_n=-t7_h;
REAL t11_h_1=t7_h_n*t8_l_n;
t11_h=t11_h_1+t11_e;
t11_l_n=t11_l_n_1+t11_e;
}
}
else if(t7_l_n>=0){
if(t8_h<=0){
REAL t11_l_n_1=t7_h*t8_l_n;
REAL t11_h_1=t7_l_n*t8_l_n;
t11_h=t11_h_1+t11_e;
t11_l_n=t11_l_n_1+t11_e;
}
else if(t8_l_n>=0){
REAL t11_l_n_1_1=t7_l_n*t8_h;
REAL t11_l_n_1_2=t8_l_n*t7_h;
REAL t11_l_n_1=max(t11_l_n_1_1,t11_l_n_1_2);
REAL t11_h_1_1=t7_l_n*t8_l_n;
REAL t11_h_1_2=t7_h*t8_h;
REAL t11_h_1=max(t11_h_1_1,t11_h_1_2);
t11_h=t11_h_1+t11_e;
t11_l_n=t11_l_n_1+t11_e;
}
else{
REAL t11_l_n_1=t7_l_n*t8_h;
REAL t11_h_1=t7_h*t8_h;
t11_h=t11_h_1+t11_e;
t11_l_n=t11_l_n_1+t11_e;
}
}
else{
if(t8_h<=0){
REAL t11_l_n_1=t8_l_n*t7_h;
REAL t8_h_n=-t8_h;
REAL t11_h_1=t8_h_n*t7_l_n;
t11_h=t11_h_1+t11_e;
t11_l_n=t11_l_n_1+t11_e;
}
else if(t8_l_n>=0){
REAL t11_l_n_1=t8_l_n*t7_h;
REAL t11_h_1=t7_h*t8_h;
t11_h=t11_h_1+t11_e;
t11_l_n=t11_l_n_1+t11_e;
}
else{
REAL t7_l=-t7_l_n;
REAL t11_l_n_1=t7_l*t8_l_n;
REAL t11_h_1=t7_h*t8_h;
t11_h=t11_h_1+t11_e;
t11_l_n=t11_l_n_1+t11_e;
}
}
REAL t12_e_1=fabs(t12);
REAL t12_e=eps*t12_e_1;
REAL  t12_h, t12_l_n;
if(t10_h<=0){
if(t9_h<=0){
REAL t10_h_n=-t10_h;
REAL t12_l_n_1=t10_h_n*t9_h;
REAL t12_h_1=t10_l_n*t9_l_n;
t12_h=t12_h_1+t12_e;
t12_l_n=t12_l_n_1+t12_e;
}
else if(t9_l_n>=0){
REAL t12_l_n_1=t10_l_n*t9_h;
REAL t12_h_1=t10_l_n*t9_l_n;
t12_h=t12_h_1+t12_e;
t12_l_n=t12_l_n_1+t12_e;
}
else{
REAL t12_l_n_1=t10_l_n*t9_h;
REAL t10_h_n=-t10_h;
REAL t12_h_1=t10_h_n*t9_l_n;
t12_h=t12_h_1+t12_e;
t12_l_n=t12_l_n_1+t12_e;
}
}
else if(t10_l_n>=0){
if(t9_h<=0){
REAL t12_l_n_1=t10_h*t9_l_n;
REAL t12_h_1=t10_l_n*t9_l_n;
t12_h=t12_h_1+t12_e;
t12_l_n=t12_l_n_1+t12_e;
}
else if(t9_l_n>=0){
REAL t12_l_n_1_1=t10_l_n*t9_h;
REAL t12_l_n_1_2=t9_l_n*t10_h;
REAL t12_l_n_1=max(t12_l_n_1_1,t12_l_n_1_2);
REAL t12_h_1_1=t10_l_n*t9_l_n;
REAL t12_h_1_2=t10_h*t9_h;
REAL t12_h_1=max(t12_h_1_1,t12_h_1_2);
t12_h=t12_h_1+t12_e;
t12_l_n=t12_l_n_1+t12_e;
}
else{
REAL t12_l_n_1=t10_l_n*t9_h;
REAL t12_h_1=t10_h*t9_h;
t12_h=t12_h_1+t12_e;
t12_l_n=t12_l_n_1+t12_e;
}
}
else{
if(t9_h<=0){
REAL t12_l_n_1=t9_l_n*t10_h;
REAL t9_h_n=-t9_h;
REAL t12_h_1=t9_h_n*t10_l_n;
t12_h=t12_h_1+t12_e;
t12_l_n=t12_l_n_1+t12_e;
}
else if(t9_l_n>=0){
REAL t12_l_n_1=t9_l_n*t10_h;
REAL t12_h_1=t10_h*t9_h;
t12_h=t12_h_1+t12_e;
t12_l_n=t12_l_n_1+t12_e;
}
else{
REAL t10_l=-t10_l_n;
REAL t12_l_n_1=t10_l*t9_l_n;
REAL t12_h_1=t10_h*t9_h;
t12_h=t12_h_1+t12_e;
t12_l_n=t12_l_n_1+t12_e;
}
}
REAL t14_e_1=fabs(t14);
REAL t14_e=eps*t14_e_1;
REAL t14_h_1=t11_h+t12_l_n;
REAL t14_h=t14_h_1+t14_e;
REAL t14_l_n_1=t12_h+t11_l_n;
REAL t14_l_n=t14_l_n_1+t14_e;
REAL t15_e_1=fabs(t15);
REAL t15_e=eps*t15_e_1;
REAL t15_h=t15+t15_e;
REAL t15_l_n=t15_e-t15;
REAL t16_e_1=fabs(t16);
REAL t16_e=eps*t16_e_1;
REAL t16_h=t16+t16_e;
REAL t16_l_n=t16_e-t16;
REAL t17_e_1=fabs(t17);
REAL t17_e=eps*t17_e_1;
REAL t17_h=t17+t17_e;
REAL t17_l_n=t17_e-t17;
REAL t18_e_1=fabs(t18);
REAL t18_e=eps*t18_e_1;
REAL  t18_h, t18_l_n;
if(t16_h<=0){
if(t9_h<=0){
REAL t16_h_n=-t16_h;
REAL t18_l_n_1=t16_h_n*t9_h;
REAL t18_h_1=t16_l_n*t9_l_n;
t18_h=t18_h_1+t18_e;
t18_l_n=t18_l_n_1+t18_e;
}
else if(t9_l_n>=0){
REAL t18_l_n_1=t16_l_n*t9_h;
REAL t18_h_1=t16_l_n*t9_l_n;
t18_h=t18_h_1+t18_e;
t18_l_n=t18_l_n_1+t18_e;
}
else{
REAL t18_l_n_1=t16_l_n*t9_h;
REAL t16_h_n=-t16_h;
REAL t18_h_1=t16_h_n*t9_l_n;
t18_h=t18_h_1+t18_e;
t18_l_n=t18_l_n_1+t18_e;
}
}
else if(t16_l_n>=0){
if(t9_h<=0){
REAL t18_l_n_1=t16_h*t9_l_n;
REAL t18_h_1=t16_l_n*t9_l_n;
t18_h=t18_h_1+t18_e;
t18_l_n=t18_l_n_1+t18_e;
}
else if(t9_l_n>=0){
REAL t18_l_n_1_1=t16_l_n*t9_h;
REAL t18_l_n_1_2=t9_l_n*t16_h;
REAL t18_l_n_1=max(t18_l_n_1_1,t18_l_n_1_2);
REAL t18_h_1_1=t16_l_n*t9_l_n;
REAL t18_h_1_2=t16_h*t9_h;
REAL t18_h_1=max(t18_h_1_1,t18_h_1_2);
t18_h=t18_h_1+t18_e;
t18_l_n=t18_l_n_1+t18_e;
}
else{
REAL t18_l_n_1=t16_l_n*t9_h;
REAL t18_h_1=t16_h*t9_h;
t18_h=t18_h_1+t18_e;
t18_l_n=t18_l_n_1+t18_e;
}
}
else{
if(t9_h<=0){
REAL t18_l_n_1=t9_l_n*t16_h;
REAL t9_h_n=-t9_h;
REAL t18_h_1=t9_h_n*t16_l_n;
t18_h=t18_h_1+t18_e;
t18_l_n=t18_l_n_1+t18_e;
}
else if(t9_l_n>=0){
REAL t18_l_n_1=t9_l_n*t16_h;
REAL t18_h_1=t16_h*t9_h;
t18_h=t18_h_1+t18_e;
t18_l_n=t18_l_n_1+t18_e;
}
else{
REAL t16_l=-t16_l_n;
REAL t18_l_n_1=t16_l*t9_l_n;
REAL t18_h_1=t16_h*t9_h;
t18_h=t18_h_1+t18_e;
t18_l_n=t18_l_n_1+t18_e;
}
}
REAL t19_e_1=fabs(t19);
REAL t19_e=eps*t19_e_1;
REAL  t19_h, t19_l_n;
if(t17_h<=0){
if(t7_h<=0){
REAL t17_h_n=-t17_h;
REAL t19_l_n_1=t17_h_n*t7_h;
REAL t19_h_1=t17_l_n*t7_l_n;
t19_h=t19_h_1+t19_e;
t19_l_n=t19_l_n_1+t19_e;
}
else if(t7_l_n>=0){
REAL t19_l_n_1=t17_l_n*t7_h;
REAL t19_h_1=t17_l_n*t7_l_n;
t19_h=t19_h_1+t19_e;
t19_l_n=t19_l_n_1+t19_e;
}
else{
REAL t19_l_n_1=t17_l_n*t7_h;
REAL t17_h_n=-t17_h;
REAL t19_h_1=t17_h_n*t7_l_n;
t19_h=t19_h_1+t19_e;
t19_l_n=t19_l_n_1+t19_e;
}
}
else if(t17_l_n>=0){
if(t7_h<=0){
REAL t19_l_n_1=t17_h*t7_l_n;
REAL t19_h_1=t17_l_n*t7_l_n;
t19_h=t19_h_1+t19_e;
t19_l_n=t19_l_n_1+t19_e;
}
else if(t7_l_n>=0){
REAL t19_l_n_1_1=t17_l_n*t7_h;
REAL t19_l_n_1_2=t7_l_n*t17_h;
REAL t19_l_n_1=max(t19_l_n_1_1,t19_l_n_1_2);
REAL t19_h_1_1=t17_l_n*t7_l_n;
REAL t19_h_1_2=t17_h*t7_h;
REAL t19_h_1=max(t19_h_1_1,t19_h_1_2);
t19_h=t19_h_1+t19_e;
t19_l_n=t19_l_n_1+t19_e;
}
else{
REAL t19_l_n_1=t17_l_n*t7_h;
REAL t19_h_1=t17_h*t7_h;
t19_h=t19_h_1+t19_e;
t19_l_n=t19_l_n_1+t19_e;
}
}
else{
if(t7_h<=0){
REAL t19_l_n_1=t7_l_n*t17_h;
REAL t7_h_n=-t7_h;
REAL t19_h_1=t7_h_n*t17_l_n;
t19_h=t19_h_1+t19_e;
t19_l_n=t19_l_n_1+t19_e;
}
else if(t7_l_n>=0){
REAL t19_l_n_1=t7_l_n*t17_h;
REAL t19_h_1=t17_h*t7_h;
t19_h=t19_h_1+t19_e;
t19_l_n=t19_l_n_1+t19_e;
}
else{
REAL t17_l=-t17_l_n;
REAL t19_l_n_1=t17_l*t7_l_n;
REAL t19_h_1=t17_h*t7_h;
t19_h=t19_h_1+t19_e;
t19_l_n=t19_l_n_1+t19_e;
}
}
REAL t21_e_1=fabs(t21);
REAL t21_e=eps*t21_e_1;
REAL t21_h_1=t18_h+t19_l_n;
REAL t21_h=t21_h_1+t21_e;
REAL t21_l_n_1=t19_h+t18_l_n;
REAL t21_l_n=t21_l_n_1+t21_e;
REAL t22_e_1=fabs(t22);
REAL t22_e=eps*t22_e_1;
REAL t22_h=t22+t22_e;
REAL t22_l_n=t22_e-t22;
REAL t23_e_1=fabs(t23);
REAL t23_e=eps*t23_e_1;
REAL  t23_h, t23_l_n;
if(t16_h<=0){
if(t8_h<=0){
REAL t16_h_n=-t16_h;
REAL t23_l_n_1=t16_h_n*t8_h;
REAL t23_h_1=t16_l_n*t8_l_n;
t23_h=t23_h_1+t23_e;
t23_l_n=t23_l_n_1+t23_e;
}
else if(t8_l_n>=0){
REAL t23_l_n_1=t16_l_n*t8_h;
REAL t23_h_1=t16_l_n*t8_l_n;
t23_h=t23_h_1+t23_e;
t23_l_n=t23_l_n_1+t23_e;
}
else{
REAL t23_l_n_1=t16_l_n*t8_h;
REAL t16_h_n=-t16_h;
REAL t23_h_1=t16_h_n*t8_l_n;
t23_h=t23_h_1+t23_e;
t23_l_n=t23_l_n_1+t23_e;
}
}
else if(t16_l_n>=0){
if(t8_h<=0){
REAL t23_l_n_1=t16_h*t8_l_n;
REAL t23_h_1=t16_l_n*t8_l_n;
t23_h=t23_h_1+t23_e;
t23_l_n=t23_l_n_1+t23_e;
}
else if(t8_l_n>=0){
REAL t23_l_n_1_1=t16_l_n*t8_h;
REAL t23_l_n_1_2=t8_l_n*t16_h;
REAL t23_l_n_1=max(t23_l_n_1_1,t23_l_n_1_2);
REAL t23_h_1_1=t16_l_n*t8_l_n;
REAL t23_h_1_2=t16_h*t8_h;
REAL t23_h_1=max(t23_h_1_1,t23_h_1_2);
t23_h=t23_h_1+t23_e;
t23_l_n=t23_l_n_1+t23_e;
}
else{
REAL t23_l_n_1=t16_l_n*t8_h;
REAL t23_h_1=t16_h*t8_h;
t23_h=t23_h_1+t23_e;
t23_l_n=t23_l_n_1+t23_e;
}
}
else{
if(t8_h<=0){
REAL t23_l_n_1=t8_l_n*t16_h;
REAL t8_h_n=-t8_h;
REAL t23_h_1=t8_h_n*t16_l_n;
t23_h=t23_h_1+t23_e;
t23_l_n=t23_l_n_1+t23_e;
}
else if(t8_l_n>=0){
REAL t23_l_n_1=t8_l_n*t16_h;
REAL t23_h_1=t16_h*t8_h;
t23_h=t23_h_1+t23_e;
t23_l_n=t23_l_n_1+t23_e;
}
else{
REAL t16_l=-t16_l_n;
REAL t23_l_n_1=t16_l*t8_l_n;
REAL t23_h_1=t16_h*t8_h;
t23_h=t23_h_1+t23_e;
t23_l_n=t23_l_n_1+t23_e;
}
}
REAL t24_e_1=fabs(t24);
REAL t24_e=eps*t24_e_1;
REAL  t24_h, t24_l_n;
if(t10_h<=0){
if(t17_h<=0){
REAL t10_h_n=-t10_h;
REAL t24_l_n_1=t10_h_n*t17_h;
REAL t24_h_1=t10_l_n*t17_l_n;
t24_h=t24_h_1+t24_e;
t24_l_n=t24_l_n_1+t24_e;
}
else if(t17_l_n>=0){
REAL t24_l_n_1=t10_l_n*t17_h;
REAL t24_h_1=t10_l_n*t17_l_n;
t24_h=t24_h_1+t24_e;
t24_l_n=t24_l_n_1+t24_e;
}
else{
REAL t24_l_n_1=t10_l_n*t17_h;
REAL t10_h_n=-t10_h;
REAL t24_h_1=t10_h_n*t17_l_n;
t24_h=t24_h_1+t24_e;
t24_l_n=t24_l_n_1+t24_e;
}
}
else if(t10_l_n>=0){
if(t17_h<=0){
REAL t24_l_n_1=t10_h*t17_l_n;
REAL t24_h_1=t10_l_n*t17_l_n;
t24_h=t24_h_1+t24_e;
t24_l_n=t24_l_n_1+t24_e;
}
else if(t17_l_n>=0){
REAL t24_l_n_1_1=t10_l_n*t17_h;
REAL t24_l_n_1_2=t17_l_n*t10_h;
REAL t24_l_n_1=max(t24_l_n_1_1,t24_l_n_1_2);
REAL t24_h_1_1=t10_l_n*t17_l_n;
REAL t24_h_1_2=t10_h*t17_h;
REAL t24_h_1=max(t24_h_1_1,t24_h_1_2);
t24_h=t24_h_1+t24_e;
t24_l_n=t24_l_n_1+t24_e;
}
else{
REAL t24_l_n_1=t10_l_n*t17_h;
REAL t24_h_1=t10_h*t17_h;
t24_h=t24_h_1+t24_e;
t24_l_n=t24_l_n_1+t24_e;
}
}
else{
if(t17_h<=0){
REAL t24_l_n_1=t17_l_n*t10_h;
REAL t17_h_n=-t17_h;
REAL t24_h_1=t17_h_n*t10_l_n;
t24_h=t24_h_1+t24_e;
t24_l_n=t24_l_n_1+t24_e;
}
else if(t17_l_n>=0){
REAL t24_l_n_1=t17_l_n*t10_h;
REAL t24_h_1=t10_h*t17_h;
t24_h=t24_h_1+t24_e;
t24_l_n=t24_l_n_1+t24_e;
}
else{
REAL t10_l=-t10_l_n;
REAL t24_l_n_1=t10_l*t17_l_n;
REAL t24_h_1=t10_h*t17_h;
t24_h=t24_h_1+t24_e;
t24_l_n=t24_l_n_1+t24_e;
}
}
REAL t26_e_1=fabs(t26);
REAL t26_e=eps*t26_e_1;
REAL t26_h_1=t23_h+t24_l_n;
REAL t26_h=t26_h_1+t26_e;
REAL t26_l_n_1=t24_h+t23_l_n;
REAL t26_l_n=t26_l_n_1+t26_e;
REAL t27_e_1=fabs(t27);
REAL t27_e=eps*t27_e_1;
REAL  t27_h, t27_l_n;
if(t10_h<=0){
if(t3_h<=0){
REAL t10_h_n=-t10_h;
REAL t27_l_n_1=t10_h_n*t3_h;
REAL t27_h_1=t10_l_n*t3_l_n;
t27_h=t27_h_1+t27_e;
t27_l_n=t27_l_n_1+t27_e;
}
else if(t3_l_n>=0){
REAL t27_l_n_1=t10_l_n*t3_h;
REAL t27_h_1=t10_l_n*t3_l_n;
t27_h=t27_h_1+t27_e;
t27_l_n=t27_l_n_1+t27_e;
}
else{
REAL t27_l_n_1=t10_l_n*t3_h;
REAL t10_h_n=-t10_h;
REAL t27_h_1=t10_h_n*t3_l_n;
t27_h=t27_h_1+t27_e;
t27_l_n=t27_l_n_1+t27_e;
}
}
else if(t10_l_n>=0){
if(t3_h<=0){
REAL t27_l_n_1=t10_h*t3_l_n;
REAL t27_h_1=t10_l_n*t3_l_n;
t27_h=t27_h_1+t27_e;
t27_l_n=t27_l_n_1+t27_e;
}
else if(t3_l_n>=0){
REAL t27_l_n_1_1=t10_l_n*t3_h;
REAL t27_l_n_1_2=t3_l_n*t10_h;
REAL t27_l_n_1=max(t27_l_n_1_1,t27_l_n_1_2);
REAL t27_h_1_1=t10_l_n*t3_l_n;
REAL t27_h_1_2=t10_h*t3_h;
REAL t27_h_1=max(t27_h_1_1,t27_h_1_2);
t27_h=t27_h_1+t27_e;
t27_l_n=t27_l_n_1+t27_e;
}
else{
REAL t27_l_n_1=t10_l_n*t3_h;
REAL t27_h_1=t10_h*t3_h;
t27_h=t27_h_1+t27_e;
t27_l_n=t27_l_n_1+t27_e;
}
}
else{
if(t3_h<=0){
REAL t27_l_n_1=t3_l_n*t10_h;
REAL t3_h_n=-t3_h;
REAL t27_h_1=t3_h_n*t10_l_n;
t27_h=t27_h_1+t27_e;
t27_l_n=t27_l_n_1+t27_e;
}
else if(t3_l_n>=0){
REAL t27_l_n_1=t3_l_n*t10_h;
REAL t27_h_1=t10_h*t3_h;
t27_h=t27_h_1+t27_e;
t27_l_n=t27_l_n_1+t27_e;
}
else{
REAL t10_l=-t10_l_n;
REAL t27_l_n_1=t10_l*t3_l_n;
REAL t27_h_1=t10_h*t3_h;
t27_h=t27_h_1+t27_e;
t27_l_n=t27_l_n_1+t27_e;
}
}
REAL t28_e_1=fabs(t28);
REAL t28_e=eps*t28_e_1;
REAL  t28_h, t28_l_n;
if(t1_h<=0){
if(t8_h<=0){
REAL t1_h_n=-t1_h;
REAL t28_l_n_1=t1_h_n*t8_h;
REAL t28_h_1=t1_l_n*t8_l_n;
t28_h=t28_h_1+t28_e;
t28_l_n=t28_l_n_1+t28_e;
}
else if(t8_l_n>=0){
REAL t28_l_n_1=t1_l_n*t8_h;
REAL t28_h_1=t1_l_n*t8_l_n;
t28_h=t28_h_1+t28_e;
t28_l_n=t28_l_n_1+t28_e;
}
else{
REAL t28_l_n_1=t1_l_n*t8_h;
REAL t1_h_n=-t1_h;
REAL t28_h_1=t1_h_n*t8_l_n;
t28_h=t28_h_1+t28_e;
t28_l_n=t28_l_n_1+t28_e;
}
}
else if(t1_l_n>=0){
if(t8_h<=0){
REAL t28_l_n_1=t1_h*t8_l_n;
REAL t28_h_1=t1_l_n*t8_l_n;
t28_h=t28_h_1+t28_e;
t28_l_n=t28_l_n_1+t28_e;
}
else if(t8_l_n>=0){
REAL t28_l_n_1_1=t1_l_n*t8_h;
REAL t28_l_n_1_2=t8_l_n*t1_h;
REAL t28_l_n_1=max(t28_l_n_1_1,t28_l_n_1_2);
REAL t28_h_1_1=t1_l_n*t8_l_n;
REAL t28_h_1_2=t1_h*t8_h;
REAL t28_h_1=max(t28_h_1_1,t28_h_1_2);
t28_h=t28_h_1+t28_e;
t28_l_n=t28_l_n_1+t28_e;
}
else{
REAL t28_l_n_1=t1_l_n*t8_h;
REAL t28_h_1=t1_h*t8_h;
t28_h=t28_h_1+t28_e;
t28_l_n=t28_l_n_1+t28_e;
}
}
else{
if(t8_h<=0){
REAL t28_l_n_1=t8_l_n*t1_h;
REAL t8_h_n=-t8_h;
REAL t28_h_1=t8_h_n*t1_l_n;
t28_h=t28_h_1+t28_e;
t28_l_n=t28_l_n_1+t28_e;
}
else if(t8_l_n>=0){
REAL t28_l_n_1=t8_l_n*t1_h;
REAL t28_h_1=t1_h*t8_h;
t28_h=t28_h_1+t28_e;
t28_l_n=t28_l_n_1+t28_e;
}
else{
REAL t1_l=-t1_l_n;
REAL t28_l_n_1=t1_l*t8_l_n;
REAL t28_h_1=t1_h*t8_h;
t28_h=t28_h_1+t28_e;
t28_l_n=t28_l_n_1+t28_e;
}
}
REAL t30_e_1=fabs(t30);
REAL t30_e=eps*t30_e_1;
REAL t30_h_1=t27_h+t28_l_n;
REAL t30_h=t30_h_1+t30_e;
REAL t30_l_n_1=t28_h+t27_l_n;
REAL t30_l_n=t30_l_n_1+t30_e;
REAL t31_e_1=fabs(t31);
REAL t31_e=eps*t31_e_1;
REAL  t31_h, t31_l_n;
if(t16_h<=0){
if(t3_h<=0){
REAL t16_h_n=-t16_h;
REAL t31_l_n_1=t16_h_n*t3_h;
REAL t31_h_1=t16_l_n*t3_l_n;
t31_h=t31_h_1+t31_e;
t31_l_n=t31_l_n_1+t31_e;
}
else if(t3_l_n>=0){
REAL t31_l_n_1=t16_l_n*t3_h;
REAL t31_h_1=t16_l_n*t3_l_n;
t31_h=t31_h_1+t31_e;
t31_l_n=t31_l_n_1+t31_e;
}
else{
REAL t31_l_n_1=t16_l_n*t3_h;
REAL t16_h_n=-t16_h;
REAL t31_h_1=t16_h_n*t3_l_n;
t31_h=t31_h_1+t31_e;
t31_l_n=t31_l_n_1+t31_e;
}
}
else if(t16_l_n>=0){
if(t3_h<=0){
REAL t31_l_n_1=t16_h*t3_l_n;
REAL t31_h_1=t16_l_n*t3_l_n;
t31_h=t31_h_1+t31_e;
t31_l_n=t31_l_n_1+t31_e;
}
else if(t3_l_n>=0){
REAL t31_l_n_1_1=t16_l_n*t3_h;
REAL t31_l_n_1_2=t3_l_n*t16_h;
REAL t31_l_n_1=max(t31_l_n_1_1,t31_l_n_1_2);
REAL t31_h_1_1=t16_l_n*t3_l_n;
REAL t31_h_1_2=t16_h*t3_h;
REAL t31_h_1=max(t31_h_1_1,t31_h_1_2);
t31_h=t31_h_1+t31_e;
t31_l_n=t31_l_n_1+t31_e;
}
else{
REAL t31_l_n_1=t16_l_n*t3_h;
REAL t31_h_1=t16_h*t3_h;
t31_h=t31_h_1+t31_e;
t31_l_n=t31_l_n_1+t31_e;
}
}
else{
if(t3_h<=0){
REAL t31_l_n_1=t3_l_n*t16_h;
REAL t3_h_n=-t3_h;
REAL t31_h_1=t3_h_n*t16_l_n;
t31_h=t31_h_1+t31_e;
t31_l_n=t31_l_n_1+t31_e;
}
else if(t3_l_n>=0){
REAL t31_l_n_1=t3_l_n*t16_h;
REAL t31_h_1=t16_h*t3_h;
t31_h=t31_h_1+t31_e;
t31_l_n=t31_l_n_1+t31_e;
}
else{
REAL t16_l=-t16_l_n;
REAL t31_l_n_1=t16_l*t3_l_n;
REAL t31_h_1=t16_h*t3_h;
t31_h=t31_h_1+t31_e;
t31_l_n=t31_l_n_1+t31_e;
}
}
REAL t32_e_1=fabs(t32);
REAL t32_e=eps*t32_e_1;
REAL  t32_h, t32_l_n;
if(t1_h<=0){
if(t17_h<=0){
REAL t1_h_n=-t1_h;
REAL t32_l_n_1=t1_h_n*t17_h;
REAL t32_h_1=t1_l_n*t17_l_n;
t32_h=t32_h_1+t32_e;
t32_l_n=t32_l_n_1+t32_e;
}
else if(t17_l_n>=0){
REAL t32_l_n_1=t1_l_n*t17_h;
REAL t32_h_1=t1_l_n*t17_l_n;
t32_h=t32_h_1+t32_e;
t32_l_n=t32_l_n_1+t32_e;
}
else{
REAL t32_l_n_1=t1_l_n*t17_h;
REAL t1_h_n=-t1_h;
REAL t32_h_1=t1_h_n*t17_l_n;
t32_h=t32_h_1+t32_e;
t32_l_n=t32_l_n_1+t32_e;
}
}
else if(t1_l_n>=0){
if(t17_h<=0){
REAL t32_l_n_1=t1_h*t17_l_n;
REAL t32_h_1=t1_l_n*t17_l_n;
t32_h=t32_h_1+t32_e;
t32_l_n=t32_l_n_1+t32_e;
}
else if(t17_l_n>=0){
REAL t32_l_n_1_1=t1_l_n*t17_h;
REAL t32_l_n_1_2=t17_l_n*t1_h;
REAL t32_l_n_1=max(t32_l_n_1_1,t32_l_n_1_2);
REAL t32_h_1_1=t1_l_n*t17_l_n;
REAL t32_h_1_2=t1_h*t17_h;
REAL t32_h_1=max(t32_h_1_1,t32_h_1_2);
t32_h=t32_h_1+t32_e;
t32_l_n=t32_l_n_1+t32_e;
}
else{
REAL t32_l_n_1=t1_l_n*t17_h;
REAL t32_h_1=t1_h*t17_h;
t32_h=t32_h_1+t32_e;
t32_l_n=t32_l_n_1+t32_e;
}
}
else{
if(t17_h<=0){
REAL t32_l_n_1=t17_l_n*t1_h;
REAL t17_h_n=-t17_h;
REAL t32_h_1=t17_h_n*t1_l_n;
t32_h=t32_h_1+t32_e;
t32_l_n=t32_l_n_1+t32_e;
}
else if(t17_l_n>=0){
REAL t32_l_n_1=t17_l_n*t1_h;
REAL t32_h_1=t1_h*t17_h;
t32_h=t32_h_1+t32_e;
t32_l_n=t32_l_n_1+t32_e;
}
else{
REAL t1_l=-t1_l_n;
REAL t32_l_n_1=t1_l*t17_l_n;
REAL t32_h_1=t1_h*t17_h;
t32_h=t32_h_1+t32_e;
t32_l_n=t32_l_n_1+t32_e;
}
}
REAL t34_e_1=fabs(t34);
REAL t34_e=eps*t34_e_1;
REAL t34_h_1=t31_h+t32_l_n;
REAL t34_h=t34_h_1+t34_e;
REAL t34_l_n_1=t32_h+t31_l_n;
REAL t34_l_n=t34_l_n_1+t34_e;
REAL t35_e_1=fabs(t35);
REAL t35_e=eps*t35_e_1;
REAL  t35_h, t35_l_n;
if(t3_h<=0){
if(t7_h<=0){
REAL t3_h_n=-t3_h;
REAL t35_l_n_1=t3_h_n*t7_h;
REAL t35_h_1=t3_l_n*t7_l_n;
t35_h=t35_h_1+t35_e;
t35_l_n=t35_l_n_1+t35_e;
}
else if(t7_l_n>=0){
REAL t35_l_n_1=t3_l_n*t7_h;
REAL t35_h_1=t3_l_n*t7_l_n;
t35_h=t35_h_1+t35_e;
t35_l_n=t35_l_n_1+t35_e;
}
else{
REAL t35_l_n_1=t3_l_n*t7_h;
REAL t3_h_n=-t3_h;
REAL t35_h_1=t3_h_n*t7_l_n;
t35_h=t35_h_1+t35_e;
t35_l_n=t35_l_n_1+t35_e;
}
}
else if(t3_l_n>=0){
if(t7_h<=0){
REAL t35_l_n_1=t3_h*t7_l_n;
REAL t35_h_1=t3_l_n*t7_l_n;
t35_h=t35_h_1+t35_e;
t35_l_n=t35_l_n_1+t35_e;
}
else if(t7_l_n>=0){
REAL t35_l_n_1_1=t3_l_n*t7_h;
REAL t35_l_n_1_2=t7_l_n*t3_h;
REAL t35_l_n_1=max(t35_l_n_1_1,t35_l_n_1_2);
REAL t35_h_1_1=t3_l_n*t7_l_n;
REAL t35_h_1_2=t3_h*t7_h;
REAL t35_h_1=max(t35_h_1_1,t35_h_1_2);
t35_h=t35_h_1+t35_e;
t35_l_n=t35_l_n_1+t35_e;
}
else{
REAL t35_l_n_1=t3_l_n*t7_h;
REAL t35_h_1=t3_h*t7_h;
t35_h=t35_h_1+t35_e;
t35_l_n=t35_l_n_1+t35_e;
}
}
else{
if(t7_h<=0){
REAL t35_l_n_1=t7_l_n*t3_h;
REAL t7_h_n=-t7_h;
REAL t35_h_1=t7_h_n*t3_l_n;
t35_h=t35_h_1+t35_e;
t35_l_n=t35_l_n_1+t35_e;
}
else if(t7_l_n>=0){
REAL t35_l_n_1=t7_l_n*t3_h;
REAL t35_h_1=t3_h*t7_h;
t35_h=t35_h_1+t35_e;
t35_l_n=t35_l_n_1+t35_e;
}
else{
REAL t3_l=-t3_l_n;
REAL t35_l_n_1=t3_l*t7_l_n;
REAL t35_h_1=t3_h*t7_h;
t35_h=t35_h_1+t35_e;
t35_l_n=t35_l_n_1+t35_e;
}
}
REAL t36_e_1=fabs(t36);
REAL t36_e=eps*t36_e_1;
REAL  t36_h, t36_l_n;
if(t1_h<=0){
if(t9_h<=0){
REAL t1_h_n=-t1_h;
REAL t36_l_n_1=t1_h_n*t9_h;
REAL t36_h_1=t1_l_n*t9_l_n;
t36_h=t36_h_1+t36_e;
t36_l_n=t36_l_n_1+t36_e;
}
else if(t9_l_n>=0){
REAL t36_l_n_1=t1_l_n*t9_h;
REAL t36_h_1=t1_l_n*t9_l_n;
t36_h=t36_h_1+t36_e;
t36_l_n=t36_l_n_1+t36_e;
}
else{
REAL t36_l_n_1=t1_l_n*t9_h;
REAL t1_h_n=-t1_h;
REAL t36_h_1=t1_h_n*t9_l_n;
t36_h=t36_h_1+t36_e;
t36_l_n=t36_l_n_1+t36_e;
}
}
else if(t1_l_n>=0){
if(t9_h<=0){
REAL t36_l_n_1=t1_h*t9_l_n;
REAL t36_h_1=t1_l_n*t9_l_n;
t36_h=t36_h_1+t36_e;
t36_l_n=t36_l_n_1+t36_e;
}
else if(t9_l_n>=0){
REAL t36_l_n_1_1=t1_l_n*t9_h;
REAL t36_l_n_1_2=t9_l_n*t1_h;
REAL t36_l_n_1=max(t36_l_n_1_1,t36_l_n_1_2);
REAL t36_h_1_1=t1_l_n*t9_l_n;
REAL t36_h_1_2=t1_h*t9_h;
REAL t36_h_1=max(t36_h_1_1,t36_h_1_2);
t36_h=t36_h_1+t36_e;
t36_l_n=t36_l_n_1+t36_e;
}
else{
REAL t36_l_n_1=t1_l_n*t9_h;
REAL t36_h_1=t1_h*t9_h;
t36_h=t36_h_1+t36_e;
t36_l_n=t36_l_n_1+t36_e;
}
}
else{
if(t9_h<=0){
REAL t36_l_n_1=t9_l_n*t1_h;
REAL t9_h_n=-t9_h;
REAL t36_h_1=t9_h_n*t1_l_n;
t36_h=t36_h_1+t36_e;
t36_l_n=t36_l_n_1+t36_e;
}
else if(t9_l_n>=0){
REAL t36_l_n_1=t9_l_n*t1_h;
REAL t36_h_1=t1_h*t9_h;
t36_h=t36_h_1+t36_e;
t36_l_n=t36_l_n_1+t36_e;
}
else{
REAL t1_l=-t1_l_n;
REAL t36_l_n_1=t1_l*t9_l_n;
REAL t36_h_1=t1_h*t9_h;
t36_h=t36_h_1+t36_e;
t36_l_n=t36_l_n_1+t36_e;
}
}
REAL t38_e_1=fabs(t38);
REAL t38_e=eps*t38_e_1;
REAL t38_h_1=t35_h+t36_l_n;
REAL t38_h=t38_h_1+t38_e;
REAL t38_l_n_1=t36_h+t35_l_n;
REAL t38_l_n=t38_l_n_1+t38_e;
REAL t39_e_1=fabs(t39);
REAL t39_e=eps*t39_e_1;
REAL  t39_h, t39_l_n;
if(t1_h<=0){
if(t1_h<=0){
REAL t1_h_n=-t1_h;
REAL t39_l_n_1=t1_h_n*t1_h;
REAL t39_h_1=t1_l_n*t1_l_n;
t39_h=t39_h_1+t39_e;
t39_l_n=t39_l_n_1+t39_e;
}
else if(t1_l_n>=0){
REAL t39_l_n_1=t1_l_n*t1_h;
REAL t39_h_1=t1_l_n*t1_l_n;
t39_h=t39_h_1+t39_e;
t39_l_n=t39_l_n_1+t39_e;
}
else{
REAL t39_l_n_1=t1_l_n*t1_h;
REAL t1_h_n=-t1_h;
REAL t39_h_1=t1_h_n*t1_l_n;
t39_h=t39_h_1+t39_e;
t39_l_n=t39_l_n_1+t39_e;
}
}
else if(t1_l_n>=0){
if(t1_h<=0){
REAL t39_l_n_1=t1_h*t1_l_n;
REAL t39_h_1=t1_l_n*t1_l_n;
t39_h=t39_h_1+t39_e;
t39_l_n=t39_l_n_1+t39_e;
}
else if(t1_l_n>=0){
REAL t39_l_n_1_1=t1_l_n*t1_h;
REAL t39_l_n_1_2=t1_l_n*t1_h;
REAL t39_l_n_1=max(t39_l_n_1_1,t39_l_n_1_2);
REAL t39_h_1_1=t1_l_n*t1_l_n;
REAL t39_h_1_2=t1_h*t1_h;
REAL t39_h_1=max(t39_h_1_1,t39_h_1_2);
t39_h=t39_h_1+t39_e;
t39_l_n=t39_l_n_1+t39_e;
}
else{
REAL t39_l_n_1=t1_l_n*t1_h;
REAL t39_h_1=t1_h*t1_h;
t39_h=t39_h_1+t39_e;
t39_l_n=t39_l_n_1+t39_e;
}
}
else{
if(t1_h<=0){
REAL t39_l_n_1=t1_l_n*t1_h;
REAL t1_h_n=-t1_h;
REAL t39_h_1=t1_h_n*t1_l_n;
t39_h=t39_h_1+t39_e;
t39_l_n=t39_l_n_1+t39_e;
}
else if(t1_l_n>=0){
REAL t39_l_n_1=t1_l_n*t1_h;
REAL t39_h_1=t1_h*t1_h;
t39_h=t39_h_1+t39_e;
t39_l_n=t39_l_n_1+t39_e;
}
else{
REAL t1_l=-t1_l_n;
REAL t39_l_n_1=t1_l*t1_l_n;
REAL t39_h_1=t1_h*t1_h;
t39_h=t39_h_1+t39_e;
t39_l_n=t39_l_n_1+t39_e;
}
}
REAL t40_e_1=fabs(t40);
REAL t40_e=eps*t40_e_1;
REAL  t40_h, t40_l_n;
if(t3_h<=0){
if(t3_h<=0){
REAL t3_h_n=-t3_h;
REAL t40_l_n_1=t3_h_n*t3_h;
REAL t40_h_1=t3_l_n*t3_l_n;
t40_h=t40_h_1+t40_e;
t40_l_n=t40_l_n_1+t40_e;
}
else if(t3_l_n>=0){
REAL t40_l_n_1=t3_l_n*t3_h;
REAL t40_h_1=t3_l_n*t3_l_n;
t40_h=t40_h_1+t40_e;
t40_l_n=t40_l_n_1+t40_e;
}
else{
REAL t40_l_n_1=t3_l_n*t3_h;
REAL t3_h_n=-t3_h;
REAL t40_h_1=t3_h_n*t3_l_n;
t40_h=t40_h_1+t40_e;
t40_l_n=t40_l_n_1+t40_e;
}
}
else if(t3_l_n>=0){
if(t3_h<=0){
REAL t40_l_n_1=t3_h*t3_l_n;
REAL t40_h_1=t3_l_n*t3_l_n;
t40_h=t40_h_1+t40_e;
t40_l_n=t40_l_n_1+t40_e;
}
else if(t3_l_n>=0){
REAL t40_l_n_1_1=t3_l_n*t3_h;
REAL t40_l_n_1_2=t3_l_n*t3_h;
REAL t40_l_n_1=max(t40_l_n_1_1,t40_l_n_1_2);
REAL t40_h_1_1=t3_l_n*t3_l_n;
REAL t40_h_1_2=t3_h*t3_h;
REAL t40_h_1=max(t40_h_1_1,t40_h_1_2);
t40_h=t40_h_1+t40_e;
t40_l_n=t40_l_n_1+t40_e;
}
else{
REAL t40_l_n_1=t3_l_n*t3_h;
REAL t40_h_1=t3_h*t3_h;
t40_h=t40_h_1+t40_e;
t40_l_n=t40_l_n_1+t40_e;
}
}
else{
if(t3_h<=0){
REAL t40_l_n_1=t3_l_n*t3_h;
REAL t3_h_n=-t3_h;
REAL t40_h_1=t3_h_n*t3_l_n;
t40_h=t40_h_1+t40_e;
t40_l_n=t40_l_n_1+t40_e;
}
else if(t3_l_n>=0){
REAL t40_l_n_1=t3_l_n*t3_h;
REAL t40_h_1=t3_h*t3_h;
t40_h=t40_h_1+t40_e;
t40_l_n=t40_l_n_1+t40_e;
}
else{
REAL t3_l=-t3_l_n;
REAL t40_l_n_1=t3_l*t3_l_n;
REAL t40_h_1=t3_h*t3_h;
t40_h=t40_h_1+t40_e;
t40_l_n=t40_l_n_1+t40_e;
}
}
REAL t41_e_1=fabs(t41);
REAL t41_e=eps*t41_e_1;
REAL  t41_h, t41_l_n;
if(t5_h<=0){
if(t5_h<=0){
REAL t5_h_n=-t5_h;
REAL t41_l_n_1=t5_h_n*t5_h;
REAL t41_h_1=t5_l_n*t5_l_n;
t41_h=t41_h_1+t41_e;
t41_l_n=t41_l_n_1+t41_e;
}
else if(t5_l_n>=0){
REAL t41_l_n_1=t5_l_n*t5_h;
REAL t41_h_1=t5_l_n*t5_l_n;
t41_h=t41_h_1+t41_e;
t41_l_n=t41_l_n_1+t41_e;
}
else{
REAL t41_l_n_1=t5_l_n*t5_h;
REAL t5_h_n=-t5_h;
REAL t41_h_1=t5_h_n*t5_l_n;
t41_h=t41_h_1+t41_e;
t41_l_n=t41_l_n_1+t41_e;
}
}
else if(t5_l_n>=0){
if(t5_h<=0){
REAL t41_l_n_1=t5_h*t5_l_n;
REAL t41_h_1=t5_l_n*t5_l_n;
t41_h=t41_h_1+t41_e;
t41_l_n=t41_l_n_1+t41_e;
}
else if(t5_l_n>=0){
REAL t41_l_n_1_1=t5_l_n*t5_h;
REAL t41_l_n_1_2=t5_l_n*t5_h;
REAL t41_l_n_1=max(t41_l_n_1_1,t41_l_n_1_2);
REAL t41_h_1_1=t5_l_n*t5_l_n;
REAL t41_h_1_2=t5_h*t5_h;
REAL t41_h_1=max(t41_h_1_1,t41_h_1_2);
t41_h=t41_h_1+t41_e;
t41_l_n=t41_l_n_1+t41_e;
}
else{
REAL t41_l_n_1=t5_l_n*t5_h;
REAL t41_h_1=t5_h*t5_h;
t41_h=t41_h_1+t41_e;
t41_l_n=t41_l_n_1+t41_e;
}
}
else{
if(t5_h<=0){
REAL t41_l_n_1=t5_l_n*t5_h;
REAL t5_h_n=-t5_h;
REAL t41_h_1=t5_h_n*t5_l_n;
t41_h=t41_h_1+t41_e;
t41_l_n=t41_l_n_1+t41_e;
}
else if(t5_l_n>=0){
REAL t41_l_n_1=t5_l_n*t5_h;
REAL t41_h_1=t5_h*t5_h;
t41_h=t41_h_1+t41_e;
t41_l_n=t41_l_n_1+t41_e;
}
else{
REAL t5_l=-t5_l_n;
REAL t41_l_n_1=t5_l*t5_l_n;
REAL t41_h_1=t5_h*t5_h;
t41_h=t41_h_1+t41_e;
t41_l_n=t41_l_n_1+t41_e;
}
}
REAL t42_e_1=fabs(t42);
REAL t42_e=eps*t42_e_1;
REAL t42_h_1=t39_h+t40_h;
REAL t42_h=t42_h_1+t42_e;
REAL t42_l_n_1=t39_l_n+t40_l_n;
REAL t42_l_n=t42_l_n_1+t42_e;
REAL t43_e_1=fabs(t43);
REAL t43_e=eps*t43_e_1;
REAL t43_h_1=t41_h+t42_h;
REAL t43_h=t43_h_1+t43_e;
REAL t43_l_n_1=t41_l_n+t42_l_n;
REAL t43_l_n=t43_l_n_1+t43_e;
REAL t44_e_1=fabs(t44);
REAL t44_e=eps*t44_e_1;
REAL  t44_h, t44_l_n;
if(t14_h<=0){
if(t6_h<=0){
REAL t14_h_n=-t14_h;
REAL t44_l_n_1=t14_h_n*t6_h;
REAL t44_h_1=t14_l_n*t6_l_n;
t44_h=t44_h_1+t44_e;
t44_l_n=t44_l_n_1+t44_e;
}
else if(t6_l_n>=0){
REAL t44_l_n_1=t14_l_n*t6_h;
REAL t44_h_1=t14_l_n*t6_l_n;
t44_h=t44_h_1+t44_e;
t44_l_n=t44_l_n_1+t44_e;
}
else{
REAL t44_l_n_1=t14_l_n*t6_h;
REAL t14_h_n=-t14_h;
REAL t44_h_1=t14_h_n*t6_l_n;
t44_h=t44_h_1+t44_e;
t44_l_n=t44_l_n_1+t44_e;
}
}
else if(t14_l_n>=0){
if(t6_h<=0){
REAL t44_l_n_1=t14_h*t6_l_n;
REAL t44_h_1=t14_l_n*t6_l_n;
t44_h=t44_h_1+t44_e;
t44_l_n=t44_l_n_1+t44_e;
}
else if(t6_l_n>=0){
REAL t44_l_n_1_1=t14_l_n*t6_h;
REAL t44_l_n_1_2=t6_l_n*t14_h;
REAL t44_l_n_1=max(t44_l_n_1_1,t44_l_n_1_2);
REAL t44_h_1_1=t14_l_n*t6_l_n;
REAL t44_h_1_2=t14_h*t6_h;
REAL t44_h_1=max(t44_h_1_1,t44_h_1_2);
t44_h=t44_h_1+t44_e;
t44_l_n=t44_l_n_1+t44_e;
}
else{
REAL t44_l_n_1=t14_l_n*t6_h;
REAL t44_h_1=t14_h*t6_h;
t44_h=t44_h_1+t44_e;
t44_l_n=t44_l_n_1+t44_e;
}
}
else{
if(t6_h<=0){
REAL t44_l_n_1=t6_l_n*t14_h;
REAL t6_h_n=-t6_h;
REAL t44_h_1=t6_h_n*t14_l_n;
t44_h=t44_h_1+t44_e;
t44_l_n=t44_l_n_1+t44_e;
}
else if(t6_l_n>=0){
REAL t44_l_n_1=t6_l_n*t14_h;
REAL t44_h_1=t14_h*t6_h;
t44_h=t44_h_1+t44_e;
t44_l_n=t44_l_n_1+t44_e;
}
else{
REAL t14_l=-t14_l_n;
REAL t44_l_n_1=t14_l*t6_l_n;
REAL t44_h_1=t14_h*t6_h;
t44_h=t44_h_1+t44_e;
t44_l_n=t44_l_n_1+t44_e;
}
}
REAL t45_e_1=fabs(t45);
REAL t45_e=eps*t45_e_1;
REAL  t45_h, t45_l_n;
if(t15_h<=0){
if(t21_h<=0){
REAL t15_h_n=-t15_h;
REAL t45_l_n_1=t15_h_n*t21_h;
REAL t45_h_1=t15_l_n*t21_l_n;
t45_h=t45_h_1+t45_e;
t45_l_n=t45_l_n_1+t45_e;
}
else if(t21_l_n>=0){
REAL t45_l_n_1=t15_l_n*t21_h;
REAL t45_h_1=t15_l_n*t21_l_n;
t45_h=t45_h_1+t45_e;
t45_l_n=t45_l_n_1+t45_e;
}
else{
REAL t45_l_n_1=t15_l_n*t21_h;
REAL t15_h_n=-t15_h;
REAL t45_h_1=t15_h_n*t21_l_n;
t45_h=t45_h_1+t45_e;
t45_l_n=t45_l_n_1+t45_e;
}
}
else if(t15_l_n>=0){
if(t21_h<=0){
REAL t45_l_n_1=t15_h*t21_l_n;
REAL t45_h_1=t15_l_n*t21_l_n;
t45_h=t45_h_1+t45_e;
t45_l_n=t45_l_n_1+t45_e;
}
else if(t21_l_n>=0){
REAL t45_l_n_1_1=t15_l_n*t21_h;
REAL t45_l_n_1_2=t21_l_n*t15_h;
REAL t45_l_n_1=max(t45_l_n_1_1,t45_l_n_1_2);
REAL t45_h_1_1=t15_l_n*t21_l_n;
REAL t45_h_1_2=t15_h*t21_h;
REAL t45_h_1=max(t45_h_1_1,t45_h_1_2);
t45_h=t45_h_1+t45_e;
t45_l_n=t45_l_n_1+t45_e;
}
else{
REAL t45_l_n_1=t15_l_n*t21_h;
REAL t45_h_1=t15_h*t21_h;
t45_h=t45_h_1+t45_e;
t45_l_n=t45_l_n_1+t45_e;
}
}
else{
if(t21_h<=0){
REAL t45_l_n_1=t21_l_n*t15_h;
REAL t21_h_n=-t21_h;
REAL t45_h_1=t21_h_n*t15_l_n;
t45_h=t45_h_1+t45_e;
t45_l_n=t45_l_n_1+t45_e;
}
else if(t21_l_n>=0){
REAL t45_l_n_1=t21_l_n*t15_h;
REAL t45_h_1=t15_h*t21_h;
t45_h=t45_h_1+t45_e;
t45_l_n=t45_l_n_1+t45_e;
}
else{
REAL t15_l=-t15_l_n;
REAL t45_l_n_1=t15_l*t21_l_n;
REAL t45_h_1=t15_h*t21_h;
t45_h=t45_h_1+t45_e;
t45_l_n=t45_l_n_1+t45_e;
}
}
REAL t46_e_1=fabs(t46);
REAL t46_e=eps*t46_e_1;
REAL  t46_h, t46_l_n;
if(t22_h<=0){
if(t26_h<=0){
REAL t22_h_n=-t22_h;
REAL t46_l_n_1=t22_h_n*t26_h;
REAL t46_h_1=t22_l_n*t26_l_n;
t46_h=t46_h_1+t46_e;
t46_l_n=t46_l_n_1+t46_e;
}
else if(t26_l_n>=0){
REAL t46_l_n_1=t22_l_n*t26_h;
REAL t46_h_1=t22_l_n*t26_l_n;
t46_h=t46_h_1+t46_e;
t46_l_n=t46_l_n_1+t46_e;
}
else{
REAL t46_l_n_1=t22_l_n*t26_h;
REAL t22_h_n=-t22_h;
REAL t46_h_1=t22_h_n*t26_l_n;
t46_h=t46_h_1+t46_e;
t46_l_n=t46_l_n_1+t46_e;
}
}
else if(t22_l_n>=0){
if(t26_h<=0){
REAL t46_l_n_1=t22_h*t26_l_n;
REAL t46_h_1=t22_l_n*t26_l_n;
t46_h=t46_h_1+t46_e;
t46_l_n=t46_l_n_1+t46_e;
}
else if(t26_l_n>=0){
REAL t46_l_n_1_1=t22_l_n*t26_h;
REAL t46_l_n_1_2=t26_l_n*t22_h;
REAL t46_l_n_1=max(t46_l_n_1_1,t46_l_n_1_2);
REAL t46_h_1_1=t22_l_n*t26_l_n;
REAL t46_h_1_2=t22_h*t26_h;
REAL t46_h_1=max(t46_h_1_1,t46_h_1_2);
t46_h=t46_h_1+t46_e;
t46_l_n=t46_l_n_1+t46_e;
}
else{
REAL t46_l_n_1=t22_l_n*t26_h;
REAL t46_h_1=t22_h*t26_h;
t46_h=t46_h_1+t46_e;
t46_l_n=t46_l_n_1+t46_e;
}
}
else{
if(t26_h<=0){
REAL t46_l_n_1=t26_l_n*t22_h;
REAL t26_h_n=-t26_h;
REAL t46_h_1=t26_h_n*t22_l_n;
t46_h=t46_h_1+t46_e;
t46_l_n=t46_l_n_1+t46_e;
}
else if(t26_l_n>=0){
REAL t46_l_n_1=t26_l_n*t22_h;
REAL t46_h_1=t22_h*t26_h;
t46_h=t46_h_1+t46_e;
t46_l_n=t46_l_n_1+t46_e;
}
else{
REAL t22_l=-t22_l_n;
REAL t46_l_n_1=t22_l*t26_l_n;
REAL t46_h_1=t22_h*t26_h;
t46_h=t46_h_1+t46_e;
t46_l_n=t46_l_n_1+t46_e;
}
}
REAL t48_e_1=fabs(t48);
REAL t48_e=eps*t48_e_1;
REAL t48_h_1=t44_h+t45_h;
REAL t48_h=t48_h_1+t48_e;
REAL t48_l_n_1=t44_l_n+t45_l_n;
REAL t48_l_n=t48_l_n_1+t48_e;
REAL t49_e_1=fabs(t49);
REAL t49_e=eps*t49_e_1;
REAL t49_h_1=t48_h+t46_l_n;
REAL t49_h=t49_h_1+t49_e;
REAL t49_l_n_1=t46_h+t48_l_n;
REAL t49_l_n=t49_l_n_1+t49_e;
REAL t50_e_1=fabs(t50);
REAL t50_e=eps*t50_e_1;
REAL  t50_h, t50_l_n;
if(t43_h<=0){
if(t49_h<=0){
REAL t43_h_n=-t43_h;
REAL t50_l_n_1=t43_h_n*t49_h;
REAL t50_h_1=t43_l_n*t49_l_n;
t50_h=t50_h_1+t50_e;
t50_l_n=t50_l_n_1+t50_e;
}
else if(t49_l_n>=0){
REAL t50_l_n_1=t43_l_n*t49_h;
REAL t50_h_1=t43_l_n*t49_l_n;
t50_h=t50_h_1+t50_e;
t50_l_n=t50_l_n_1+t50_e;
}
else{
REAL t50_l_n_1=t43_l_n*t49_h;
REAL t43_h_n=-t43_h;
REAL t50_h_1=t43_h_n*t49_l_n;
t50_h=t50_h_1+t50_e;
t50_l_n=t50_l_n_1+t50_e;
}
}
else if(t43_l_n>=0){
if(t49_h<=0){
REAL t50_l_n_1=t43_h*t49_l_n;
REAL t50_h_1=t43_l_n*t49_l_n;
t50_h=t50_h_1+t50_e;
t50_l_n=t50_l_n_1+t50_e;
}
else if(t49_l_n>=0){
REAL t50_l_n_1_1=t43_l_n*t49_h;
REAL t50_l_n_1_2=t49_l_n*t43_h;
REAL t50_l_n_1=max(t50_l_n_1_1,t50_l_n_1_2);
REAL t50_h_1_1=t43_l_n*t49_l_n;
REAL t50_h_1_2=t43_h*t49_h;
REAL t50_h_1=max(t50_h_1_1,t50_h_1_2);
t50_h=t50_h_1+t50_e;
t50_l_n=t50_l_n_1+t50_e;
}
else{
REAL t50_l_n_1=t43_l_n*t49_h;
REAL t50_h_1=t43_h*t49_h;
t50_h=t50_h_1+t50_e;
t50_l_n=t50_l_n_1+t50_e;
}
}
else{
if(t49_h<=0){
REAL t50_l_n_1=t49_l_n*t43_h;
REAL t49_h_n=-t49_h;
REAL t50_h_1=t49_h_n*t43_l_n;
t50_h=t50_h_1+t50_e;
t50_l_n=t50_l_n_1+t50_e;
}
else if(t49_l_n>=0){
REAL t50_l_n_1=t49_l_n*t43_h;
REAL t50_h_1=t43_h*t49_h;
t50_h=t50_h_1+t50_e;
t50_l_n=t50_l_n_1+t50_e;
}
else{
REAL t43_l=-t43_l_n;
REAL t50_l_n_1=t43_l*t49_l_n;
REAL t50_h_1=t43_h*t49_h;
t50_h=t50_h_1+t50_e;
t50_l_n=t50_l_n_1+t50_e;
}
}
REAL t51_e_1=fabs(t51);
REAL t51_e=eps*t51_e_1;
REAL  t51_h, t51_l_n;
if(t22_h<=0){
if(t22_h<=0){
REAL t22_h_n=-t22_h;
REAL t51_l_n_1=t22_h_n*t22_h;
REAL t51_h_1=t22_l_n*t22_l_n;
t51_h=t51_h_1+t51_e;
t51_l_n=t51_l_n_1+t51_e;
}
else if(t22_l_n>=0){
REAL t51_l_n_1=t22_l_n*t22_h;
REAL t51_h_1=t22_l_n*t22_l_n;
t51_h=t51_h_1+t51_e;
t51_l_n=t51_l_n_1+t51_e;
}
else{
REAL t51_l_n_1=t22_l_n*t22_h;
REAL t22_h_n=-t22_h;
REAL t51_h_1=t22_h_n*t22_l_n;
t51_h=t51_h_1+t51_e;
t51_l_n=t51_l_n_1+t51_e;
}
}
else if(t22_l_n>=0){
if(t22_h<=0){
REAL t51_l_n_1=t22_h*t22_l_n;
REAL t51_h_1=t22_l_n*t22_l_n;
t51_h=t51_h_1+t51_e;
t51_l_n=t51_l_n_1+t51_e;
}
else if(t22_l_n>=0){
REAL t51_l_n_1_1=t22_l_n*t22_h;
REAL t51_l_n_1_2=t22_l_n*t22_h;
REAL t51_l_n_1=max(t51_l_n_1_1,t51_l_n_1_2);
REAL t51_h_1_1=t22_l_n*t22_l_n;
REAL t51_h_1_2=t22_h*t22_h;
REAL t51_h_1=max(t51_h_1_1,t51_h_1_2);
t51_h=t51_h_1+t51_e;
t51_l_n=t51_l_n_1+t51_e;
}
else{
REAL t51_l_n_1=t22_l_n*t22_h;
REAL t51_h_1=t22_h*t22_h;
t51_h=t51_h_1+t51_e;
t51_l_n=t51_l_n_1+t51_e;
}
}
else{
if(t22_h<=0){
REAL t51_l_n_1=t22_l_n*t22_h;
REAL t22_h_n=-t22_h;
REAL t51_h_1=t22_h_n*t22_l_n;
t51_h=t51_h_1+t51_e;
t51_l_n=t51_l_n_1+t51_e;
}
else if(t22_l_n>=0){
REAL t51_l_n_1=t22_l_n*t22_h;
REAL t51_h_1=t22_h*t22_h;
t51_h=t51_h_1+t51_e;
t51_l_n=t51_l_n_1+t51_e;
}
else{
REAL t22_l=-t22_l_n;
REAL t51_l_n_1=t22_l*t22_l_n;
REAL t51_h_1=t22_h*t22_h;
t51_h=t51_h_1+t51_e;
t51_l_n=t51_l_n_1+t51_e;
}
}
REAL t52_e_1=fabs(t52);
REAL t52_e=eps*t52_e_1;
REAL  t52_h, t52_l_n;
if(t7_h<=0){
if(t7_h<=0){
REAL t7_h_n=-t7_h;
REAL t52_l_n_1=t7_h_n*t7_h;
REAL t52_h_1=t7_l_n*t7_l_n;
t52_h=t52_h_1+t52_e;
t52_l_n=t52_l_n_1+t52_e;
}
else if(t7_l_n>=0){
REAL t52_l_n_1=t7_l_n*t7_h;
REAL t52_h_1=t7_l_n*t7_l_n;
t52_h=t52_h_1+t52_e;
t52_l_n=t52_l_n_1+t52_e;
}
else{
REAL t52_l_n_1=t7_l_n*t7_h;
REAL t7_h_n=-t7_h;
REAL t52_h_1=t7_h_n*t7_l_n;
t52_h=t52_h_1+t52_e;
t52_l_n=t52_l_n_1+t52_e;
}
}
else if(t7_l_n>=0){
if(t7_h<=0){
REAL t52_l_n_1=t7_h*t7_l_n;
REAL t52_h_1=t7_l_n*t7_l_n;
t52_h=t52_h_1+t52_e;
t52_l_n=t52_l_n_1+t52_e;
}
else if(t7_l_n>=0){
REAL t52_l_n_1_1=t7_l_n*t7_h;
REAL t52_l_n_1_2=t7_l_n*t7_h;
REAL t52_l_n_1=max(t52_l_n_1_1,t52_l_n_1_2);
REAL t52_h_1_1=t7_l_n*t7_l_n;
REAL t52_h_1_2=t7_h*t7_h;
REAL t52_h_1=max(t52_h_1_1,t52_h_1_2);
t52_h=t52_h_1+t52_e;
t52_l_n=t52_l_n_1+t52_e;
}
else{
REAL t52_l_n_1=t7_l_n*t7_h;
REAL t52_h_1=t7_h*t7_h;
t52_h=t52_h_1+t52_e;
t52_l_n=t52_l_n_1+t52_e;
}
}
else{
if(t7_h<=0){
REAL t52_l_n_1=t7_l_n*t7_h;
REAL t7_h_n=-t7_h;
REAL t52_h_1=t7_h_n*t7_l_n;
t52_h=t52_h_1+t52_e;
t52_l_n=t52_l_n_1+t52_e;
}
else if(t7_l_n>=0){
REAL t52_l_n_1=t7_l_n*t7_h;
REAL t52_h_1=t7_h*t7_h;
t52_h=t52_h_1+t52_e;
t52_l_n=t52_l_n_1+t52_e;
}
else{
REAL t7_l=-t7_l_n;
REAL t52_l_n_1=t7_l*t7_l_n;
REAL t52_h_1=t7_h*t7_h;
t52_h=t52_h_1+t52_e;
t52_l_n=t52_l_n_1+t52_e;
}
}
REAL t53_e_1=fabs(t53);
REAL t53_e=eps*t53_e_1;
REAL  t53_h, t53_l_n;
if(t9_h<=0){
if(t9_h<=0){
REAL t9_h_n=-t9_h;
REAL t53_l_n_1=t9_h_n*t9_h;
REAL t53_h_1=t9_l_n*t9_l_n;
t53_h=t53_h_1+t53_e;
t53_l_n=t53_l_n_1+t53_e;
}
else if(t9_l_n>=0){
REAL t53_l_n_1=t9_l_n*t9_h;
REAL t53_h_1=t9_l_n*t9_l_n;
t53_h=t53_h_1+t53_e;
t53_l_n=t53_l_n_1+t53_e;
}
else{
REAL t53_l_n_1=t9_l_n*t9_h;
REAL t9_h_n=-t9_h;
REAL t53_h_1=t9_h_n*t9_l_n;
t53_h=t53_h_1+t53_e;
t53_l_n=t53_l_n_1+t53_e;
}
}
else if(t9_l_n>=0){
if(t9_h<=0){
REAL t53_l_n_1=t9_h*t9_l_n;
REAL t53_h_1=t9_l_n*t9_l_n;
t53_h=t53_h_1+t53_e;
t53_l_n=t53_l_n_1+t53_e;
}
else if(t9_l_n>=0){
REAL t53_l_n_1_1=t9_l_n*t9_h;
REAL t53_l_n_1_2=t9_l_n*t9_h;
REAL t53_l_n_1=max(t53_l_n_1_1,t53_l_n_1_2);
REAL t53_h_1_1=t9_l_n*t9_l_n;
REAL t53_h_1_2=t9_h*t9_h;
REAL t53_h_1=max(t53_h_1_1,t53_h_1_2);
t53_h=t53_h_1+t53_e;
t53_l_n=t53_l_n_1+t53_e;
}
else{
REAL t53_l_n_1=t9_l_n*t9_h;
REAL t53_h_1=t9_h*t9_h;
t53_h=t53_h_1+t53_e;
t53_l_n=t53_l_n_1+t53_e;
}
}
else{
if(t9_h<=0){
REAL t53_l_n_1=t9_l_n*t9_h;
REAL t9_h_n=-t9_h;
REAL t53_h_1=t9_h_n*t9_l_n;
t53_h=t53_h_1+t53_e;
t53_l_n=t53_l_n_1+t53_e;
}
else if(t9_l_n>=0){
REAL t53_l_n_1=t9_l_n*t9_h;
REAL t53_h_1=t9_h*t9_h;
t53_h=t53_h_1+t53_e;
t53_l_n=t53_l_n_1+t53_e;
}
else{
REAL t9_l=-t9_l_n;
REAL t53_l_n_1=t9_l*t9_l_n;
REAL t53_h_1=t9_h*t9_h;
t53_h=t53_h_1+t53_e;
t53_l_n=t53_l_n_1+t53_e;
}
}
REAL t54_e_1=fabs(t54);
REAL t54_e=eps*t54_e_1;
REAL t54_h_1=t51_h+t52_h;
REAL t54_h=t54_h_1+t54_e;
REAL t54_l_n_1=t51_l_n+t52_l_n;
REAL t54_l_n=t54_l_n_1+t54_e;
REAL t55_e_1=fabs(t55);
REAL t55_e=eps*t55_e_1;
REAL t55_h_1=t53_h+t54_h;
REAL t55_h=t55_h_1+t55_e;
REAL t55_l_n_1=t53_l_n+t54_l_n;
REAL t55_l_n=t55_l_n_1+t55_e;
REAL t56_e_1=fabs(t56);
REAL t56_e=eps*t56_e_1;
REAL  t56_h, t56_l_n;
if(t26_h<=0){
if(t5_h<=0){
REAL t26_h_n=-t26_h;
REAL t56_l_n_1=t26_h_n*t5_h;
REAL t56_h_1=t26_l_n*t5_l_n;
t56_h=t56_h_1+t56_e;
t56_l_n=t56_l_n_1+t56_e;
}
else if(t5_l_n>=0){
REAL t56_l_n_1=t26_l_n*t5_h;
REAL t56_h_1=t26_l_n*t5_l_n;
t56_h=t56_h_1+t56_e;
t56_l_n=t56_l_n_1+t56_e;
}
else{
REAL t56_l_n_1=t26_l_n*t5_h;
REAL t26_h_n=-t26_h;
REAL t56_h_1=t26_h_n*t5_l_n;
t56_h=t56_h_1+t56_e;
t56_l_n=t56_l_n_1+t56_e;
}
}
else if(t26_l_n>=0){
if(t5_h<=0){
REAL t56_l_n_1=t26_h*t5_l_n;
REAL t56_h_1=t26_l_n*t5_l_n;
t56_h=t56_h_1+t56_e;
t56_l_n=t56_l_n_1+t56_e;
}
else if(t5_l_n>=0){
REAL t56_l_n_1_1=t26_l_n*t5_h;
REAL t56_l_n_1_2=t5_l_n*t26_h;
REAL t56_l_n_1=max(t56_l_n_1_1,t56_l_n_1_2);
REAL t56_h_1_1=t26_l_n*t5_l_n;
REAL t56_h_1_2=t26_h*t5_h;
REAL t56_h_1=max(t56_h_1_1,t56_h_1_2);
t56_h=t56_h_1+t56_e;
t56_l_n=t56_l_n_1+t56_e;
}
else{
REAL t56_l_n_1=t26_l_n*t5_h;
REAL t56_h_1=t26_h*t5_h;
t56_h=t56_h_1+t56_e;
t56_l_n=t56_l_n_1+t56_e;
}
}
else{
if(t5_h<=0){
REAL t56_l_n_1=t5_l_n*t26_h;
REAL t5_h_n=-t5_h;
REAL t56_h_1=t5_h_n*t26_l_n;
t56_h=t56_h_1+t56_e;
t56_l_n=t56_l_n_1+t56_e;
}
else if(t5_l_n>=0){
REAL t56_l_n_1=t5_l_n*t26_h;
REAL t56_h_1=t26_h*t5_h;
t56_h=t56_h_1+t56_e;
t56_l_n=t56_l_n_1+t56_e;
}
else{
REAL t26_l=-t26_l_n;
REAL t56_l_n_1=t26_l*t5_l_n;
REAL t56_h_1=t26_h*t5_h;
t56_h=t56_h_1+t56_e;
t56_l_n=t56_l_n_1+t56_e;
}
}
REAL t57_e_1=fabs(t57);
REAL t57_e=eps*t57_e_1;
REAL  t57_h, t57_l_n;
if(t30_h<=0){
if(t6_h<=0){
REAL t30_h_n=-t30_h;
REAL t57_l_n_1=t30_h_n*t6_h;
REAL t57_h_1=t30_l_n*t6_l_n;
t57_h=t57_h_1+t57_e;
t57_l_n=t57_l_n_1+t57_e;
}
else if(t6_l_n>=0){
REAL t57_l_n_1=t30_l_n*t6_h;
REAL t57_h_1=t30_l_n*t6_l_n;
t57_h=t57_h_1+t57_e;
t57_l_n=t57_l_n_1+t57_e;
}
else{
REAL t57_l_n_1=t30_l_n*t6_h;
REAL t30_h_n=-t30_h;
REAL t57_h_1=t30_h_n*t6_l_n;
t57_h=t57_h_1+t57_e;
t57_l_n=t57_l_n_1+t57_e;
}
}
else if(t30_l_n>=0){
if(t6_h<=0){
REAL t57_l_n_1=t30_h*t6_l_n;
REAL t57_h_1=t30_l_n*t6_l_n;
t57_h=t57_h_1+t57_e;
t57_l_n=t57_l_n_1+t57_e;
}
else if(t6_l_n>=0){
REAL t57_l_n_1_1=t30_l_n*t6_h;
REAL t57_l_n_1_2=t6_l_n*t30_h;
REAL t57_l_n_1=max(t57_l_n_1_1,t57_l_n_1_2);
REAL t57_h_1_1=t30_l_n*t6_l_n;
REAL t57_h_1_2=t30_h*t6_h;
REAL t57_h_1=max(t57_h_1_1,t57_h_1_2);
t57_h=t57_h_1+t57_e;
t57_l_n=t57_l_n_1+t57_e;
}
else{
REAL t57_l_n_1=t30_l_n*t6_h;
REAL t57_h_1=t30_h*t6_h;
t57_h=t57_h_1+t57_e;
t57_l_n=t57_l_n_1+t57_e;
}
}
else{
if(t6_h<=0){
REAL t57_l_n_1=t6_l_n*t30_h;
REAL t6_h_n=-t6_h;
REAL t57_h_1=t6_h_n*t30_l_n;
t57_h=t57_h_1+t57_e;
t57_l_n=t57_l_n_1+t57_e;
}
else if(t6_l_n>=0){
REAL t57_l_n_1=t6_l_n*t30_h;
REAL t57_h_1=t30_h*t6_h;
t57_h=t57_h_1+t57_e;
t57_l_n=t57_l_n_1+t57_e;
}
else{
REAL t30_l=-t30_l_n;
REAL t57_l_n_1=t30_l*t6_l_n;
REAL t57_h_1=t30_h*t6_h;
t57_h=t57_h_1+t57_e;
t57_l_n=t57_l_n_1+t57_e;
}
}
REAL t58_e_1=fabs(t58);
REAL t58_e=eps*t58_e_1;
REAL  t58_h, t58_l_n;
if(t15_h<=0){
if(t34_h<=0){
REAL t15_h_n=-t15_h;
REAL t58_l_n_1=t15_h_n*t34_h;
REAL t58_h_1=t15_l_n*t34_l_n;
t58_h=t58_h_1+t58_e;
t58_l_n=t58_l_n_1+t58_e;
}
else if(t34_l_n>=0){
REAL t58_l_n_1=t15_l_n*t34_h;
REAL t58_h_1=t15_l_n*t34_l_n;
t58_h=t58_h_1+t58_e;
t58_l_n=t58_l_n_1+t58_e;
}
else{
REAL t58_l_n_1=t15_l_n*t34_h;
REAL t15_h_n=-t15_h;
REAL t58_h_1=t15_h_n*t34_l_n;
t58_h=t58_h_1+t58_e;
t58_l_n=t58_l_n_1+t58_e;
}
}
else if(t15_l_n>=0){
if(t34_h<=0){
REAL t58_l_n_1=t15_h*t34_l_n;
REAL t58_h_1=t15_l_n*t34_l_n;
t58_h=t58_h_1+t58_e;
t58_l_n=t58_l_n_1+t58_e;
}
else if(t34_l_n>=0){
REAL t58_l_n_1_1=t15_l_n*t34_h;
REAL t58_l_n_1_2=t34_l_n*t15_h;
REAL t58_l_n_1=max(t58_l_n_1_1,t58_l_n_1_2);
REAL t58_h_1_1=t15_l_n*t34_l_n;
REAL t58_h_1_2=t15_h*t34_h;
REAL t58_h_1=max(t58_h_1_1,t58_h_1_2);
t58_h=t58_h_1+t58_e;
t58_l_n=t58_l_n_1+t58_e;
}
else{
REAL t58_l_n_1=t15_l_n*t34_h;
REAL t58_h_1=t15_h*t34_h;
t58_h=t58_h_1+t58_e;
t58_l_n=t58_l_n_1+t58_e;
}
}
else{
if(t34_h<=0){
REAL t58_l_n_1=t34_l_n*t15_h;
REAL t34_h_n=-t34_h;
REAL t58_h_1=t34_h_n*t15_l_n;
t58_h=t58_h_1+t58_e;
t58_l_n=t58_l_n_1+t58_e;
}
else if(t34_l_n>=0){
REAL t58_l_n_1=t34_l_n*t15_h;
REAL t58_h_1=t15_h*t34_h;
t58_h=t58_h_1+t58_e;
t58_l_n=t58_l_n_1+t58_e;
}
else{
REAL t15_l=-t15_l_n;
REAL t58_l_n_1=t15_l*t34_l_n;
REAL t58_h_1=t15_h*t34_h;
t58_h=t58_h_1+t58_e;
t58_l_n=t58_l_n_1+t58_e;
}
}
REAL t60_e_1=fabs(t60);
REAL t60_e=eps*t60_e_1;
REAL t60_h_1=t56_h+t57_h;
REAL t60_h=t60_h_1+t60_e;
REAL t60_l_n_1=t56_l_n+t57_l_n;
REAL t60_l_n=t60_l_n_1+t60_e;
REAL t61_e_1=fabs(t61);
REAL t61_e=eps*t61_e_1;
REAL t61_h_1=t60_h+t58_l_n;
REAL t61_h=t61_h_1+t61_e;
REAL t61_l_n_1=t58_h+t60_l_n;
REAL t61_l_n=t61_l_n_1+t61_e;
REAL t62_e_1=fabs(t62);
REAL t62_e=eps*t62_e_1;
REAL  t62_h, t62_l_n;
if(t55_h<=0){
if(t61_h<=0){
REAL t55_h_n=-t55_h;
REAL t62_l_n_1=t55_h_n*t61_h;
REAL t62_h_1=t55_l_n*t61_l_n;
t62_h=t62_h_1+t62_e;
t62_l_n=t62_l_n_1+t62_e;
}
else if(t61_l_n>=0){
REAL t62_l_n_1=t55_l_n*t61_h;
REAL t62_h_1=t55_l_n*t61_l_n;
t62_h=t62_h_1+t62_e;
t62_l_n=t62_l_n_1+t62_e;
}
else{
REAL t62_l_n_1=t55_l_n*t61_h;
REAL t55_h_n=-t55_h;
REAL t62_h_1=t55_h_n*t61_l_n;
t62_h=t62_h_1+t62_e;
t62_l_n=t62_l_n_1+t62_e;
}
}
else if(t55_l_n>=0){
if(t61_h<=0){
REAL t62_l_n_1=t55_h*t61_l_n;
REAL t62_h_1=t55_l_n*t61_l_n;
t62_h=t62_h_1+t62_e;
t62_l_n=t62_l_n_1+t62_e;
}
else if(t61_l_n>=0){
REAL t62_l_n_1_1=t55_l_n*t61_h;
REAL t62_l_n_1_2=t61_l_n*t55_h;
REAL t62_l_n_1=max(t62_l_n_1_1,t62_l_n_1_2);
REAL t62_h_1_1=t55_l_n*t61_l_n;
REAL t62_h_1_2=t55_h*t61_h;
REAL t62_h_1=max(t62_h_1_1,t62_h_1_2);
t62_h=t62_h_1+t62_e;
t62_l_n=t62_l_n_1+t62_e;
}
else{
REAL t62_l_n_1=t55_l_n*t61_h;
REAL t62_h_1=t55_h*t61_h;
t62_h=t62_h_1+t62_e;
t62_l_n=t62_l_n_1+t62_e;
}
}
else{
if(t61_h<=0){
REAL t62_l_n_1=t61_l_n*t55_h;
REAL t61_h_n=-t61_h;
REAL t62_h_1=t61_h_n*t55_l_n;
t62_h=t62_h_1+t62_e;
t62_l_n=t62_l_n_1+t62_e;
}
else if(t61_l_n>=0){
REAL t62_l_n_1=t61_l_n*t55_h;
REAL t62_h_1=t55_h*t61_h;
t62_h=t62_h_1+t62_e;
t62_l_n=t62_l_n_1+t62_e;
}
else{
REAL t55_l=-t55_l_n;
REAL t62_l_n_1=t55_l*t61_l_n;
REAL t62_h_1=t55_h*t61_h;
t62_h=t62_h_1+t62_e;
t62_l_n=t62_l_n_1+t62_e;
}
}
REAL t63_e_1=fabs(t63);
REAL t63_e=eps*t63_e_1;
REAL  t63_h, t63_l_n;
if(t10_h<=0){
if(t10_h<=0){
REAL t10_h_n=-t10_h;
REAL t63_l_n_1=t10_h_n*t10_h;
REAL t63_h_1=t10_l_n*t10_l_n;
t63_h=t63_h_1+t63_e;
t63_l_n=t63_l_n_1+t63_e;
}
else if(t10_l_n>=0){
REAL t63_l_n_1=t10_l_n*t10_h;
REAL t63_h_1=t10_l_n*t10_l_n;
t63_h=t63_h_1+t63_e;
t63_l_n=t63_l_n_1+t63_e;
}
else{
REAL t63_l_n_1=t10_l_n*t10_h;
REAL t10_h_n=-t10_h;
REAL t63_h_1=t10_h_n*t10_l_n;
t63_h=t63_h_1+t63_e;
t63_l_n=t63_l_n_1+t63_e;
}
}
else if(t10_l_n>=0){
if(t10_h<=0){
REAL t63_l_n_1=t10_h*t10_l_n;
REAL t63_h_1=t10_l_n*t10_l_n;
t63_h=t63_h_1+t63_e;
t63_l_n=t63_l_n_1+t63_e;
}
else if(t10_l_n>=0){
REAL t63_l_n_1_1=t10_l_n*t10_h;
REAL t63_l_n_1_2=t10_l_n*t10_h;
REAL t63_l_n_1=max(t63_l_n_1_1,t63_l_n_1_2);
REAL t63_h_1_1=t10_l_n*t10_l_n;
REAL t63_h_1_2=t10_h*t10_h;
REAL t63_h_1=max(t63_h_1_1,t63_h_1_2);
t63_h=t63_h_1+t63_e;
t63_l_n=t63_l_n_1+t63_e;
}
else{
REAL t63_l_n_1=t10_l_n*t10_h;
REAL t63_h_1=t10_h*t10_h;
t63_h=t63_h_1+t63_e;
t63_l_n=t63_l_n_1+t63_e;
}
}
else{
if(t10_h<=0){
REAL t63_l_n_1=t10_l_n*t10_h;
REAL t10_h_n=-t10_h;
REAL t63_h_1=t10_h_n*t10_l_n;
t63_h=t63_h_1+t63_e;
t63_l_n=t63_l_n_1+t63_e;
}
else if(t10_l_n>=0){
REAL t63_l_n_1=t10_l_n*t10_h;
REAL t63_h_1=t10_h*t10_h;
t63_h=t63_h_1+t63_e;
t63_l_n=t63_l_n_1+t63_e;
}
else{
REAL t10_l=-t10_l_n;
REAL t63_l_n_1=t10_l*t10_l_n;
REAL t63_h_1=t10_h*t10_h;
t63_h=t63_h_1+t63_e;
t63_l_n=t63_l_n_1+t63_e;
}
}
REAL t64_e_1=fabs(t64);
REAL t64_e=eps*t64_e_1;
REAL  t64_h, t64_l_n;
if(t15_h<=0){
if(t15_h<=0){
REAL t15_h_n=-t15_h;
REAL t64_l_n_1=t15_h_n*t15_h;
REAL t64_h_1=t15_l_n*t15_l_n;
t64_h=t64_h_1+t64_e;
t64_l_n=t64_l_n_1+t64_e;
}
else if(t15_l_n>=0){
REAL t64_l_n_1=t15_l_n*t15_h;
REAL t64_h_1=t15_l_n*t15_l_n;
t64_h=t64_h_1+t64_e;
t64_l_n=t64_l_n_1+t64_e;
}
else{
REAL t64_l_n_1=t15_l_n*t15_h;
REAL t15_h_n=-t15_h;
REAL t64_h_1=t15_h_n*t15_l_n;
t64_h=t64_h_1+t64_e;
t64_l_n=t64_l_n_1+t64_e;
}
}
else if(t15_l_n>=0){
if(t15_h<=0){
REAL t64_l_n_1=t15_h*t15_l_n;
REAL t64_h_1=t15_l_n*t15_l_n;
t64_h=t64_h_1+t64_e;
t64_l_n=t64_l_n_1+t64_e;
}
else if(t15_l_n>=0){
REAL t64_l_n_1_1=t15_l_n*t15_h;
REAL t64_l_n_1_2=t15_l_n*t15_h;
REAL t64_l_n_1=max(t64_l_n_1_1,t64_l_n_1_2);
REAL t64_h_1_1=t15_l_n*t15_l_n;
REAL t64_h_1_2=t15_h*t15_h;
REAL t64_h_1=max(t64_h_1_1,t64_h_1_2);
t64_h=t64_h_1+t64_e;
t64_l_n=t64_l_n_1+t64_e;
}
else{
REAL t64_l_n_1=t15_l_n*t15_h;
REAL t64_h_1=t15_h*t15_h;
t64_h=t64_h_1+t64_e;
t64_l_n=t64_l_n_1+t64_e;
}
}
else{
if(t15_h<=0){
REAL t64_l_n_1=t15_l_n*t15_h;
REAL t15_h_n=-t15_h;
REAL t64_h_1=t15_h_n*t15_l_n;
t64_h=t64_h_1+t64_e;
t64_l_n=t64_l_n_1+t64_e;
}
else if(t15_l_n>=0){
REAL t64_l_n_1=t15_l_n*t15_h;
REAL t64_h_1=t15_h*t15_h;
t64_h=t64_h_1+t64_e;
t64_l_n=t64_l_n_1+t64_e;
}
else{
REAL t15_l=-t15_l_n;
REAL t64_l_n_1=t15_l*t15_l_n;
REAL t64_h_1=t15_h*t15_h;
t64_h=t64_h_1+t64_e;
t64_l_n=t64_l_n_1+t64_e;
}
}
REAL t65_e_1=fabs(t65);
REAL t65_e=eps*t65_e_1;
REAL  t65_h, t65_l_n;
if(t8_h<=0){
if(t8_h<=0){
REAL t8_h_n=-t8_h;
REAL t65_l_n_1=t8_h_n*t8_h;
REAL t65_h_1=t8_l_n*t8_l_n;
t65_h=t65_h_1+t65_e;
t65_l_n=t65_l_n_1+t65_e;
}
else if(t8_l_n>=0){
REAL t65_l_n_1=t8_l_n*t8_h;
REAL t65_h_1=t8_l_n*t8_l_n;
t65_h=t65_h_1+t65_e;
t65_l_n=t65_l_n_1+t65_e;
}
else{
REAL t65_l_n_1=t8_l_n*t8_h;
REAL t8_h_n=-t8_h;
REAL t65_h_1=t8_h_n*t8_l_n;
t65_h=t65_h_1+t65_e;
t65_l_n=t65_l_n_1+t65_e;
}
}
else if(t8_l_n>=0){
if(t8_h<=0){
REAL t65_l_n_1=t8_h*t8_l_n;
REAL t65_h_1=t8_l_n*t8_l_n;
t65_h=t65_h_1+t65_e;
t65_l_n=t65_l_n_1+t65_e;
}
else if(t8_l_n>=0){
REAL t65_l_n_1_1=t8_l_n*t8_h;
REAL t65_l_n_1_2=t8_l_n*t8_h;
REAL t65_l_n_1=max(t65_l_n_1_1,t65_l_n_1_2);
REAL t65_h_1_1=t8_l_n*t8_l_n;
REAL t65_h_1_2=t8_h*t8_h;
REAL t65_h_1=max(t65_h_1_1,t65_h_1_2);
t65_h=t65_h_1+t65_e;
t65_l_n=t65_l_n_1+t65_e;
}
else{
REAL t65_l_n_1=t8_l_n*t8_h;
REAL t65_h_1=t8_h*t8_h;
t65_h=t65_h_1+t65_e;
t65_l_n=t65_l_n_1+t65_e;
}
}
else{
if(t8_h<=0){
REAL t65_l_n_1=t8_l_n*t8_h;
REAL t8_h_n=-t8_h;
REAL t65_h_1=t8_h_n*t8_l_n;
t65_h=t65_h_1+t65_e;
t65_l_n=t65_l_n_1+t65_e;
}
else if(t8_l_n>=0){
REAL t65_l_n_1=t8_l_n*t8_h;
REAL t65_h_1=t8_h*t8_h;
t65_h=t65_h_1+t65_e;
t65_l_n=t65_l_n_1+t65_e;
}
else{
REAL t8_l=-t8_l_n;
REAL t65_l_n_1=t8_l*t8_l_n;
REAL t65_h_1=t8_h*t8_h;
t65_h=t65_h_1+t65_e;
t65_l_n=t65_l_n_1+t65_e;
}
}
REAL t66_e_1=fabs(t66);
REAL t66_e=eps*t66_e_1;
REAL t66_h_1=t63_h+t64_h;
REAL t66_h=t66_h_1+t66_e;
REAL t66_l_n_1=t63_l_n+t64_l_n;
REAL t66_l_n=t66_l_n_1+t66_e;
REAL t67_e_1=fabs(t67);
REAL t67_e=eps*t67_e_1;
REAL t67_h_1=t65_h+t66_h;
REAL t67_h=t67_h_1+t67_e;
REAL t67_l_n_1=t65_l_n+t66_l_n;
REAL t67_l_n=t67_l_n_1+t67_e;
REAL t68_e_1=fabs(t68);
REAL t68_e=eps*t68_e_1;
REAL  t68_h, t68_l_n;
if(t21_h<=0){
if(t5_h<=0){
REAL t21_h_n=-t21_h;
REAL t68_l_n_1=t21_h_n*t5_h;
REAL t68_h_1=t21_l_n*t5_l_n;
t68_h=t68_h_1+t68_e;
t68_l_n=t68_l_n_1+t68_e;
}
else if(t5_l_n>=0){
REAL t68_l_n_1=t21_l_n*t5_h;
REAL t68_h_1=t21_l_n*t5_l_n;
t68_h=t68_h_1+t68_e;
t68_l_n=t68_l_n_1+t68_e;
}
else{
REAL t68_l_n_1=t21_l_n*t5_h;
REAL t21_h_n=-t21_h;
REAL t68_h_1=t21_h_n*t5_l_n;
t68_h=t68_h_1+t68_e;
t68_l_n=t68_l_n_1+t68_e;
}
}
else if(t21_l_n>=0){
if(t5_h<=0){
REAL t68_l_n_1=t21_h*t5_l_n;
REAL t68_h_1=t21_l_n*t5_l_n;
t68_h=t68_h_1+t68_e;
t68_l_n=t68_l_n_1+t68_e;
}
else if(t5_l_n>=0){
REAL t68_l_n_1_1=t21_l_n*t5_h;
REAL t68_l_n_1_2=t5_l_n*t21_h;
REAL t68_l_n_1=max(t68_l_n_1_1,t68_l_n_1_2);
REAL t68_h_1_1=t21_l_n*t5_l_n;
REAL t68_h_1_2=t21_h*t5_h;
REAL t68_h_1=max(t68_h_1_1,t68_h_1_2);
t68_h=t68_h_1+t68_e;
t68_l_n=t68_l_n_1+t68_e;
}
else{
REAL t68_l_n_1=t21_l_n*t5_h;
REAL t68_h_1=t21_h*t5_h;
t68_h=t68_h_1+t68_e;
t68_l_n=t68_l_n_1+t68_e;
}
}
else{
if(t5_h<=0){
REAL t68_l_n_1=t5_l_n*t21_h;
REAL t5_h_n=-t5_h;
REAL t68_h_1=t5_h_n*t21_l_n;
t68_h=t68_h_1+t68_e;
t68_l_n=t68_l_n_1+t68_e;
}
else if(t5_l_n>=0){
REAL t68_l_n_1=t5_l_n*t21_h;
REAL t68_h_1=t21_h*t5_h;
t68_h=t68_h_1+t68_e;
t68_l_n=t68_l_n_1+t68_e;
}
else{
REAL t21_l=-t21_l_n;
REAL t68_l_n_1=t21_l*t5_l_n;
REAL t68_h_1=t21_h*t5_h;
t68_h=t68_h_1+t68_e;
t68_l_n=t68_l_n_1+t68_e;
}
}
REAL t69_e_1=fabs(t69);
REAL t69_e=eps*t69_e_1;
REAL  t69_h, t69_l_n;
if(t38_h<=0){
if(t6_h<=0){
REAL t38_h_n=-t38_h;
REAL t69_l_n_1=t38_h_n*t6_h;
REAL t69_h_1=t38_l_n*t6_l_n;
t69_h=t69_h_1+t69_e;
t69_l_n=t69_l_n_1+t69_e;
}
else if(t6_l_n>=0){
REAL t69_l_n_1=t38_l_n*t6_h;
REAL t69_h_1=t38_l_n*t6_l_n;
t69_h=t69_h_1+t69_e;
t69_l_n=t69_l_n_1+t69_e;
}
else{
REAL t69_l_n_1=t38_l_n*t6_h;
REAL t38_h_n=-t38_h;
REAL t69_h_1=t38_h_n*t6_l_n;
t69_h=t69_h_1+t69_e;
t69_l_n=t69_l_n_1+t69_e;
}
}
else if(t38_l_n>=0){
if(t6_h<=0){
REAL t69_l_n_1=t38_h*t6_l_n;
REAL t69_h_1=t38_l_n*t6_l_n;
t69_h=t69_h_1+t69_e;
t69_l_n=t69_l_n_1+t69_e;
}
else if(t6_l_n>=0){
REAL t69_l_n_1_1=t38_l_n*t6_h;
REAL t69_l_n_1_2=t6_l_n*t38_h;
REAL t69_l_n_1=max(t69_l_n_1_1,t69_l_n_1_2);
REAL t69_h_1_1=t38_l_n*t6_l_n;
REAL t69_h_1_2=t38_h*t6_h;
REAL t69_h_1=max(t69_h_1_1,t69_h_1_2);
t69_h=t69_h_1+t69_e;
t69_l_n=t69_l_n_1+t69_e;
}
else{
REAL t69_l_n_1=t38_l_n*t6_h;
REAL t69_h_1=t38_h*t6_h;
t69_h=t69_h_1+t69_e;
t69_l_n=t69_l_n_1+t69_e;
}
}
else{
if(t6_h<=0){
REAL t69_l_n_1=t6_l_n*t38_h;
REAL t6_h_n=-t6_h;
REAL t69_h_1=t6_h_n*t38_l_n;
t69_h=t69_h_1+t69_e;
t69_l_n=t69_l_n_1+t69_e;
}
else if(t6_l_n>=0){
REAL t69_l_n_1=t6_l_n*t38_h;
REAL t69_h_1=t38_h*t6_h;
t69_h=t69_h_1+t69_e;
t69_l_n=t69_l_n_1+t69_e;
}
else{
REAL t38_l=-t38_l_n;
REAL t69_l_n_1=t38_l*t6_l_n;
REAL t69_h_1=t38_h*t6_h;
t69_h=t69_h_1+t69_e;
t69_l_n=t69_l_n_1+t69_e;
}
}
REAL t70_e_1=fabs(t70);
REAL t70_e=eps*t70_e_1;
REAL  t70_h, t70_l_n;
if(t22_h<=0){
if(t34_h<=0){
REAL t22_h_n=-t22_h;
REAL t70_l_n_1=t22_h_n*t34_h;
REAL t70_h_1=t22_l_n*t34_l_n;
t70_h=t70_h_1+t70_e;
t70_l_n=t70_l_n_1+t70_e;
}
else if(t34_l_n>=0){
REAL t70_l_n_1=t22_l_n*t34_h;
REAL t70_h_1=t22_l_n*t34_l_n;
t70_h=t70_h_1+t70_e;
t70_l_n=t70_l_n_1+t70_e;
}
else{
REAL t70_l_n_1=t22_l_n*t34_h;
REAL t22_h_n=-t22_h;
REAL t70_h_1=t22_h_n*t34_l_n;
t70_h=t70_h_1+t70_e;
t70_l_n=t70_l_n_1+t70_e;
}
}
else if(t22_l_n>=0){
if(t34_h<=0){
REAL t70_l_n_1=t22_h*t34_l_n;
REAL t70_h_1=t22_l_n*t34_l_n;
t70_h=t70_h_1+t70_e;
t70_l_n=t70_l_n_1+t70_e;
}
else if(t34_l_n>=0){
REAL t70_l_n_1_1=t22_l_n*t34_h;
REAL t70_l_n_1_2=t34_l_n*t22_h;
REAL t70_l_n_1=max(t70_l_n_1_1,t70_l_n_1_2);
REAL t70_h_1_1=t22_l_n*t34_l_n;
REAL t70_h_1_2=t22_h*t34_h;
REAL t70_h_1=max(t70_h_1_1,t70_h_1_2);
t70_h=t70_h_1+t70_e;
t70_l_n=t70_l_n_1+t70_e;
}
else{
REAL t70_l_n_1=t22_l_n*t34_h;
REAL t70_h_1=t22_h*t34_h;
t70_h=t70_h_1+t70_e;
t70_l_n=t70_l_n_1+t70_e;
}
}
else{
if(t34_h<=0){
REAL t70_l_n_1=t34_l_n*t22_h;
REAL t34_h_n=-t34_h;
REAL t70_h_1=t34_h_n*t22_l_n;
t70_h=t70_h_1+t70_e;
t70_l_n=t70_l_n_1+t70_e;
}
else if(t34_l_n>=0){
REAL t70_l_n_1=t34_l_n*t22_h;
REAL t70_h_1=t22_h*t34_h;
t70_h=t70_h_1+t70_e;
t70_l_n=t70_l_n_1+t70_e;
}
else{
REAL t22_l=-t22_l_n;
REAL t70_l_n_1=t22_l*t34_l_n;
REAL t70_h_1=t22_h*t34_h;
t70_h=t70_h_1+t70_e;
t70_l_n=t70_l_n_1+t70_e;
}
}
REAL t72_e_1=fabs(t72);
REAL t72_e=eps*t72_e_1;
REAL t72_h_1=t68_h+t69_h;
REAL t72_h=t72_h_1+t72_e;
REAL t72_l_n_1=t68_l_n+t69_l_n;
REAL t72_l_n=t72_l_n_1+t72_e;
REAL t73_e_1=fabs(t73);
REAL t73_e=eps*t73_e_1;
REAL t73_h_1=t72_h+t70_l_n;
REAL t73_h=t73_h_1+t73_e;
REAL t73_l_n_1=t70_h+t72_l_n;
REAL t73_l_n=t73_l_n_1+t73_e;
REAL t74_e_1=fabs(t74);
REAL t74_e=eps*t74_e_1;
REAL  t74_h, t74_l_n;
if(t67_h<=0){
if(t73_h<=0){
REAL t67_h_n=-t67_h;
REAL t74_l_n_1=t67_h_n*t73_h;
REAL t74_h_1=t67_l_n*t73_l_n;
t74_h=t74_h_1+t74_e;
t74_l_n=t74_l_n_1+t74_e;
}
else if(t73_l_n>=0){
REAL t74_l_n_1=t67_l_n*t73_h;
REAL t74_h_1=t67_l_n*t73_l_n;
t74_h=t74_h_1+t74_e;
t74_l_n=t74_l_n_1+t74_e;
}
else{
REAL t74_l_n_1=t67_l_n*t73_h;
REAL t67_h_n=-t67_h;
REAL t74_h_1=t67_h_n*t73_l_n;
t74_h=t74_h_1+t74_e;
t74_l_n=t74_l_n_1+t74_e;
}
}
else if(t67_l_n>=0){
if(t73_h<=0){
REAL t74_l_n_1=t67_h*t73_l_n;
REAL t74_h_1=t67_l_n*t73_l_n;
t74_h=t74_h_1+t74_e;
t74_l_n=t74_l_n_1+t74_e;
}
else if(t73_l_n>=0){
REAL t74_l_n_1_1=t67_l_n*t73_h;
REAL t74_l_n_1_2=t73_l_n*t67_h;
REAL t74_l_n_1=max(t74_l_n_1_1,t74_l_n_1_2);
REAL t74_h_1_1=t67_l_n*t73_l_n;
REAL t74_h_1_2=t67_h*t73_h;
REAL t74_h_1=max(t74_h_1_1,t74_h_1_2);
t74_h=t74_h_1+t74_e;
t74_l_n=t74_l_n_1+t74_e;
}
else{
REAL t74_l_n_1=t67_l_n*t73_h;
REAL t74_h_1=t67_h*t73_h;
t74_h=t74_h_1+t74_e;
t74_l_n=t74_l_n_1+t74_e;
}
}
else{
if(t73_h<=0){
REAL t74_l_n_1=t73_l_n*t67_h;
REAL t73_h_n=-t73_h;
REAL t74_h_1=t73_h_n*t67_l_n;
t74_h=t74_h_1+t74_e;
t74_l_n=t74_l_n_1+t74_e;
}
else if(t73_l_n>=0){
REAL t74_l_n_1=t73_l_n*t67_h;
REAL t74_h_1=t67_h*t73_h;
t74_h=t74_h_1+t74_e;
t74_l_n=t74_l_n_1+t74_e;
}
else{
REAL t67_l=-t67_l_n;
REAL t74_l_n_1=t67_l*t73_l_n;
REAL t74_h_1=t67_h*t73_h;
t74_h=t74_h_1+t74_e;
t74_l_n=t74_l_n_1+t74_e;
}
}
REAL t76_e_1=fabs(t76);
REAL t76_e=eps*t76_e_1;
REAL  t76_h, t76_l_n;
if(t16_h<=0){
if(t16_h<=0){
REAL t16_h_n=-t16_h;
REAL t76_l_n_1=t16_h_n*t16_h;
REAL t76_h_1=t16_l_n*t16_l_n;
t76_h=t76_h_1+t76_e;
t76_l_n=t76_l_n_1+t76_e;
}
else if(t16_l_n>=0){
REAL t76_l_n_1=t16_l_n*t16_h;
REAL t76_h_1=t16_l_n*t16_l_n;
t76_h=t76_h_1+t76_e;
t76_l_n=t76_l_n_1+t76_e;
}
else{
REAL t76_l_n_1=t16_l_n*t16_h;
REAL t16_h_n=-t16_h;
REAL t76_h_1=t16_h_n*t16_l_n;
t76_h=t76_h_1+t76_e;
t76_l_n=t76_l_n_1+t76_e;
}
}
else if(t16_l_n>=0){
if(t16_h<=0){
REAL t76_l_n_1=t16_h*t16_l_n;
REAL t76_h_1=t16_l_n*t16_l_n;
t76_h=t76_h_1+t76_e;
t76_l_n=t76_l_n_1+t76_e;
}
else if(t16_l_n>=0){
REAL t76_l_n_1_1=t16_l_n*t16_h;
REAL t76_l_n_1_2=t16_l_n*t16_h;
REAL t76_l_n_1=max(t76_l_n_1_1,t76_l_n_1_2);
REAL t76_h_1_1=t16_l_n*t16_l_n;
REAL t76_h_1_2=t16_h*t16_h;
REAL t76_h_1=max(t76_h_1_1,t76_h_1_2);
t76_h=t76_h_1+t76_e;
t76_l_n=t76_l_n_1+t76_e;
}
else{
REAL t76_l_n_1=t16_l_n*t16_h;
REAL t76_h_1=t16_h*t16_h;
t76_h=t76_h_1+t76_e;
t76_l_n=t76_l_n_1+t76_e;
}
}
else{
if(t16_h<=0){
REAL t76_l_n_1=t16_l_n*t16_h;
REAL t16_h_n=-t16_h;
REAL t76_h_1=t16_h_n*t16_l_n;
t76_h=t76_h_1+t76_e;
t76_l_n=t76_l_n_1+t76_e;
}
else if(t16_l_n>=0){
REAL t76_l_n_1=t16_l_n*t16_h;
REAL t76_h_1=t16_h*t16_h;
t76_h=t76_h_1+t76_e;
t76_l_n=t76_l_n_1+t76_e;
}
else{
REAL t16_l=-t16_l_n;
REAL t76_l_n_1=t16_l*t16_l_n;
REAL t76_h_1=t16_h*t16_h;
t76_h=t76_h_1+t76_e;
t76_l_n=t76_l_n_1+t76_e;
}
}
REAL t77_e_1=fabs(t77);
REAL t77_e=eps*t77_e_1;
REAL  t77_h, t77_l_n;
if(t17_h<=0){
if(t17_h<=0){
REAL t17_h_n=-t17_h;
REAL t77_l_n_1=t17_h_n*t17_h;
REAL t77_h_1=t17_l_n*t17_l_n;
t77_h=t77_h_1+t77_e;
t77_l_n=t77_l_n_1+t77_e;
}
else if(t17_l_n>=0){
REAL t77_l_n_1=t17_l_n*t17_h;
REAL t77_h_1=t17_l_n*t17_l_n;
t77_h=t77_h_1+t77_e;
t77_l_n=t77_l_n_1+t77_e;
}
else{
REAL t77_l_n_1=t17_l_n*t17_h;
REAL t17_h_n=-t17_h;
REAL t77_h_1=t17_h_n*t17_l_n;
t77_h=t77_h_1+t77_e;
t77_l_n=t77_l_n_1+t77_e;
}
}
else if(t17_l_n>=0){
if(t17_h<=0){
REAL t77_l_n_1=t17_h*t17_l_n;
REAL t77_h_1=t17_l_n*t17_l_n;
t77_h=t77_h_1+t77_e;
t77_l_n=t77_l_n_1+t77_e;
}
else if(t17_l_n>=0){
REAL t77_l_n_1_1=t17_l_n*t17_h;
REAL t77_l_n_1_2=t17_l_n*t17_h;
REAL t77_l_n_1=max(t77_l_n_1_1,t77_l_n_1_2);
REAL t77_h_1_1=t17_l_n*t17_l_n;
REAL t77_h_1_2=t17_h*t17_h;
REAL t77_h_1=max(t77_h_1_1,t77_h_1_2);
t77_h=t77_h_1+t77_e;
t77_l_n=t77_l_n_1+t77_e;
}
else{
REAL t77_l_n_1=t17_l_n*t17_h;
REAL t77_h_1=t17_h*t17_h;
t77_h=t77_h_1+t77_e;
t77_l_n=t77_l_n_1+t77_e;
}
}
else{
if(t17_h<=0){
REAL t77_l_n_1=t17_l_n*t17_h;
REAL t17_h_n=-t17_h;
REAL t77_h_1=t17_h_n*t17_l_n;
t77_h=t77_h_1+t77_e;
t77_l_n=t77_l_n_1+t77_e;
}
else if(t17_l_n>=0){
REAL t77_l_n_1=t17_l_n*t17_h;
REAL t77_h_1=t17_h*t17_h;
t77_h=t77_h_1+t77_e;
t77_l_n=t77_l_n_1+t77_e;
}
else{
REAL t17_l=-t17_l_n;
REAL t77_l_n_1=t17_l*t17_l_n;
REAL t77_h_1=t17_h*t17_h;
t77_h=t77_h_1+t77_e;
t77_l_n=t77_l_n_1+t77_e;
}
}
REAL t78_e_1=fabs(t78);
REAL t78_e=eps*t78_e_1;
REAL  t78_h, t78_l_n;
if(t6_h<=0){
if(t6_h<=0){
REAL t6_h_n=-t6_h;
REAL t78_l_n_1=t6_h_n*t6_h;
REAL t78_h_1=t6_l_n*t6_l_n;
t78_h=t78_h_1+t78_e;
t78_l_n=t78_l_n_1+t78_e;
}
else if(t6_l_n>=0){
REAL t78_l_n_1=t6_l_n*t6_h;
REAL t78_h_1=t6_l_n*t6_l_n;
t78_h=t78_h_1+t78_e;
t78_l_n=t78_l_n_1+t78_e;
}
else{
REAL t78_l_n_1=t6_l_n*t6_h;
REAL t6_h_n=-t6_h;
REAL t78_h_1=t6_h_n*t6_l_n;
t78_h=t78_h_1+t78_e;
t78_l_n=t78_l_n_1+t78_e;
}
}
else if(t6_l_n>=0){
if(t6_h<=0){
REAL t78_l_n_1=t6_h*t6_l_n;
REAL t78_h_1=t6_l_n*t6_l_n;
t78_h=t78_h_1+t78_e;
t78_l_n=t78_l_n_1+t78_e;
}
else if(t6_l_n>=0){
REAL t78_l_n_1_1=t6_l_n*t6_h;
REAL t78_l_n_1_2=t6_l_n*t6_h;
REAL t78_l_n_1=max(t78_l_n_1_1,t78_l_n_1_2);
REAL t78_h_1_1=t6_l_n*t6_l_n;
REAL t78_h_1_2=t6_h*t6_h;
REAL t78_h_1=max(t78_h_1_1,t78_h_1_2);
t78_h=t78_h_1+t78_e;
t78_l_n=t78_l_n_1+t78_e;
}
else{
REAL t78_l_n_1=t6_l_n*t6_h;
REAL t78_h_1=t6_h*t6_h;
t78_h=t78_h_1+t78_e;
t78_l_n=t78_l_n_1+t78_e;
}
}
else{
if(t6_h<=0){
REAL t78_l_n_1=t6_l_n*t6_h;
REAL t6_h_n=-t6_h;
REAL t78_h_1=t6_h_n*t6_l_n;
t78_h=t78_h_1+t78_e;
t78_l_n=t78_l_n_1+t78_e;
}
else if(t6_l_n>=0){
REAL t78_l_n_1=t6_l_n*t6_h;
REAL t78_h_1=t6_h*t6_h;
t78_h=t78_h_1+t78_e;
t78_l_n=t78_l_n_1+t78_e;
}
else{
REAL t6_l=-t6_l_n;
REAL t78_l_n_1=t6_l*t6_l_n;
REAL t78_h_1=t6_h*t6_h;
t78_h=t78_h_1+t78_e;
t78_l_n=t78_l_n_1+t78_e;
}
}
REAL t79_e_1=fabs(t79);
REAL t79_e=eps*t79_e_1;
REAL t79_h_1=t76_h+t77_h;
REAL t79_h=t79_h_1+t79_e;
REAL t79_l_n_1=t76_l_n+t77_l_n;
REAL t79_l_n=t79_l_n_1+t79_e;
REAL t80_e_1=fabs(t80);
REAL t80_e=eps*t80_e_1;
REAL t80_h_1=t78_h+t79_h;
REAL t80_h=t80_h_1+t80_e;
REAL t80_l_n_1=t78_l_n+t79_l_n;
REAL t80_l_n=t80_l_n_1+t80_e;
REAL t81_e_1=fabs(t81);
REAL t81_e=eps*t81_e_1;
REAL  t81_h, t81_l_n;
if(t14_h<=0){
if(t5_h<=0){
REAL t14_h_n=-t14_h;
REAL t81_l_n_1=t14_h_n*t5_h;
REAL t81_h_1=t14_l_n*t5_l_n;
t81_h=t81_h_1+t81_e;
t81_l_n=t81_l_n_1+t81_e;
}
else if(t5_l_n>=0){
REAL t81_l_n_1=t14_l_n*t5_h;
REAL t81_h_1=t14_l_n*t5_l_n;
t81_h=t81_h_1+t81_e;
t81_l_n=t81_l_n_1+t81_e;
}
else{
REAL t81_l_n_1=t14_l_n*t5_h;
REAL t14_h_n=-t14_h;
REAL t81_h_1=t14_h_n*t5_l_n;
t81_h=t81_h_1+t81_e;
t81_l_n=t81_l_n_1+t81_e;
}
}
else if(t14_l_n>=0){
if(t5_h<=0){
REAL t81_l_n_1=t14_h*t5_l_n;
REAL t81_h_1=t14_l_n*t5_l_n;
t81_h=t81_h_1+t81_e;
t81_l_n=t81_l_n_1+t81_e;
}
else if(t5_l_n>=0){
REAL t81_l_n_1_1=t14_l_n*t5_h;
REAL t81_l_n_1_2=t5_l_n*t14_h;
REAL t81_l_n_1=max(t81_l_n_1_1,t81_l_n_1_2);
REAL t81_h_1_1=t14_l_n*t5_l_n;
REAL t81_h_1_2=t14_h*t5_h;
REAL t81_h_1=max(t81_h_1_1,t81_h_1_2);
t81_h=t81_h_1+t81_e;
t81_l_n=t81_l_n_1+t81_e;
}
else{
REAL t81_l_n_1=t14_l_n*t5_h;
REAL t81_h_1=t14_h*t5_h;
t81_h=t81_h_1+t81_e;
t81_l_n=t81_l_n_1+t81_e;
}
}
else{
if(t5_h<=0){
REAL t81_l_n_1=t5_l_n*t14_h;
REAL t5_h_n=-t5_h;
REAL t81_h_1=t5_h_n*t14_l_n;
t81_h=t81_h_1+t81_e;
t81_l_n=t81_l_n_1+t81_e;
}
else if(t5_l_n>=0){
REAL t81_l_n_1=t5_l_n*t14_h;
REAL t81_h_1=t14_h*t5_h;
t81_h=t81_h_1+t81_e;
t81_l_n=t81_l_n_1+t81_e;
}
else{
REAL t14_l=-t14_l_n;
REAL t81_l_n_1=t14_l*t5_l_n;
REAL t81_h_1=t14_h*t5_h;
t81_h=t81_h_1+t81_e;
t81_l_n=t81_l_n_1+t81_e;
}
}
REAL t82_e_1=fabs(t82);
REAL t82_e=eps*t82_e_1;
REAL  t82_h, t82_l_n;
if(t22_h<=0){
if(t30_h<=0){
REAL t22_h_n=-t22_h;
REAL t82_l_n_1=t22_h_n*t30_h;
REAL t82_h_1=t22_l_n*t30_l_n;
t82_h=t82_h_1+t82_e;
t82_l_n=t82_l_n_1+t82_e;
}
else if(t30_l_n>=0){
REAL t82_l_n_1=t22_l_n*t30_h;
REAL t82_h_1=t22_l_n*t30_l_n;
t82_h=t82_h_1+t82_e;
t82_l_n=t82_l_n_1+t82_e;
}
else{
REAL t82_l_n_1=t22_l_n*t30_h;
REAL t22_h_n=-t22_h;
REAL t82_h_1=t22_h_n*t30_l_n;
t82_h=t82_h_1+t82_e;
t82_l_n=t82_l_n_1+t82_e;
}
}
else if(t22_l_n>=0){
if(t30_h<=0){
REAL t82_l_n_1=t22_h*t30_l_n;
REAL t82_h_1=t22_l_n*t30_l_n;
t82_h=t82_h_1+t82_e;
t82_l_n=t82_l_n_1+t82_e;
}
else if(t30_l_n>=0){
REAL t82_l_n_1_1=t22_l_n*t30_h;
REAL t82_l_n_1_2=t30_l_n*t22_h;
REAL t82_l_n_1=max(t82_l_n_1_1,t82_l_n_1_2);
REAL t82_h_1_1=t22_l_n*t30_l_n;
REAL t82_h_1_2=t22_h*t30_h;
REAL t82_h_1=max(t82_h_1_1,t82_h_1_2);
t82_h=t82_h_1+t82_e;
t82_l_n=t82_l_n_1+t82_e;
}
else{
REAL t82_l_n_1=t22_l_n*t30_h;
REAL t82_h_1=t22_h*t30_h;
t82_h=t82_h_1+t82_e;
t82_l_n=t82_l_n_1+t82_e;
}
}
else{
if(t30_h<=0){
REAL t82_l_n_1=t30_l_n*t22_h;
REAL t30_h_n=-t30_h;
REAL t82_h_1=t30_h_n*t22_l_n;
t82_h=t82_h_1+t82_e;
t82_l_n=t82_l_n_1+t82_e;
}
else if(t30_l_n>=0){
REAL t82_l_n_1=t30_l_n*t22_h;
REAL t82_h_1=t22_h*t30_h;
t82_h=t82_h_1+t82_e;
t82_l_n=t82_l_n_1+t82_e;
}
else{
REAL t22_l=-t22_l_n;
REAL t82_l_n_1=t22_l*t30_l_n;
REAL t82_h_1=t22_h*t30_h;
t82_h=t82_h_1+t82_e;
t82_l_n=t82_l_n_1+t82_e;
}
}
REAL t83_e_1=fabs(t83);
REAL t83_e=eps*t83_e_1;
REAL  t83_h, t83_l_n;
if(t15_h<=0){
if(t38_h<=0){
REAL t15_h_n=-t15_h;
REAL t83_l_n_1=t15_h_n*t38_h;
REAL t83_h_1=t15_l_n*t38_l_n;
t83_h=t83_h_1+t83_e;
t83_l_n=t83_l_n_1+t83_e;
}
else if(t38_l_n>=0){
REAL t83_l_n_1=t15_l_n*t38_h;
REAL t83_h_1=t15_l_n*t38_l_n;
t83_h=t83_h_1+t83_e;
t83_l_n=t83_l_n_1+t83_e;
}
else{
REAL t83_l_n_1=t15_l_n*t38_h;
REAL t15_h_n=-t15_h;
REAL t83_h_1=t15_h_n*t38_l_n;
t83_h=t83_h_1+t83_e;
t83_l_n=t83_l_n_1+t83_e;
}
}
else if(t15_l_n>=0){
if(t38_h<=0){
REAL t83_l_n_1=t15_h*t38_l_n;
REAL t83_h_1=t15_l_n*t38_l_n;
t83_h=t83_h_1+t83_e;
t83_l_n=t83_l_n_1+t83_e;
}
else if(t38_l_n>=0){
REAL t83_l_n_1_1=t15_l_n*t38_h;
REAL t83_l_n_1_2=t38_l_n*t15_h;
REAL t83_l_n_1=max(t83_l_n_1_1,t83_l_n_1_2);
REAL t83_h_1_1=t15_l_n*t38_l_n;
REAL t83_h_1_2=t15_h*t38_h;
REAL t83_h_1=max(t83_h_1_1,t83_h_1_2);
t83_h=t83_h_1+t83_e;
t83_l_n=t83_l_n_1+t83_e;
}
else{
REAL t83_l_n_1=t15_l_n*t38_h;
REAL t83_h_1=t15_h*t38_h;
t83_h=t83_h_1+t83_e;
t83_l_n=t83_l_n_1+t83_e;
}
}
else{
if(t38_h<=0){
REAL t83_l_n_1=t38_l_n*t15_h;
REAL t38_h_n=-t38_h;
REAL t83_h_1=t38_h_n*t15_l_n;
t83_h=t83_h_1+t83_e;
t83_l_n=t83_l_n_1+t83_e;
}
else if(t38_l_n>=0){
REAL t83_l_n_1=t38_l_n*t15_h;
REAL t83_h_1=t15_h*t38_h;
t83_h=t83_h_1+t83_e;
t83_l_n=t83_l_n_1+t83_e;
}
else{
REAL t15_l=-t15_l_n;
REAL t83_l_n_1=t15_l*t38_l_n;
REAL t83_h_1=t15_h*t38_h;
t83_h=t83_h_1+t83_e;
t83_l_n=t83_l_n_1+t83_e;
}
}
REAL t85_e_1=fabs(t85);
REAL t85_e=eps*t85_e_1;
REAL t85_h_1=t81_h+t82_h;
REAL t85_h=t85_h_1+t85_e;
REAL t85_l_n_1=t81_l_n+t82_l_n;
REAL t85_l_n=t85_l_n_1+t85_e;
REAL t86_e_1=fabs(t86);
REAL t86_e=eps*t86_e_1;
REAL t86_h_1=t85_h+t83_l_n;
REAL t86_h=t86_h_1+t86_e;
REAL t86_l_n_1=t83_h+t85_l_n;
REAL t86_l_n=t86_l_n_1+t86_e;
REAL t87_e_1=fabs(t87);
REAL t87_e=eps*t87_e_1;
REAL  t87_h, t87_l_n;
if(t80_h<=0){
if(t86_h<=0){
REAL t80_h_n=-t80_h;
REAL t87_l_n_1=t80_h_n*t86_h;
REAL t87_h_1=t80_l_n*t86_l_n;
t87_h=t87_h_1+t87_e;
t87_l_n=t87_l_n_1+t87_e;
}
else if(t86_l_n>=0){
REAL t87_l_n_1=t80_l_n*t86_h;
REAL t87_h_1=t80_l_n*t86_l_n;
t87_h=t87_h_1+t87_e;
t87_l_n=t87_l_n_1+t87_e;
}
else{
REAL t87_l_n_1=t80_l_n*t86_h;
REAL t80_h_n=-t80_h;
REAL t87_h_1=t80_h_n*t86_l_n;
t87_h=t87_h_1+t87_e;
t87_l_n=t87_l_n_1+t87_e;
}
}
else if(t80_l_n>=0){
if(t86_h<=0){
REAL t87_l_n_1=t80_h*t86_l_n;
REAL t87_h_1=t80_l_n*t86_l_n;
t87_h=t87_h_1+t87_e;
t87_l_n=t87_l_n_1+t87_e;
}
else if(t86_l_n>=0){
REAL t87_l_n_1_1=t80_l_n*t86_h;
REAL t87_l_n_1_2=t86_l_n*t80_h;
REAL t87_l_n_1=max(t87_l_n_1_1,t87_l_n_1_2);
REAL t87_h_1_1=t80_l_n*t86_l_n;
REAL t87_h_1_2=t80_h*t86_h;
REAL t87_h_1=max(t87_h_1_1,t87_h_1_2);
t87_h=t87_h_1+t87_e;
t87_l_n=t87_l_n_1+t87_e;
}
else{
REAL t87_l_n_1=t80_l_n*t86_h;
REAL t87_h_1=t80_h*t86_h;
t87_h=t87_h_1+t87_e;
t87_l_n=t87_l_n_1+t87_e;
}
}
else{
if(t86_h<=0){
REAL t87_l_n_1=t86_l_n*t80_h;
REAL t86_h_n=-t86_h;
REAL t87_h_1=t86_h_n*t80_l_n;
t87_h=t87_h_1+t87_e;
t87_l_n=t87_l_n_1+t87_e;
}
else if(t86_l_n>=0){
REAL t87_l_n_1=t86_l_n*t80_h;
REAL t87_h_1=t80_h*t86_h;
t87_h=t87_h_1+t87_e;
t87_l_n=t87_l_n_1+t87_e;
}
else{
REAL t80_l=-t80_l_n;
REAL t87_l_n_1=t80_l*t86_l_n;
REAL t87_h_1=t80_h*t86_h;
t87_h=t87_h_1+t87_e;
t87_l_n=t87_l_n_1+t87_e;
}
}
REAL t89_e_1=fabs(t89);
REAL t89_e=eps*t89_e_1;
REAL t89_h_1=t50_h+t62_h;
REAL t89_h=t89_h_1+t89_e;
REAL t89_l_n_1=t50_l_n+t62_l_n;
REAL t89_l_n=t89_l_n_1+t89_e;
REAL t90_e_1=fabs(t90);
REAL t90_e=eps*t90_e_1;
REAL t90_h_1=t89_h+t74_l_n;
REAL t90_h=t90_h_1+t90_e;
REAL t90_l_n_1=t74_h+t89_l_n;
REAL t90_l_n=t90_l_n_1+t90_e;
REAL t91_e_1=fabs(t91);
REAL t91_e=eps*t91_e_1;
REAL t91_h_1=t90_h+t87_l_n;
REAL t91_h=t91_h_1+t91_e;
REAL t91_l_n_1=t87_h+t90_l_n;
REAL t91_l_n=t91_l_n_1+t91_e;
fesetround(roundmode);
*result=t91;
if (inferr==0) {
*hbound=t91_h;
*lbound=-t91_l_n;
}
else{
*hbound=1/0;
*lbound=1/0;
}
}
//Main test Function
int main(){
int sample=10000;
int numtest=100;
FILE *save;
save=fopen("SingleSave.txt","w");
FILE *report;
report=fopen("SingleReport.txt","w");
REAL range=1000000;
REAL varlist;
REAL res, err,herr,lerr;
REAL FPStotaltime=0; 
REAL IAtotaltime=0;
REAL Manualtotaltime=0;
REAL FPStime=0; 
REAL IAtime=0;
REAL Manualtime=0;
clock_t before, after;
REAL ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez;
int i,j;
for(j=0; j<sample; j++){
FPStime=0;
IAtime=0;
Manualtime=0;

ax=randfl_uniReal(range);
ay=randfl_uniReal(range);
az=randfl_uniReal(range);
bx=randfl_uniReal(range);
by=randfl_uniReal(range);
bz=randfl_uniReal(range);
cx=randfl_uniReal(range);
cy=randfl_uniReal(range);
cz=randfl_uniReal(range);
dx=randfl_uniReal(range);
dy=randfl_uniReal(range);
dz=randfl_uniReal(range);
ex=randfl_uniReal(range);
ey=randfl_uniReal(range);
ez=randfl_uniReal(range);

before = clock();
FPSyn_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
FPSyn_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
FPSyn_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
FPSyn_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
FPSyn_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
FPSyn_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
FPSyn_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
FPSyn_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
FPSyn_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
FPSyn_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
FPSyn_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
FPSyn_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
FPSyn_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
FPSyn_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
FPSyn_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
FPSyn_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
FPSyn_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
FPSyn_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
FPSyn_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
FPSyn_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
FPSyn_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
FPSyn_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
FPSyn_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
FPSyn_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
FPSyn_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
FPSyn_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
FPSyn_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
FPSyn_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
FPSyn_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
FPSyn_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
FPSyn_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
FPSyn_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
FPSyn_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
FPSyn_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
FPSyn_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
FPSyn_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
FPSyn_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
FPSyn_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
FPSyn_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
FPSyn_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
FPSyn_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
FPSyn_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
FPSyn_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
FPSyn_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
FPSyn_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
FPSyn_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
FPSyn_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
FPSyn_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
FPSyn_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
FPSyn_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
FPSyn_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
FPSyn_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
FPSyn_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
FPSyn_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
FPSyn_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
FPSyn_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
FPSyn_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
FPSyn_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
FPSyn_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
FPSyn_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
FPSyn_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
FPSyn_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
FPSyn_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
FPSyn_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
FPSyn_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
FPSyn_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
FPSyn_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
FPSyn_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
FPSyn_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
FPSyn_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
FPSyn_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
FPSyn_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
FPSyn_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
FPSyn_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
FPSyn_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
FPSyn_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
FPSyn_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
FPSyn_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
FPSyn_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
FPSyn_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
FPSyn_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
FPSyn_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
FPSyn_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
FPSyn_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
FPSyn_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
FPSyn_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
FPSyn_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
FPSyn_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
FPSyn_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
FPSyn_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
FPSyn_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
FPSyn_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
FPSyn_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
FPSyn_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
FPSyn_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
FPSyn_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
FPSyn_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
FPSyn_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
FPSyn_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
FPSyn_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
FPSyn_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
FPSyn_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
FPSyn_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
FPSyn_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
FPSyn_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
FPSyn_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
FPSyn_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
FPSyn_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
FPSyn_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
FPSyn_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
FPSyn_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
FPSyn_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
FPSyn_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
FPSyn_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
FPSyn_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
FPSyn_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
FPSyn_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
FPSyn_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
FPSyn_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
FPSyn_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
FPSyn_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
FPSyn_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
FPSyn_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
FPSyn_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
FPSyn_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
FPSyn_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
FPSyn_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
FPSyn_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
FPSyn_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
FPSyn_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
FPSyn_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
FPSyn_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
FPSyn_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
FPSyn_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
FPSyn_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
FPSyn_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
FPSyn_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
FPSyn_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
FPSyn_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
FPSyn_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
FPSyn_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
FPSyn_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
FPSyn_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
FPSyn_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
FPSyn_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
FPSyn_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
FPSyn_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
FPSyn_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
FPSyn_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
FPSyn_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
FPSyn_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
FPSyn_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
FPSyn_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
FPSyn_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
FPSyn_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
FPSyn_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
FPSyn_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
FPSyn_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
FPSyn_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
FPSyn_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
FPSyn_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
FPSyn_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
FPSyn_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
FPSyn_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
FPSyn_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
FPSyn_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
FPSyn_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
FPSyn_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
FPSyn_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
FPSyn_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
FPSyn_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
FPSyn_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
FPSyn_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
FPSyn_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
FPSyn_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
FPSyn_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
FPSyn_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
FPSyn_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
FPSyn_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
FPSyn_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
FPSyn_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
FPSyn_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
FPSyn_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
FPSyn_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
FPSyn_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
FPSyn_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
FPSyn_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
FPSyn_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
FPSyn_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
FPSyn_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
FPSyn_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
FPSyn_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
FPSyn_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
FPSyn_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
FPSyn_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
FPSyn_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
FPSyn_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
FPSyn_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
FPSyn_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
FPSyn_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
FPSyn_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
FPSyn_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
FPSyn_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
FPSyn_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
FPSyn_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
FPSyn_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
FPSyn_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
FPSyn_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
FPSyn_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
FPSyn_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
FPSyn_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
FPSyn_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
FPSyn_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
FPSyn_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
FPSyn_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
FPSyn_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
FPSyn_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
FPSyn_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
FPSyn_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
FPSyn_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
FPSyn_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
FPSyn_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
FPSyn_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
FPSyn_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
FPSyn_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
FPSyn_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
FPSyn_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
FPSyn_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
FPSyn_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
FPSyn_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
FPSyn_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
FPSyn_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
FPSyn_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
FPSyn_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
FPSyn_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
FPSyn_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
FPSyn_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
FPSyn_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
FPSyn_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
FPSyn_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
FPSyn_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
FPSyn_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
FPSyn_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
FPSyn_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
FPSyn_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
FPSyn_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
FPSyn_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
FPSyn_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
FPSyn_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
FPSyn_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
FPSyn_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
FPSyn_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
FPSyn_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
FPSyn_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
FPSyn_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
FPSyn_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
FPSyn_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
FPSyn_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
FPSyn_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
FPSyn_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
FPSyn_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
FPSyn_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
FPSyn_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
FPSyn_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
FPSyn_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
FPSyn_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
FPSyn_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
FPSyn_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
FPSyn_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
FPSyn_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
FPSyn_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
FPSyn_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
FPSyn_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
FPSyn_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
FPSyn_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
FPSyn_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
FPSyn_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
FPSyn_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
FPSyn_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
FPSyn_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
FPSyn_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
FPSyn_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
FPSyn_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
FPSyn_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
FPSyn_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
FPSyn_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
FPSyn_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
FPSyn_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
FPSyn_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
FPSyn_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
FPSyn_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
FPSyn_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
FPSyn_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
FPSyn_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
FPSyn_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
FPSyn_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
FPSyn_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
FPSyn_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
FPSyn_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
FPSyn_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
FPSyn_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
FPSyn_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
FPSyn_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
FPSyn_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
FPSyn_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
FPSyn_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
FPSyn_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
FPSyn_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
FPSyn_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
FPSyn_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
FPSyn_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
FPSyn_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
FPSyn_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
FPSyn_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
FPSyn_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
FPSyn_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
FPSyn_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
FPSyn_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
FPSyn_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
FPSyn_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
FPSyn_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
FPSyn_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
FPSyn_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
FPSyn_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
FPSyn_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
FPSyn_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
FPSyn_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
FPSyn_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
FPSyn_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
FPSyn_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
FPSyn_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
FPSyn_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
FPSyn_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
FPSyn_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
FPSyn_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
FPSyn_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
FPSyn_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
FPSyn_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
FPSyn_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
FPSyn_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
FPSyn_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
FPSyn_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
FPSyn_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
FPSyn_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
FPSyn_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
FPSyn_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
FPSyn_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
FPSyn_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
FPSyn_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
FPSyn_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
FPSyn_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
FPSyn_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
FPSyn_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
FPSyn_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
FPSyn_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
FPSyn_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
FPSyn_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
FPSyn_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
FPSyn_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
FPSyn_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
FPSyn_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
FPSyn_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
FPSyn_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
FPSyn_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
FPSyn_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
FPSyn_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
FPSyn_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
FPSyn_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
FPSyn_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
FPSyn_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
FPSyn_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
FPSyn_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
FPSyn_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
FPSyn_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
FPSyn_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
FPSyn_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
FPSyn_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
FPSyn_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
FPSyn_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
FPSyn_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
FPSyn_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
FPSyn_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
FPSyn_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
FPSyn_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
FPSyn_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
FPSyn_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
FPSyn_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
FPSyn_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
FPSyn_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
FPSyn_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
FPSyn_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
FPSyn_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
FPSyn_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
FPSyn_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
FPSyn_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
FPSyn_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
FPSyn_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
FPSyn_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
FPSyn_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
FPSyn_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
FPSyn_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
FPSyn_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
FPSyn_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
FPSyn_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
FPSyn_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
FPSyn_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
FPSyn_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
FPSyn_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
FPSyn_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
FPSyn_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
FPSyn_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
FPSyn_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
FPSyn_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
FPSyn_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
FPSyn_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
FPSyn_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
FPSyn_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
FPSyn_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
FPSyn_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
FPSyn_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
FPSyn_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
FPSyn_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
FPSyn_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
FPSyn_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
FPSyn_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
FPSyn_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
FPSyn_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
FPSyn_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
FPSyn_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
FPSyn_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
FPSyn_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
FPSyn_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
FPSyn_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
FPSyn_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
FPSyn_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
FPSyn_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
FPSyn_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
FPSyn_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
FPSyn_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
FPSyn_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
FPSyn_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
FPSyn_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
FPSyn_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
FPSyn_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
FPSyn_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
FPSyn_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
FPSyn_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
FPSyn_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
FPSyn_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
FPSyn_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
FPSyn_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
FPSyn_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
FPSyn_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
FPSyn_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
FPSyn_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
FPSyn_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
FPSyn_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
FPSyn_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
FPSyn_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
FPSyn_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
FPSyn_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
FPSyn_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
FPSyn_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
FPSyn_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
FPSyn_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
FPSyn_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
FPSyn_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
FPSyn_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
FPSyn_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
FPSyn_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
FPSyn_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
FPSyn_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
FPSyn_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
FPSyn_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
FPSyn_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
FPSyn_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
FPSyn_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
FPSyn_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
FPSyn_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
FPSyn_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
FPSyn_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
FPSyn_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
FPSyn_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
FPSyn_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
FPSyn_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
FPSyn_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
FPSyn_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
FPSyn_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
FPSyn_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
FPSyn_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
FPSyn_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
FPSyn_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
FPSyn_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
FPSyn_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
FPSyn_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
FPSyn_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
FPSyn_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
FPSyn_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
FPSyn_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
FPSyn_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
FPSyn_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
FPSyn_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
FPSyn_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
FPSyn_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
FPSyn_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
FPSyn_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
FPSyn_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
FPSyn_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
FPSyn_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
FPSyn_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
FPSyn_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
FPSyn_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
FPSyn_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
FPSyn_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
FPSyn_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
FPSyn_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
FPSyn_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
FPSyn_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
FPSyn_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
FPSyn_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
FPSyn_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
FPSyn_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
FPSyn_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
FPSyn_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
FPSyn_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
FPSyn_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
FPSyn_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
FPSyn_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
FPSyn_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
FPSyn_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
FPSyn_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
FPSyn_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
FPSyn_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
FPSyn_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
FPSyn_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
FPSyn_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
FPSyn_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
FPSyn_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
FPSyn_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
FPSyn_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
FPSyn_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
FPSyn_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
FPSyn_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
FPSyn_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
FPSyn_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
FPSyn_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
FPSyn_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
FPSyn_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
FPSyn_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
FPSyn_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
FPSyn_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
FPSyn_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
FPSyn_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
FPSyn_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
FPSyn_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
FPSyn_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
FPSyn_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
FPSyn_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
FPSyn_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
FPSyn_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
FPSyn_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
FPSyn_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
FPSyn_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
FPSyn_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
FPSyn_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
FPSyn_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
FPSyn_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
FPSyn_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
FPSyn_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
FPSyn_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
FPSyn_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
FPSyn_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
FPSyn_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
FPSyn_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
FPSyn_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
FPSyn_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
FPSyn_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
FPSyn_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
FPSyn_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
FPSyn_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
FPSyn_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
FPSyn_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
FPSyn_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
FPSyn_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
FPSyn_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
FPSyn_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
FPSyn_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
FPSyn_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
FPSyn_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
FPSyn_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
FPSyn_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
FPSyn_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
FPSyn_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
FPSyn_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
FPSyn_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
FPSyn_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
FPSyn_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
FPSyn_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
FPSyn_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
FPSyn_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
FPSyn_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
FPSyn_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
FPSyn_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
FPSyn_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
FPSyn_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
FPSyn_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
FPSyn_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
FPSyn_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
FPSyn_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
FPSyn_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
FPSyn_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
FPSyn_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
FPSyn_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
FPSyn_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
FPSyn_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
FPSyn_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
FPSyn_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
FPSyn_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
FPSyn_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
FPSyn_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
FPSyn_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
FPSyn_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
FPSyn_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
FPSyn_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
FPSyn_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
FPSyn_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
FPSyn_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
FPSyn_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
FPSyn_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
FPSyn_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
FPSyn_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
FPSyn_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
FPSyn_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
FPSyn_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
FPSyn_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
FPSyn_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
FPSyn_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
FPSyn_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
FPSyn_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
FPSyn_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
FPSyn_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
FPSyn_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
FPSyn_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
FPSyn_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
FPSyn_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
FPSyn_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
FPSyn_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
FPSyn_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
FPSyn_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
FPSyn_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
FPSyn_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
FPSyn_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
FPSyn_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
FPSyn_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
FPSyn_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
FPSyn_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
FPSyn_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
FPSyn_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
FPSyn_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
FPSyn_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
FPSyn_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
FPSyn_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
FPSyn_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
FPSyn_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
FPSyn_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
FPSyn_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
FPSyn_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
FPSyn_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
FPSyn_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
FPSyn_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
FPSyn_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
FPSyn_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
FPSyn_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
FPSyn_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
FPSyn_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
FPSyn_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
FPSyn_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
FPSyn_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
FPSyn_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
FPSyn_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
FPSyn_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
FPSyn_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
FPSyn_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
FPSyn_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
FPSyn_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
FPSyn_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
FPSyn_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
FPSyn_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
FPSyn_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
FPSyn_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
FPSyn_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
FPSyn_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
FPSyn_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
FPSyn_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
FPSyn_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
FPSyn_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
FPSyn_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
FPSyn_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
FPSyn_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
FPSyn_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
FPSyn_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
FPSyn_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
FPSyn_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
FPSyn_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
FPSyn_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
FPSyn_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
FPSyn_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
FPSyn_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
FPSyn_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
FPSyn_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
FPSyn_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
FPSyn_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
FPSyn_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
FPSyn_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
FPSyn_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
FPSyn_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
FPSyn_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
FPSyn_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
FPSyn_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
FPSyn_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
FPSyn_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
FPSyn_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
FPSyn_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
FPSyn_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
FPSyn_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
FPSyn_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
FPSyn_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
FPSyn_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
FPSyn_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
FPSyn_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
FPSyn_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
FPSyn_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
FPSyn_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
FPSyn_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
FPSyn_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
FPSyn_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
FPSyn_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
FPSyn_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
FPSyn_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
FPSyn_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
FPSyn_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
FPSyn_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
FPSyn_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
FPSyn_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
FPSyn_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
FPSyn_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
FPSyn_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
FPSyn_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
FPSyn_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
FPSyn_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
FPSyn_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
FPSyn_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
FPSyn_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
FPSyn_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
FPSyn_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
FPSyn_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
FPSyn_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
FPSyn_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
FPSyn_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
FPSyn_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
FPSyn_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
FPSyn_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
FPSyn_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
FPSyn_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
FPSyn_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
FPSyn_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
FPSyn_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
FPSyn_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
FPSyn_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
FPSyn_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
FPSyn_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
FPSyn_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
FPSyn_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
FPSyn_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
FPSyn_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
FPSyn_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
FPSyn_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
FPSyn_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
FPSyn_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
FPSyn_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
FPSyn_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
FPSyn_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
FPSyn_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
FPSyn_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
FPSyn_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
FPSyn_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
FPSyn_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
FPSyn_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
FPSyn_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
FPSyn_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
FPSyn_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
FPSyn_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
FPSyn_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
FPSyn_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
FPSyn_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
FPSyn_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
FPSyn_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
FPSyn_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
FPSyn_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
FPSyn_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
FPSyn_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
FPSyn_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
FPSyn_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
FPSyn_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
FPSyn_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
FPSyn_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
FPSyn_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
FPSyn_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
FPSyn_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
FPSyn_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
FPSyn_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
FPSyn_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
FPSyn_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
FPSyn_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
FPSyn_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
FPSyn_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
FPSyn_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
FPSyn_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
FPSyn_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
FPSyn_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
FPSyn_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
FPSyn_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
FPSyn_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
FPSyn_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
FPSyn_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
FPSyn_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
FPSyn_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
FPSyn_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
FPSyn_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
FPSyn_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
FPSyn_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
FPSyn_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
FPSyn_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
FPSyn_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
FPSyn_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
FPSyn_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
FPSyn_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
FPSyn_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
FPSyn_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
FPSyn_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
FPSyn_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
FPSyn_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
FPSyn_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
FPSyn_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
FPSyn_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
FPSyn_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
FPSyn_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
FPSyn_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
FPSyn_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
FPSyn_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
FPSyn_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
FPSyn_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
FPSyn_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
FPSyn_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
FPSyn_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
FPSyn_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
FPSyn_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
FPSyn_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
FPSyn_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
FPSyn_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
FPSyn_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
FPSyn_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
FPSyn_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
FPSyn_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
FPSyn_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
FPSyn_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
FPSyn_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
FPSyn_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
FPSyn_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
FPSyn_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
FPSyn_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
FPSyn_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
FPSyn_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
FPSyn_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
FPSyn_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
FPSyn_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
FPSyn_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
FPSyn_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
FPSyn_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
FPSyn_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
FPSyn_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
FPSyn_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
FPSyn_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
FPSyn_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
FPSyn_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
FPSyn_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
FPSyn_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
FPSyn_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
FPSyn_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
FPSyn_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
FPSyn_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
FPSyn_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
FPSyn_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
FPSyn_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
FPSyn_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
FPSyn_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
FPSyn_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
FPSyn_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
FPSyn_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
FPSyn_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
FPSyn_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
FPSyn_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
FPSyn_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
FPSyn_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
FPSyn_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
FPSyn_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
FPSyn_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
FPSyn_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
FPSyn_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
FPSyn_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
FPSyn_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
FPSyn_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
FPSyn_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
FPSyn_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
FPSyn_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
FPSyn_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
FPSyn_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
FPSyn_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
FPSyn_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
FPSyn_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
FPSyn_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
FPSyn_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
FPSyn_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
FPSyn_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
FPSyn_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
FPSyn_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
FPSyn_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
FPSyn_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
FPSyn_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
FPSyn_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
FPSyn_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
FPSyn_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
FPSyn_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
FPSyn_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
FPSyn_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
FPSyn_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
FPSyn_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
FPSyn_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
FPSyn_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
FPSyn_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
FPSyn_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
FPSyn_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
FPSyn_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
FPSyn_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
FPSyn_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
FPSyn_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
FPSyn_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
FPSyn_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
FPSyn_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
FPSyn_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
FPSyn_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
FPSyn_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
FPSyn_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
FPSyn_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
FPSyn_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
FPSyn_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
FPSyn_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
FPSyn_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
FPSyn_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
FPSyn_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
FPSyn_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
FPSyn_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
FPSyn_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
FPSyn_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
FPSyn_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
FPSyn_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
FPSyn_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
FPSyn_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
FPSyn_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
FPSyn_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
FPSyn_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
FPSyn_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
FPSyn_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
FPSyn_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
FPSyn_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
FPSyn_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
FPSyn_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
FPSyn_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
FPSyn_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
FPSyn_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
FPSyn_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
FPSyn_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
FPSyn_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
FPSyn_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
FPSyn_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
FPSyn_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
FPSyn_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
FPSyn_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
FPSyn_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
FPSyn_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
FPSyn_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
FPSyn_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
FPSyn_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
FPSyn_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
FPSyn_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
FPSyn_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
FPSyn_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
FPSyn_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
FPSyn_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
FPSyn_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
after = clock();
FPStime = ((double)(after - before))/CLOCKS_PER_SEC;
FPStotaltime = FPStotaltime + FPStime;

before = clock();
Manual_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
Manual_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
Manual_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
Manual_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
Manual_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
Manual_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
Manual_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
Manual_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
Manual_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
Manual_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
Manual_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
Manual_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
Manual_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
Manual_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
Manual_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
Manual_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
Manual_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
Manual_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
Manual_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
Manual_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
Manual_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
Manual_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
Manual_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
Manual_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
Manual_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
Manual_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
Manual_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
Manual_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
Manual_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
Manual_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
Manual_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
Manual_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
Manual_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
Manual_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
Manual_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
Manual_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
Manual_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
Manual_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
Manual_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
Manual_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
Manual_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
Manual_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
Manual_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
Manual_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
Manual_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
Manual_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
Manual_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
Manual_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
Manual_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
Manual_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
Manual_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
Manual_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
Manual_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
Manual_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
Manual_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
Manual_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
Manual_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
Manual_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
Manual_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
Manual_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
Manual_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
Manual_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
Manual_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
Manual_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
Manual_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
Manual_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
Manual_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
Manual_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
Manual_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
Manual_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
Manual_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
Manual_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
Manual_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
Manual_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
Manual_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
Manual_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
Manual_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
Manual_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
Manual_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
Manual_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
Manual_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
Manual_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
Manual_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
Manual_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
Manual_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
Manual_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
Manual_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
Manual_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
Manual_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
Manual_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
Manual_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
Manual_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
Manual_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
Manual_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
Manual_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
Manual_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
Manual_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
Manual_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
Manual_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
Manual_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
Manual_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
Manual_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
Manual_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
Manual_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
Manual_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
Manual_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
Manual_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
Manual_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
Manual_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
Manual_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
Manual_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
Manual_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
Manual_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
Manual_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
Manual_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
Manual_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
Manual_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
Manual_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
Manual_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
Manual_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
Manual_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
Manual_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
Manual_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
Manual_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
Manual_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
Manual_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
Manual_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
Manual_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
Manual_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
Manual_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
Manual_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
Manual_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
Manual_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
Manual_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
Manual_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
Manual_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
Manual_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
Manual_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
Manual_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
Manual_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
Manual_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
Manual_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
Manual_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
Manual_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
Manual_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
Manual_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
Manual_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
Manual_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
Manual_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
Manual_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
Manual_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
Manual_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
Manual_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
Manual_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
Manual_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
Manual_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
Manual_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
Manual_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
Manual_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
Manual_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
Manual_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
Manual_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
Manual_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
Manual_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
Manual_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
Manual_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
Manual_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
Manual_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
Manual_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
Manual_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
Manual_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
Manual_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
Manual_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
Manual_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
Manual_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
Manual_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
Manual_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
Manual_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
Manual_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
Manual_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
Manual_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
Manual_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
Manual_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
Manual_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
Manual_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
Manual_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
Manual_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
Manual_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
Manual_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
Manual_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
Manual_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
Manual_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
Manual_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
Manual_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
Manual_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
Manual_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
Manual_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
Manual_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
Manual_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
Manual_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
Manual_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
Manual_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
Manual_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
Manual_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
Manual_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
Manual_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
Manual_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
Manual_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
Manual_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
Manual_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
Manual_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
Manual_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
Manual_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
Manual_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
Manual_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
Manual_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
Manual_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
Manual_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
Manual_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
Manual_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
Manual_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
Manual_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
Manual_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
Manual_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
Manual_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
Manual_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
Manual_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
Manual_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
Manual_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
Manual_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
Manual_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
Manual_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
Manual_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
Manual_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
Manual_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
Manual_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
Manual_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
Manual_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
Manual_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
Manual_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
Manual_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
Manual_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
Manual_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
Manual_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
Manual_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
Manual_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
Manual_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
Manual_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
Manual_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
Manual_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
Manual_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
Manual_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
Manual_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
Manual_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
Manual_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
Manual_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
Manual_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
Manual_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
Manual_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
Manual_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
Manual_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
Manual_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
Manual_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
Manual_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
Manual_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
Manual_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
Manual_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
Manual_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
Manual_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
Manual_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
Manual_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
Manual_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
Manual_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
Manual_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
Manual_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
Manual_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
Manual_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
Manual_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
Manual_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
Manual_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
Manual_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
Manual_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
Manual_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
Manual_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
Manual_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
Manual_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
Manual_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
Manual_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
Manual_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
Manual_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
Manual_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
Manual_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
Manual_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
Manual_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
Manual_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
Manual_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
Manual_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
Manual_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
Manual_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
Manual_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
Manual_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
Manual_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
Manual_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
Manual_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
Manual_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
Manual_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
Manual_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
Manual_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
Manual_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
Manual_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
Manual_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
Manual_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
Manual_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
Manual_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
Manual_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
Manual_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
Manual_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
Manual_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
Manual_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
Manual_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
Manual_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
Manual_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
Manual_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
Manual_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
Manual_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
Manual_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
Manual_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
Manual_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
Manual_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
Manual_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
Manual_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
Manual_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
Manual_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
Manual_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
Manual_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
Manual_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
Manual_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
Manual_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
Manual_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
Manual_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
Manual_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
Manual_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
Manual_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
Manual_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
Manual_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
Manual_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
Manual_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
Manual_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
Manual_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
Manual_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
Manual_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
Manual_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
Manual_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
Manual_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
Manual_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
Manual_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
Manual_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
Manual_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
Manual_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
Manual_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
Manual_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
Manual_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
Manual_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
Manual_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
Manual_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
Manual_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
Manual_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
Manual_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
Manual_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
Manual_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
Manual_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
Manual_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
Manual_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
Manual_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
Manual_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
Manual_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
Manual_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
Manual_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
Manual_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
Manual_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
Manual_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
Manual_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
Manual_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
Manual_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
Manual_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
Manual_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
Manual_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
Manual_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
Manual_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
Manual_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
Manual_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
Manual_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
Manual_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
Manual_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
Manual_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
Manual_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
Manual_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
Manual_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
Manual_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
Manual_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
Manual_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
Manual_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
Manual_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
Manual_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
Manual_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
Manual_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
Manual_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
Manual_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
Manual_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
Manual_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
Manual_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
Manual_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
Manual_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
Manual_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
Manual_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
Manual_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
Manual_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
Manual_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
Manual_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
Manual_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
Manual_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
Manual_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
Manual_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
Manual_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
Manual_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
Manual_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
Manual_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
Manual_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
Manual_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
Manual_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
Manual_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
Manual_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
Manual_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
Manual_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
Manual_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
Manual_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
Manual_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
Manual_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
Manual_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
Manual_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
Manual_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
Manual_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
Manual_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
Manual_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
Manual_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
Manual_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
Manual_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
Manual_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
Manual_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
Manual_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
Manual_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
Manual_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
Manual_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
Manual_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
Manual_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
Manual_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
Manual_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
Manual_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
Manual_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
Manual_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
Manual_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
Manual_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
Manual_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
Manual_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
Manual_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
Manual_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
Manual_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
Manual_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
Manual_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
Manual_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
Manual_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
Manual_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
Manual_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
Manual_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
Manual_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
Manual_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
Manual_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
Manual_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
Manual_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
Manual_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
Manual_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
Manual_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
Manual_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
Manual_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
Manual_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
Manual_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
Manual_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
Manual_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
Manual_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
Manual_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
Manual_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
Manual_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
Manual_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
Manual_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
Manual_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
Manual_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
Manual_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
Manual_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
Manual_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
Manual_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
Manual_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
Manual_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
Manual_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
Manual_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
Manual_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
Manual_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
Manual_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
Manual_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
Manual_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
Manual_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
Manual_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
Manual_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
Manual_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
Manual_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
Manual_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
Manual_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
Manual_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
Manual_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
Manual_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
Manual_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
Manual_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
Manual_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
Manual_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
Manual_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
Manual_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
Manual_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
Manual_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
Manual_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
Manual_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
Manual_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
Manual_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
Manual_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
Manual_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
Manual_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
Manual_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
Manual_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
Manual_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
Manual_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
Manual_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
Manual_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
Manual_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
Manual_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
Manual_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
Manual_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
Manual_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
Manual_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
Manual_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
Manual_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
Manual_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
Manual_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
Manual_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
Manual_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
Manual_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
Manual_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
Manual_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
Manual_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
Manual_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
Manual_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
Manual_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
Manual_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
Manual_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
Manual_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
Manual_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
Manual_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
Manual_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
Manual_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
Manual_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
Manual_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
Manual_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
Manual_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
Manual_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
Manual_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
Manual_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
Manual_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
Manual_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
Manual_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
Manual_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
Manual_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
Manual_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
Manual_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
Manual_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
Manual_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
Manual_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
Manual_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
Manual_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
Manual_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
Manual_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
Manual_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
Manual_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
Manual_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
Manual_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
Manual_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
Manual_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
Manual_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
Manual_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
Manual_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
Manual_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
Manual_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
Manual_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
Manual_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
Manual_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
Manual_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
Manual_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
Manual_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
Manual_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
Manual_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
Manual_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
Manual_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
Manual_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
Manual_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
Manual_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
Manual_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
Manual_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
Manual_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
Manual_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
Manual_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
Manual_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
Manual_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
Manual_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
Manual_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
Manual_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
Manual_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
Manual_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
Manual_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
Manual_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
Manual_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
Manual_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
Manual_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
Manual_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
Manual_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
Manual_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
Manual_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
Manual_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
Manual_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
Manual_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
Manual_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
Manual_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
Manual_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
Manual_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
Manual_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
Manual_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
Manual_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
Manual_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
Manual_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
Manual_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
Manual_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
Manual_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
Manual_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
Manual_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
Manual_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
Manual_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
Manual_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
Manual_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
Manual_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
Manual_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
Manual_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
Manual_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
Manual_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
Manual_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
Manual_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
Manual_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
Manual_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
Manual_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
Manual_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
Manual_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
Manual_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
Manual_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
Manual_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
Manual_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
Manual_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
Manual_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
Manual_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
Manual_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
Manual_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
Manual_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
Manual_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
Manual_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
Manual_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
Manual_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
Manual_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
Manual_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
Manual_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
Manual_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
Manual_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
Manual_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
Manual_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
Manual_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
Manual_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
Manual_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
Manual_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
Manual_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
Manual_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
Manual_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
Manual_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
Manual_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
Manual_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
Manual_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
Manual_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
Manual_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
Manual_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
Manual_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
Manual_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
Manual_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
Manual_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
Manual_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
Manual_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
Manual_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
Manual_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
Manual_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
Manual_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
Manual_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
Manual_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
Manual_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
Manual_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
Manual_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
Manual_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
Manual_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
Manual_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
Manual_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
Manual_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
Manual_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
Manual_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
Manual_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
Manual_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
Manual_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
Manual_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
Manual_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
Manual_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
Manual_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
Manual_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
Manual_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
Manual_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
Manual_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
Manual_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
Manual_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
Manual_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
Manual_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
Manual_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
Manual_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
Manual_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
Manual_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
Manual_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
Manual_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
Manual_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
Manual_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
Manual_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
Manual_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
Manual_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
Manual_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
Manual_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
Manual_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
Manual_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
Manual_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
Manual_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
Manual_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
Manual_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
Manual_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
Manual_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
Manual_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
Manual_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
Manual_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
Manual_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
Manual_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
Manual_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
Manual_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
Manual_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
Manual_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
Manual_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
Manual_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
Manual_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
Manual_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
Manual_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
Manual_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
Manual_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
Manual_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
Manual_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
Manual_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
Manual_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
Manual_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
Manual_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
Manual_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
Manual_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
Manual_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
Manual_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
Manual_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
Manual_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
Manual_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
Manual_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
Manual_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
Manual_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
Manual_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
Manual_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
Manual_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
Manual_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
Manual_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
Manual_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
Manual_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
Manual_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
Manual_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
Manual_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
Manual_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
Manual_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
Manual_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
Manual_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
Manual_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
Manual_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
Manual_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
Manual_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
Manual_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
Manual_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
Manual_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
Manual_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
Manual_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
Manual_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
Manual_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
Manual_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
Manual_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
Manual_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
Manual_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
Manual_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
Manual_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
Manual_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
Manual_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
Manual_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
Manual_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
Manual_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
Manual_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
Manual_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
Manual_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
Manual_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
Manual_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
Manual_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
Manual_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
Manual_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
Manual_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
Manual_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
Manual_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
Manual_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
Manual_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
Manual_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
Manual_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
Manual_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
Manual_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
Manual_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
Manual_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
Manual_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
Manual_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
Manual_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
Manual_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
Manual_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
Manual_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
Manual_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
Manual_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
Manual_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
Manual_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
Manual_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
Manual_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
Manual_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
Manual_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
Manual_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
Manual_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
Manual_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
Manual_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
Manual_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
Manual_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
Manual_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
Manual_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
Manual_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
Manual_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
Manual_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
Manual_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
Manual_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
Manual_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
Manual_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
Manual_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
Manual_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
Manual_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
Manual_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
Manual_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
Manual_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
Manual_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
Manual_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
Manual_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
Manual_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
Manual_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
Manual_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
Manual_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
Manual_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
Manual_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
Manual_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
Manual_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
Manual_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
Manual_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
Manual_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
Manual_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
Manual_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
Manual_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
Manual_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
Manual_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
Manual_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
Manual_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
Manual_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
Manual_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
Manual_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
Manual_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
Manual_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
Manual_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
Manual_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
Manual_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
Manual_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
Manual_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
Manual_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
Manual_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
Manual_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
Manual_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
Manual_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
Manual_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
Manual_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
Manual_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
Manual_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
Manual_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
Manual_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
Manual_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
Manual_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
Manual_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
Manual_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
Manual_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
Manual_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
Manual_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
Manual_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
Manual_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
Manual_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
Manual_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
Manual_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
Manual_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
Manual_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
Manual_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
Manual_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
Manual_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
Manual_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
Manual_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
Manual_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
Manual_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
Manual_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
Manual_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
Manual_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
Manual_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
Manual_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
Manual_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
Manual_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
Manual_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
Manual_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
Manual_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
Manual_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
Manual_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
Manual_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
Manual_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
Manual_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
Manual_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
Manual_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
Manual_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
Manual_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
Manual_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
Manual_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
Manual_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
Manual_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
Manual_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
Manual_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
Manual_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
Manual_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
Manual_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
Manual_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
Manual_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
Manual_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
Manual_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
Manual_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
Manual_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
Manual_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
Manual_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
Manual_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
Manual_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
Manual_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
Manual_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
Manual_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
Manual_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
Manual_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
Manual_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
Manual_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
Manual_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
Manual_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
Manual_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
Manual_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
Manual_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
Manual_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
Manual_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
Manual_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
Manual_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
Manual_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
Manual_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
Manual_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&err);
after = clock();
Manualtime = ((double)(after - before))/CLOCKS_PER_SEC;
Manualtotaltime = Manualtotaltime + Manualtime;

before = clock();
IA_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&herr,&lerr);
IA_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&herr,&lerr);
IA_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&herr,&lerr);
IA_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&herr,&lerr);
IA_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&herr,&lerr);
IA_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&herr,&lerr);
IA_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&herr,&lerr);
IA_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&herr,&lerr);
IA_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&herr,&lerr);
IA_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&herr,&lerr);
IA_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&herr,&lerr);
IA_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&herr,&lerr);
IA_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&herr,&lerr);
IA_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&herr,&lerr);
IA_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&herr,&lerr);
IA_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&herr,&lerr);
IA_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&herr,&lerr);
IA_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&herr,&lerr);
IA_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&herr,&lerr);
IA_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&herr,&lerr);
IA_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&herr,&lerr);
IA_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&herr,&lerr);
IA_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&herr,&lerr);
IA_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&herr,&lerr);
IA_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&herr,&lerr);
IA_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&herr,&lerr);
IA_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&herr,&lerr);
IA_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&herr,&lerr);
IA_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&herr,&lerr);
IA_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&herr,&lerr);
IA_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&herr,&lerr);
IA_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&herr,&lerr);
IA_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&herr,&lerr);
IA_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&herr,&lerr);
IA_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&herr,&lerr);
IA_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&herr,&lerr);
IA_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&herr,&lerr);
IA_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&herr,&lerr);
IA_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&herr,&lerr);
IA_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&herr,&lerr);
IA_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&herr,&lerr);
IA_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&herr,&lerr);
IA_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&herr,&lerr);
IA_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&herr,&lerr);
IA_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&herr,&lerr);
IA_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&herr,&lerr);
IA_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&herr,&lerr);
IA_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&herr,&lerr);
IA_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&herr,&lerr);
IA_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&herr,&lerr);
IA_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&herr,&lerr);
IA_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&herr,&lerr);
IA_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&herr,&lerr);
IA_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&herr,&lerr);
IA_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&herr,&lerr);
IA_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&herr,&lerr);
IA_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&herr,&lerr);
IA_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&herr,&lerr);
IA_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&herr,&lerr);
IA_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&herr,&lerr);
IA_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&herr,&lerr);
IA_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&herr,&lerr);
IA_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&herr,&lerr);
IA_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&herr,&lerr);
IA_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&herr,&lerr);
IA_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&herr,&lerr);
IA_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&herr,&lerr);
IA_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&herr,&lerr);
IA_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&herr,&lerr);
IA_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&herr,&lerr);
IA_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&herr,&lerr);
IA_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&herr,&lerr);
IA_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&herr,&lerr);
IA_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&herr,&lerr);
IA_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&herr,&lerr);
IA_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&herr,&lerr);
IA_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&herr,&lerr);
IA_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&herr,&lerr);
IA_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&herr,&lerr);
IA_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&herr,&lerr);
IA_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&herr,&lerr);
IA_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&herr,&lerr);
IA_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&herr,&lerr);
IA_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&herr,&lerr);
IA_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&herr,&lerr);
IA_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&herr,&lerr);
IA_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&herr,&lerr);
IA_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&herr,&lerr);
IA_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&herr,&lerr);
IA_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&herr,&lerr);
IA_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&herr,&lerr);
IA_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&herr,&lerr);
IA_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&herr,&lerr);
IA_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&herr,&lerr);
IA_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&herr,&lerr);
IA_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&herr,&lerr);
IA_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&herr,&lerr);
IA_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&herr,&lerr);
IA_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&herr,&lerr);
IA_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&herr,&lerr);
IA_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&herr,&lerr);
IA_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&herr,&lerr);
IA_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&herr,&lerr);
IA_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&herr,&lerr);
IA_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&herr,&lerr);
IA_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&herr,&lerr);
IA_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&herr,&lerr);
IA_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&herr,&lerr);
IA_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&herr,&lerr);
IA_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&herr,&lerr);
IA_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&herr,&lerr);
IA_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&herr,&lerr);
IA_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&herr,&lerr);
IA_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&herr,&lerr);
IA_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&herr,&lerr);
IA_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&herr,&lerr);
IA_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&herr,&lerr);
IA_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&herr,&lerr);
IA_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&herr,&lerr);
IA_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&herr,&lerr);
IA_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&herr,&lerr);
IA_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&herr,&lerr);
IA_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&herr,&lerr);
IA_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&herr,&lerr);
IA_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&herr,&lerr);
IA_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&herr,&lerr);
IA_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&herr,&lerr);
IA_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&herr,&lerr);
IA_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&herr,&lerr);
IA_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&herr,&lerr);
IA_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&herr,&lerr);
IA_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&herr,&lerr);
IA_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&herr,&lerr);
IA_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&herr,&lerr);
IA_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&herr,&lerr);
IA_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&herr,&lerr);
IA_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&herr,&lerr);
IA_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&herr,&lerr);
IA_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&herr,&lerr);
IA_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&herr,&lerr);
IA_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&herr,&lerr);
IA_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&herr,&lerr);
IA_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&herr,&lerr);
IA_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&herr,&lerr);
IA_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&herr,&lerr);
IA_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&herr,&lerr);
IA_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&herr,&lerr);
IA_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&herr,&lerr);
IA_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&herr,&lerr);
IA_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&herr,&lerr);
IA_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&herr,&lerr);
IA_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&herr,&lerr);
IA_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&herr,&lerr);
IA_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&herr,&lerr);
IA_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&herr,&lerr);
IA_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&herr,&lerr);
IA_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&herr,&lerr);
IA_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&herr,&lerr);
IA_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&herr,&lerr);
IA_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&herr,&lerr);
IA_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&herr,&lerr);
IA_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&herr,&lerr);
IA_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&herr,&lerr);
IA_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&herr,&lerr);
IA_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&herr,&lerr);
IA_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&herr,&lerr);
IA_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&herr,&lerr);
IA_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&herr,&lerr);
IA_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&herr,&lerr);
IA_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&herr,&lerr);
IA_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&herr,&lerr);
IA_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&herr,&lerr);
IA_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&herr,&lerr);
IA_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&herr,&lerr);
IA_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&herr,&lerr);
IA_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&herr,&lerr);
IA_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&herr,&lerr);
IA_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&herr,&lerr);
IA_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&herr,&lerr);
IA_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&herr,&lerr);
IA_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&herr,&lerr);
IA_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&herr,&lerr);
IA_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&herr,&lerr);
IA_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&herr,&lerr);
IA_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&herr,&lerr);
IA_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&herr,&lerr);
IA_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&herr,&lerr);
IA_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&herr,&lerr);
IA_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&herr,&lerr);
IA_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&herr,&lerr);
IA_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&herr,&lerr);
IA_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&herr,&lerr);
IA_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&herr,&lerr);
IA_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&herr,&lerr);
IA_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&herr,&lerr);
IA_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&herr,&lerr);
IA_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&herr,&lerr);
IA_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&herr,&lerr);
IA_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&herr,&lerr);
IA_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&herr,&lerr);
IA_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&herr,&lerr);
IA_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&herr,&lerr);
IA_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&herr,&lerr);
IA_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&herr,&lerr);
IA_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&herr,&lerr);
IA_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&herr,&lerr);
IA_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&herr,&lerr);
IA_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&herr,&lerr);
IA_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&herr,&lerr);
IA_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&herr,&lerr);
IA_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&herr,&lerr);
IA_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&herr,&lerr);
IA_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&herr,&lerr);
IA_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&herr,&lerr);
IA_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&herr,&lerr);
IA_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&herr,&lerr);
IA_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&herr,&lerr);
IA_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&herr,&lerr);
IA_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&herr,&lerr);
IA_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&herr,&lerr);
IA_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&herr,&lerr);
IA_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&herr,&lerr);
IA_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&herr,&lerr);
IA_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&herr,&lerr);
IA_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&herr,&lerr);
IA_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&herr,&lerr);
IA_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&herr,&lerr);
IA_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&herr,&lerr);
IA_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&herr,&lerr);
IA_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&herr,&lerr);
IA_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&herr,&lerr);
IA_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&herr,&lerr);
IA_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&herr,&lerr);
IA_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&herr,&lerr);
IA_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&herr,&lerr);
IA_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&herr,&lerr);
IA_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&herr,&lerr);
IA_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&herr,&lerr);
IA_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&herr,&lerr);
IA_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&herr,&lerr);
IA_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&herr,&lerr);
IA_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&herr,&lerr);
IA_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&herr,&lerr);
IA_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&herr,&lerr);
IA_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&herr,&lerr);
IA_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&herr,&lerr);
IA_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&herr,&lerr);
IA_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&herr,&lerr);
IA_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&herr,&lerr);
IA_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&herr,&lerr);
IA_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&herr,&lerr);
IA_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&herr,&lerr);
IA_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&herr,&lerr);
IA_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&herr,&lerr);
IA_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&herr,&lerr);
IA_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&herr,&lerr);
IA_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&herr,&lerr);
IA_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&herr,&lerr);
IA_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&herr,&lerr);
IA_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&herr,&lerr);
IA_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&herr,&lerr);
IA_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&herr,&lerr);
IA_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&herr,&lerr);
IA_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&herr,&lerr);
IA_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&herr,&lerr);
IA_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&herr,&lerr);
IA_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&herr,&lerr);
IA_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&herr,&lerr);
IA_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&herr,&lerr);
IA_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&herr,&lerr);
IA_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&herr,&lerr);
IA_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&herr,&lerr);
IA_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&herr,&lerr);
IA_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&herr,&lerr);
IA_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&herr,&lerr);
IA_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&herr,&lerr);
IA_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&herr,&lerr);
IA_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&herr,&lerr);
IA_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&herr,&lerr);
IA_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&herr,&lerr);
IA_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&herr,&lerr);
IA_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&herr,&lerr);
IA_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&herr,&lerr);
IA_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&herr,&lerr);
IA_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&herr,&lerr);
IA_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&herr,&lerr);
IA_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&herr,&lerr);
IA_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&herr,&lerr);
IA_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&herr,&lerr);
IA_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&herr,&lerr);
IA_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&herr,&lerr);
IA_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&herr,&lerr);
IA_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&herr,&lerr);
IA_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&herr,&lerr);
IA_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&herr,&lerr);
IA_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&herr,&lerr);
IA_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&herr,&lerr);
IA_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&herr,&lerr);
IA_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&herr,&lerr);
IA_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&herr,&lerr);
IA_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&herr,&lerr);
IA_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&herr,&lerr);
IA_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&herr,&lerr);
IA_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&herr,&lerr);
IA_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&herr,&lerr);
IA_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&herr,&lerr);
IA_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&herr,&lerr);
IA_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&herr,&lerr);
IA_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&herr,&lerr);
IA_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&herr,&lerr);
IA_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&herr,&lerr);
IA_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&herr,&lerr);
IA_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&herr,&lerr);
IA_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&herr,&lerr);
IA_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&herr,&lerr);
IA_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&herr,&lerr);
IA_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&herr,&lerr);
IA_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&herr,&lerr);
IA_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&herr,&lerr);
IA_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&herr,&lerr);
IA_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&herr,&lerr);
IA_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&herr,&lerr);
IA_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&herr,&lerr);
IA_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&herr,&lerr);
IA_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&herr,&lerr);
IA_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&herr,&lerr);
IA_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&herr,&lerr);
IA_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&herr,&lerr);
IA_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&herr,&lerr);
IA_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&herr,&lerr);
IA_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&herr,&lerr);
IA_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&herr,&lerr);
IA_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&herr,&lerr);
IA_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&herr,&lerr);
IA_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&herr,&lerr);
IA_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&herr,&lerr);
IA_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&herr,&lerr);
IA_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&herr,&lerr);
IA_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&herr,&lerr);
IA_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&herr,&lerr);
IA_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&herr,&lerr);
IA_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&herr,&lerr);
IA_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&herr,&lerr);
IA_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&herr,&lerr);
IA_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&herr,&lerr);
IA_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&herr,&lerr);
IA_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&herr,&lerr);
IA_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&herr,&lerr);
IA_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&herr,&lerr);
IA_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&herr,&lerr);
IA_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&herr,&lerr);
IA_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&herr,&lerr);
IA_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&herr,&lerr);
IA_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&herr,&lerr);
IA_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&herr,&lerr);
IA_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&herr,&lerr);
IA_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&herr,&lerr);
IA_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&herr,&lerr);
IA_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&herr,&lerr);
IA_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&herr,&lerr);
IA_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&herr,&lerr);
IA_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&herr,&lerr);
IA_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&herr,&lerr);
IA_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&herr,&lerr);
IA_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&herr,&lerr);
IA_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&herr,&lerr);
IA_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&herr,&lerr);
IA_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&herr,&lerr);
IA_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&herr,&lerr);
IA_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&herr,&lerr);
IA_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&herr,&lerr);
IA_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&herr,&lerr);
IA_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&herr,&lerr);
IA_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&herr,&lerr);
IA_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&herr,&lerr);
IA_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&herr,&lerr);
IA_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&herr,&lerr);
IA_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&herr,&lerr);
IA_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&herr,&lerr);
IA_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&herr,&lerr);
IA_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&herr,&lerr);
IA_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&herr,&lerr);
IA_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&herr,&lerr);
IA_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&herr,&lerr);
IA_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&herr,&lerr);
IA_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&herr,&lerr);
IA_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&herr,&lerr);
IA_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&herr,&lerr);
IA_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&herr,&lerr);
IA_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&herr,&lerr);
IA_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&herr,&lerr);
IA_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&herr,&lerr);
IA_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&herr,&lerr);
IA_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&herr,&lerr);
IA_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&herr,&lerr);
IA_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&herr,&lerr);
IA_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&herr,&lerr);
IA_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&herr,&lerr);
IA_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&herr,&lerr);
IA_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&herr,&lerr);
IA_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&herr,&lerr);
IA_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&herr,&lerr);
IA_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&herr,&lerr);
IA_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&herr,&lerr);
IA_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&herr,&lerr);
IA_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&herr,&lerr);
IA_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&herr,&lerr);
IA_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&herr,&lerr);
IA_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&herr,&lerr);
IA_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&herr,&lerr);
IA_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&herr,&lerr);
IA_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&herr,&lerr);
IA_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&herr,&lerr);
IA_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&herr,&lerr);
IA_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&herr,&lerr);
IA_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&herr,&lerr);
IA_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&herr,&lerr);
IA_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&herr,&lerr);
IA_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&herr,&lerr);
IA_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&herr,&lerr);
IA_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&herr,&lerr);
IA_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&herr,&lerr);
IA_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&herr,&lerr);
IA_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&herr,&lerr);
IA_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&herr,&lerr);
IA_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&herr,&lerr);
IA_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&herr,&lerr);
IA_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&herr,&lerr);
IA_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&herr,&lerr);
IA_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&herr,&lerr);
IA_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&herr,&lerr);
IA_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&herr,&lerr);
IA_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&herr,&lerr);
IA_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&herr,&lerr);
IA_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&herr,&lerr);
IA_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&herr,&lerr);
IA_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&herr,&lerr);
IA_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&herr,&lerr);
IA_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&herr,&lerr);
IA_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&herr,&lerr);
IA_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&herr,&lerr);
IA_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&herr,&lerr);
IA_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&herr,&lerr);
IA_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&herr,&lerr);
IA_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&herr,&lerr);
IA_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&herr,&lerr);
IA_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&herr,&lerr);
IA_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&herr,&lerr);
IA_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&herr,&lerr);
IA_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&herr,&lerr);
IA_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&herr,&lerr);
IA_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&herr,&lerr);
IA_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&herr,&lerr);
IA_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&herr,&lerr);
IA_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&herr,&lerr);
IA_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&herr,&lerr);
IA_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&herr,&lerr);
IA_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&herr,&lerr);
IA_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&herr,&lerr);
IA_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&herr,&lerr);
IA_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&herr,&lerr);
IA_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&herr,&lerr);
IA_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&herr,&lerr);
IA_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&herr,&lerr);
IA_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&herr,&lerr);
IA_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&herr,&lerr);
IA_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&herr,&lerr);
IA_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&herr,&lerr);
IA_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&herr,&lerr);
IA_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&herr,&lerr);
IA_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&herr,&lerr);
IA_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&herr,&lerr);
IA_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&herr,&lerr);
IA_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&herr,&lerr);
IA_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&herr,&lerr);
IA_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&herr,&lerr);
IA_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&herr,&lerr);
IA_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&herr,&lerr);
IA_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&herr,&lerr);
IA_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&herr,&lerr);
IA_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&herr,&lerr);
IA_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&herr,&lerr);
IA_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&herr,&lerr);
IA_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&herr,&lerr);
IA_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&herr,&lerr);
IA_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&herr,&lerr);
IA_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&herr,&lerr);
IA_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&herr,&lerr);
IA_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&herr,&lerr);
IA_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&herr,&lerr);
IA_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&herr,&lerr);
IA_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&herr,&lerr);
IA_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&herr,&lerr);
IA_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&herr,&lerr);
IA_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&herr,&lerr);
IA_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&herr,&lerr);
IA_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&herr,&lerr);
IA_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&herr,&lerr);
IA_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&herr,&lerr);
IA_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&herr,&lerr);
IA_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&herr,&lerr);
IA_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&herr,&lerr);
IA_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&herr,&lerr);
IA_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&herr,&lerr);
IA_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&herr,&lerr);
IA_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&herr,&lerr);
IA_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&herr,&lerr);
IA_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&herr,&lerr);
IA_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&herr,&lerr);
IA_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&herr,&lerr);
IA_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&herr,&lerr);
IA_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&herr,&lerr);
IA_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&herr,&lerr);
IA_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&herr,&lerr);
IA_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&herr,&lerr);
IA_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&herr,&lerr);
IA_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&herr,&lerr);
IA_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&herr,&lerr);
IA_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&herr,&lerr);
IA_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&herr,&lerr);
IA_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&herr,&lerr);
IA_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&herr,&lerr);
IA_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&herr,&lerr);
IA_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&herr,&lerr);
IA_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&herr,&lerr);
IA_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&herr,&lerr);
IA_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&herr,&lerr);
IA_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&herr,&lerr);
IA_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&herr,&lerr);
IA_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&herr,&lerr);
IA_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&herr,&lerr);
IA_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&herr,&lerr);
IA_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&herr,&lerr);
IA_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&herr,&lerr);
IA_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&herr,&lerr);
IA_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&herr,&lerr);
IA_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&herr,&lerr);
IA_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&herr,&lerr);
IA_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&herr,&lerr);
IA_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&herr,&lerr);
IA_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&herr,&lerr);
IA_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&herr,&lerr);
IA_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&herr,&lerr);
IA_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&herr,&lerr);
IA_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&herr,&lerr);
IA_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&herr,&lerr);
IA_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&herr,&lerr);
IA_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&herr,&lerr);
IA_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&herr,&lerr);
IA_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&herr,&lerr);
IA_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&herr,&lerr);
IA_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&herr,&lerr);
IA_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&herr,&lerr);
IA_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&herr,&lerr);
IA_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&herr,&lerr);
IA_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&herr,&lerr);
IA_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&herr,&lerr);
IA_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&herr,&lerr);
IA_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&herr,&lerr);
IA_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&herr,&lerr);
IA_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&herr,&lerr);
IA_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&herr,&lerr);
IA_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&herr,&lerr);
IA_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&herr,&lerr);
IA_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&herr,&lerr);
IA_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&herr,&lerr);
IA_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&herr,&lerr);
IA_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&herr,&lerr);
IA_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&herr,&lerr);
IA_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&herr,&lerr);
IA_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&herr,&lerr);
IA_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&herr,&lerr);
IA_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&herr,&lerr);
IA_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&herr,&lerr);
IA_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&herr,&lerr);
IA_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&herr,&lerr);
IA_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&herr,&lerr);
IA_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&herr,&lerr);
IA_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&herr,&lerr);
IA_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&herr,&lerr);
IA_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&herr,&lerr);
IA_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&herr,&lerr);
IA_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&herr,&lerr);
IA_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&herr,&lerr);
IA_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&herr,&lerr);
IA_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&herr,&lerr);
IA_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&herr,&lerr);
IA_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&herr,&lerr);
IA_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&herr,&lerr);
IA_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&herr,&lerr);
IA_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&herr,&lerr);
IA_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&herr,&lerr);
IA_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&herr,&lerr);
IA_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&herr,&lerr);
IA_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&herr,&lerr);
IA_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&herr,&lerr);
IA_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&herr,&lerr);
IA_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&herr,&lerr);
IA_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&herr,&lerr);
IA_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&herr,&lerr);
IA_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&herr,&lerr);
IA_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&herr,&lerr);
IA_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&herr,&lerr);
IA_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&herr,&lerr);
IA_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&herr,&lerr);
IA_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&herr,&lerr);
IA_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&herr,&lerr);
IA_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&herr,&lerr);
IA_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&herr,&lerr);
IA_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&herr,&lerr);
IA_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&herr,&lerr);
IA_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&herr,&lerr);
IA_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&herr,&lerr);
IA_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&herr,&lerr);
IA_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&herr,&lerr);
IA_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&herr,&lerr);
IA_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&herr,&lerr);
IA_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&herr,&lerr);
IA_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&herr,&lerr);
IA_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&herr,&lerr);
IA_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&herr,&lerr);
IA_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&herr,&lerr);
IA_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&herr,&lerr);
IA_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&herr,&lerr);
IA_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&herr,&lerr);
IA_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&herr,&lerr);
IA_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&herr,&lerr);
IA_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&herr,&lerr);
IA_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&herr,&lerr);
IA_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&herr,&lerr);
IA_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&herr,&lerr);
IA_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&herr,&lerr);
IA_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&herr,&lerr);
IA_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&herr,&lerr);
IA_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&herr,&lerr);
IA_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&herr,&lerr);
IA_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&herr,&lerr);
IA_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&herr,&lerr);
IA_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&herr,&lerr);
IA_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&herr,&lerr);
IA_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&herr,&lerr);
IA_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&herr,&lerr);
IA_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&herr,&lerr);
IA_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&herr,&lerr);
IA_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&herr,&lerr);
IA_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&herr,&lerr);
IA_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&herr,&lerr);
IA_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&herr,&lerr);
IA_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&herr,&lerr);
IA_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&herr,&lerr);
IA_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&herr,&lerr);
IA_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&herr,&lerr);
IA_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&herr,&lerr);
IA_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&herr,&lerr);
IA_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&herr,&lerr);
IA_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&herr,&lerr);
IA_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&herr,&lerr);
IA_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&herr,&lerr);
IA_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&herr,&lerr);
IA_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&herr,&lerr);
IA_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&herr,&lerr);
IA_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&herr,&lerr);
IA_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&herr,&lerr);
IA_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&herr,&lerr);
IA_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&herr,&lerr);
IA_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&herr,&lerr);
IA_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&herr,&lerr);
IA_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&herr,&lerr);
IA_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&herr,&lerr);
IA_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&herr,&lerr);
IA_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&herr,&lerr);
IA_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&herr,&lerr);
IA_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&herr,&lerr);
IA_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&herr,&lerr);
IA_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&herr,&lerr);
IA_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&herr,&lerr);
IA_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&herr,&lerr);
IA_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&herr,&lerr);
IA_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&herr,&lerr);
IA_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&herr,&lerr);
IA_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&herr,&lerr);
IA_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&herr,&lerr);
IA_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&herr,&lerr);
IA_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&herr,&lerr);
IA_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&herr,&lerr);
IA_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&herr,&lerr);
IA_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&herr,&lerr);
IA_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&herr,&lerr);
IA_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&herr,&lerr);
IA_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&herr,&lerr);
IA_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&herr,&lerr);
IA_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&herr,&lerr);
IA_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&herr,&lerr);
IA_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&herr,&lerr);
IA_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&herr,&lerr);
IA_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&herr,&lerr);
IA_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&herr,&lerr);
IA_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&herr,&lerr);
IA_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&herr,&lerr);
IA_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&herr,&lerr);
IA_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&herr,&lerr);
IA_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&herr,&lerr);
IA_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&herr,&lerr);
IA_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&herr,&lerr);
IA_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&herr,&lerr);
IA_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&herr,&lerr);
IA_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&herr,&lerr);
IA_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&herr,&lerr);
IA_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&herr,&lerr);
IA_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&herr,&lerr);
IA_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&herr,&lerr);
IA_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&herr,&lerr);
IA_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&herr,&lerr);
IA_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&herr,&lerr);
IA_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&herr,&lerr);
IA_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&herr,&lerr);
IA_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&herr,&lerr);
IA_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&herr,&lerr);
IA_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&herr,&lerr);
IA_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&herr,&lerr);
IA_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&herr,&lerr);
IA_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&herr,&lerr);
IA_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&herr,&lerr);
IA_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&herr,&lerr);
IA_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&herr,&lerr);
IA_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&herr,&lerr);
IA_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&herr,&lerr);
IA_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&herr,&lerr);
IA_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&herr,&lerr);
IA_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&herr,&lerr);
IA_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&herr,&lerr);
IA_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&herr,&lerr);
IA_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&herr,&lerr);
IA_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&herr,&lerr);
IA_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&herr,&lerr);
IA_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&herr,&lerr);
IA_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&herr,&lerr);
IA_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&herr,&lerr);
IA_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&herr,&lerr);
IA_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&herr,&lerr);
IA_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&herr,&lerr);
IA_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&herr,&lerr);
IA_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&herr,&lerr);
IA_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&herr,&lerr);
IA_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&herr,&lerr);
IA_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&herr,&lerr);
IA_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&herr,&lerr);
IA_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&herr,&lerr);
IA_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&herr,&lerr);
IA_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&herr,&lerr);
IA_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&herr,&lerr);
IA_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&herr,&lerr);
IA_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&herr,&lerr);
IA_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&herr,&lerr);
IA_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&herr,&lerr);
IA_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&herr,&lerr);
IA_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&herr,&lerr);
IA_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&herr,&lerr);
IA_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&herr,&lerr);
IA_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&herr,&lerr);
IA_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&herr,&lerr);
IA_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&herr,&lerr);
IA_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&herr,&lerr);
IA_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&herr,&lerr);
IA_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&herr,&lerr);
IA_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&herr,&lerr);
IA_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&herr,&lerr);
IA_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&herr,&lerr);
IA_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&herr,&lerr);
IA_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&herr,&lerr);
IA_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&herr,&lerr);
IA_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&herr,&lerr);
IA_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&herr,&lerr);
IA_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&herr,&lerr);
IA_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&herr,&lerr);
IA_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&herr,&lerr);
IA_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&herr,&lerr);
IA_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&herr,&lerr);
IA_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&herr,&lerr);
IA_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&herr,&lerr);
IA_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&herr,&lerr);
IA_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&herr,&lerr);
IA_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&herr,&lerr);
IA_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&herr,&lerr);
IA_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&herr,&lerr);
IA_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&herr,&lerr);
IA_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&herr,&lerr);
IA_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&herr,&lerr);
IA_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&herr,&lerr);
IA_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&herr,&lerr);
IA_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&herr,&lerr);
IA_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&herr,&lerr);
IA_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&herr,&lerr);
IA_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&herr,&lerr);
IA_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&herr,&lerr);
IA_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&herr,&lerr);
IA_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&herr,&lerr);
IA_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&herr,&lerr);
IA_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&herr,&lerr);
IA_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&herr,&lerr);
IA_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&herr,&lerr);
IA_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&herr,&lerr);
IA_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&herr,&lerr);
IA_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&herr,&lerr);
IA_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&herr,&lerr);
IA_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&herr,&lerr);
IA_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&herr,&lerr);
IA_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&herr,&lerr);
IA_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&herr,&lerr);
IA_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&herr,&lerr);
IA_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&herr,&lerr);
IA_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&herr,&lerr);
IA_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&herr,&lerr);
IA_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&herr,&lerr);
IA_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&herr,&lerr);
IA_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&herr,&lerr);
IA_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&herr,&lerr);
IA_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&herr,&lerr);
IA_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&herr,&lerr);
IA_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&herr,&lerr);
IA_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&herr,&lerr);
IA_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&herr,&lerr);
IA_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&herr,&lerr);
IA_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&herr,&lerr);
IA_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&herr,&lerr);
IA_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&herr,&lerr);
IA_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&herr,&lerr);
IA_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&herr,&lerr);
IA_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&herr,&lerr);
IA_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&herr,&lerr);
IA_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&herr,&lerr);
IA_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&herr,&lerr);
IA_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&herr,&lerr);
IA_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&herr,&lerr);
IA_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&herr,&lerr);
IA_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&herr,&lerr);
IA_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&herr,&lerr);
IA_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&herr,&lerr);
IA_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&herr,&lerr);
IA_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&herr,&lerr);
IA_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&herr,&lerr);
IA_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&herr,&lerr);
IA_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&herr,&lerr);
IA_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&herr,&lerr);
IA_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&herr,&lerr);
IA_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&herr,&lerr);
IA_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&herr,&lerr);
IA_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&herr,&lerr);
IA_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&herr,&lerr);
IA_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&herr,&lerr);
IA_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&herr,&lerr);
IA_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&herr,&lerr);
IA_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&herr,&lerr);
IA_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&herr,&lerr);
IA_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&herr,&lerr);
IA_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&herr,&lerr);
IA_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&herr,&lerr);
IA_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&herr,&lerr);
IA_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&herr,&lerr);
IA_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&herr,&lerr);
IA_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&herr,&lerr);
IA_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&herr,&lerr);
IA_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&herr,&lerr);
IA_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&herr,&lerr);
IA_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&herr,&lerr);
IA_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&herr,&lerr);
IA_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&herr,&lerr);
IA_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&herr,&lerr);
IA_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&herr,&lerr);
IA_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&herr,&lerr);
IA_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&herr,&lerr);
IA_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&herr,&lerr);
IA_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&herr,&lerr);
IA_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&herr,&lerr);
IA_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&herr,&lerr);
IA_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&herr,&lerr);
IA_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&herr,&lerr);
IA_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&herr,&lerr);
IA_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&herr,&lerr);
IA_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&herr,&lerr);
IA_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&herr,&lerr);
IA_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&herr,&lerr);
IA_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&herr,&lerr);
IA_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&herr,&lerr);
IA_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&herr,&lerr);
IA_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&herr,&lerr);
IA_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&herr,&lerr);
IA_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&herr,&lerr);
IA_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&herr,&lerr);
IA_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&herr,&lerr);
IA_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&herr,&lerr);
IA_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&herr,&lerr);
IA_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&herr,&lerr);
IA_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&herr,&lerr);
IA_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&herr,&lerr);
IA_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&herr,&lerr);
IA_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&herr,&lerr);
IA_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&herr,&lerr);
IA_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&herr,&lerr);
IA_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&herr,&lerr);
IA_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&herr,&lerr);
IA_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&herr,&lerr);
IA_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&herr,&lerr);
IA_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&herr,&lerr);
IA_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&herr,&lerr);
IA_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&herr,&lerr);
IA_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&herr,&lerr);
IA_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&herr,&lerr);
IA_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&herr,&lerr);
IA_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&herr,&lerr);
IA_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&herr,&lerr);
IA_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&herr,&lerr);
IA_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&herr,&lerr);
IA_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&herr,&lerr);
IA_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&herr,&lerr);
IA_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&herr,&lerr);
IA_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&herr,&lerr);
IA_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&herr,&lerr);
IA_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&herr,&lerr);
IA_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&herr,&lerr);
IA_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&herr,&lerr);
IA_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&herr,&lerr);
IA_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&herr,&lerr);
IA_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&herr,&lerr);
IA_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&herr,&lerr);
IA_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&herr,&lerr);
IA_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&herr,&lerr);
IA_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&herr,&lerr);
IA_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&herr,&lerr);
IA_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&herr,&lerr);
IA_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&herr,&lerr);
IA_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&herr,&lerr);
IA_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&herr,&lerr);
IA_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&herr,&lerr);
IA_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&herr,&lerr);
IA_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&herr,&lerr);
IA_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&herr,&lerr);
IA_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&herr,&lerr);
IA_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&herr,&lerr);
IA_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&herr,&lerr);
IA_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&herr,&lerr);
IA_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&herr,&lerr);
IA_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&herr,&lerr);
IA_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&herr,&lerr);
IA_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&herr,&lerr);
IA_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&herr,&lerr);
IA_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&herr,&lerr);
IA_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&herr,&lerr);
IA_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&herr,&lerr);
IA_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&herr,&lerr);
IA_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&herr,&lerr);
IA_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&herr,&lerr);
IA_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&herr,&lerr);
IA_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&herr,&lerr);
IA_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&herr,&lerr);
IA_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&herr,&lerr);
IA_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&herr,&lerr);
IA_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&herr,&lerr);
IA_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&herr,&lerr);
IA_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&herr,&lerr);
IA_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&herr,&lerr);
IA_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&herr,&lerr);
IA_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&herr,&lerr);
IA_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&herr,&lerr);
IA_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&herr,&lerr);
IA_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&herr,&lerr);
IA_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&herr,&lerr);
IA_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&herr,&lerr);
IA_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&herr,&lerr);
IA_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&herr,&lerr);
IA_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&herr,&lerr);
IA_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&herr,&lerr);
IA_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&herr,&lerr);
IA_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&herr,&lerr);
IA_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&herr,&lerr);
IA_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&herr,&lerr);
IA_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&herr,&lerr);
IA_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&herr,&lerr);
IA_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&herr,&lerr);
IA_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&herr,&lerr);
IA_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&herr,&lerr);
IA_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&herr,&lerr);
IA_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&herr,&lerr);
IA_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&herr,&lerr);
IA_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&herr,&lerr);
IA_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&herr,&lerr);
IA_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&herr,&lerr);
IA_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&herr,&lerr);
IA_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&herr,&lerr);
IA_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&herr,&lerr);
IA_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&herr,&lerr);
IA_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&herr,&lerr);
IA_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&herr,&lerr);
IA_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&herr,&lerr);
IA_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&herr,&lerr);
IA_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&herr,&lerr);
IA_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&herr,&lerr);
IA_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&herr,&lerr);
IA_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&herr,&lerr);
IA_InSphere(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez,&res,&herr,&lerr);
after = clock();
IAtime = ((double)(after - before))/CLOCKS_PER_SEC;
IAtotaltime = IAtotaltime + IAtime;


printf("FPS: %e \n", FPStime);
printf("Manual: %e \n", Manualtime);
printf("IA: %e \n\n", IAtime);

fprintf(save,"%e\n",FPStime);
fprintf(save,"%e\n",Manualtime);
fprintf(save,"%e\n",IAtime);
}

printf("Average FPS time : %f\n",FPStotaltime/sample);
printf("Average Manual time : %f\n",Manualtotaltime/sample);
printf("Average IA time : %f\n",IAtotaltime/sample);


fprintf(report,"Average FPS time : %f\n",FPStotaltime/sample);
fprintf(report,"Average Manual time : %f\n",Manualtotaltime/sample);
fprintf(report,"Average IA time : %f\n",IAtotaltime/sample);


fclose(save);
fclose(report);
}