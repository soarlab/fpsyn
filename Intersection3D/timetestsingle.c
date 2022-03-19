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
void FPSyn_Intersection3D(REAL ax1,REAL ax2,REAL ax3,REAL ay1,REAL ay2,REAL ay3,REAL az1,REAL az2,REAL az3,REAL bx1,REAL bx2,REAL bx3,REAL by1,REAL by2,REAL by3,REAL bz1,REAL bz2,REAL bz3,REAL cx1,REAL cx2,REAL cx3,REAL cy1,REAL cy2,REAL cy3,REAL cz1,REAL cz2,REAL cz3,REAL* result, REAL* error_bound){
//Calculation steps for input: ;
//Calculation steps for input: ;
REAL t1=ay1-by1;
REAL t3=az1-cz1;
REAL t5=ay1-cy1;
REAL t7=az1-bz1;
REAL t8=t1*t3;
REAL t9=t5*t7;
REAL t11=t8-t9;
REAL t13=ax1-bx1;
REAL t15=ax1-cx1;
REAL t16=t13*t5;
REAL t17=t1*t15;
REAL t19=t16-t17;
REAL t20=t13*t3;
REAL t21=t15*t7;
REAL t23=t20-t21;
REAL t25=bx3-ax3;
REAL t27=cz3-az3;
REAL t28=t25*t27;
REAL t29=cx3-ax3;
REAL t30=bz3-az3;
REAL t31=t29*t30;
REAL t33=t28-t31;
REAL t35=ay2-by2;
REAL t37=az2-cz2;
REAL t39=ay2-cy2;
REAL t41=az2-bz2;
REAL t42=t35*t37;
REAL t43=t39*t41;
REAL t45=t42-t43;
REAL t47=bx2-ax2;
REAL t49=cz2-az2;
REAL t50=t47*t49;
REAL t51=cx2-ax2;
REAL t52=bz2-az2;
REAL t53=t51*t52;
REAL t55=t50-t53;
REAL t57=ay3-by3;
REAL t59=az3-cz3;
REAL t61=ay3-cy3;
REAL t63=az3-bz3;
REAL t64=t57*t59;
REAL t65=t61*t63;
REAL t67=t64-t65;
REAL t69=ax3-bx3;
REAL t71=ax3-cx3;
REAL t72=ax3*t67;
REAL t73=t61*t69;
REAL t74=t57*t71;
REAL t76=t73-t74;
REAL t77=az3*t76;
REAL t78=t72+t77;
REAL t80=ax2-bx2;
REAL t82=ax2-cx2;
REAL t83=ax2*t45;
REAL t84=t39*t80;
REAL t85=t35*t82;
REAL t87=t84-t85;
REAL t88=az2*t87;
REAL t89=t83+t88;
REAL t92=-ay2;
REAL t93=cy2+t92;
REAL t94=by2+t92;
REAL t95=t47*t93;
REAL t96=t51*t94;
REAL t98=t95-t96;
REAL t99=-ay3;
REAL t100=by3+t99;
REAL t101=cy3+t99;
REAL t102=t100*t27;
REAL t103=t101*t30;
REAL t105=t102-t103;
REAL t106=t101*t25;
REAL t107=t100*t29;
REAL t109=t106-t107;
REAL t110=t49*t94;
REAL t111=t52*t93;
REAL t113=t110-t111;
REAL t115=t53-t50;
REAL t117=t31-t28;
REAL t118=t59*t69;
REAL t119=t63*t71;
REAL t121=t118-t119;
REAL t122=t121*t99;
REAL t123=t122+t78;
REAL t124=t123*t45;
REAL t125=t37*t80;
REAL t126=t41*t82;
REAL t128=t125-t126;
REAL t129=t128*t92;
REAL t130=t129+t89;
REAL t131=t130*t67;
REAL t133=t124-t131;
REAL t134=t133*t23;
REAL t135=t33*t99;
REAL t136=t135+t78;
REAL t137=t136*t55;
REAL t138=t55*t92;
REAL t139=t138+t89;
REAL t140=t139*t33;
REAL t142=t137-t140;
REAL t143=t11*t142;
REAL t145=t33*t45;
REAL t146=t55*t67;
REAL t148=t145-t146;
REAL t149=ax1*t11;
REAL t150=az1*t19;
REAL t151=ay1*t23;
REAL t153=t149+t150;
REAL t154=t153-t151;
REAL t155=t148*t154;
REAL t157=t134-t143;
REAL t158=t157-t155;
REAL t159=t117*t98;
REAL t160=t109*t115;
REAL t162=t159-t160;
REAL t163=t11*t162;
REAL t164=-t163;
REAL t165=t105*t115;
REAL t166=t113*t117;
REAL t168=t165-t166;
REAL t169=t168*t19;
REAL t171=t105*t98;
REAL t172=t109*t113;
REAL t174=t171-t172;
REAL t175=bx1-ax1;
REAL t176=cz1-az1;
REAL t177=t175*t176;
REAL t178=bz1-az1;
REAL t179=cx1-ax1;
REAL t180=t178*t179;
REAL t182=t177-t180;
REAL t183=t174*t182;
REAL t185=t164-t169;
REAL t186=t185-t183;
REAL t187=t158/t186;
//Calculation steps for Error-bound: ;
t74=fabs(t74);
t21=fabs(t21);
t136=fabs(t136);
az3=fabs(az3);
t99=fabs(t99);
t72=fabs(t72);
t129=fabs(t129);
t150=fabs(t150);
t154=fabs(t154);
t50=fabs(t50);
t131=fabs(t131);
t137=fabs(t137);
t139=fabs(t139);
t31=fabs(t31);
t158=fabs(t158);
t8=fabs(t8);
t84=fabs(t84);
t123=fabs(t123);
t145=fabs(t145);
t92=fabs(t92);
ax1=fabs(ax1);
t119=fabs(t119);
t9=fabs(t9);
t83=fabs(t83);
t142=fabs(t142);
t155=fabs(t155);
t125=fabs(t125);
t140=fabs(t140);
t134=fabs(t134);
t64=fabs(t64);
t138=fabs(t138);
t122=fabs(t122);
ax2=fabs(ax2);
t17=fabs(t17);
az2=fabs(az2);
t126=fabs(t126);
t143=fabs(t143);
az1=fabs(az1);
t130=fabs(t130);
t135=fabs(t135);
t28=fabs(t28);
t20=fabs(t20);
t65=fabs(t65);
ay1=fabs(ay1);
t73=fabs(t73);
t151=fabs(t151);
t149=fabs(t149);
t124=fabs(t124);
t42=fabs(t42);
t88=fabs(t88);
t148=fabs(t148);
t118=fabs(t118);
ax3=fabs(ax3);
t16=fabs(t16);
t77=fabs(t77);
t133=fabs(t133);
t53=fabs(t53);
t85=fabs(t85);
t146=fabs(t146);
t43=fabs(t43);
t163=fabs(t163);
t110=fabs(t110);
t168=fabs(t168);
t96=fabs(t96);
t111=fabs(t111);
t186=fabs(t186);
t160=fabs(t160);
t95=fabs(t95);
t180=fabs(t180);
t171=fabs(t171);
t107=fabs(t107);
t174=fabs(t174);
t102=fabs(t102);
t106=fabs(t106);
t183=fabs(t183);
t159=fabs(t159);
t169=fabs(t169);
t177=fabs(t177);
t103=fabs(t103);
t162=fabs(t162);
t165=fabs(t165);
t166=fabs(t166);
t172=fabs(t172);
REAL temp0=t20+t21;
REAL temp1=t133*temp0;
REAL temp2=t8+t9;
REAL temp3=t142*temp2;
REAL temp4=temp1+temp3;
REAL temp5=7.7715611723761e-16*temp4;
REAL temp6=2.22044604925031e-16*t150;
REAL temp7=2.22044604925031e-16*t154;
REAL temp8=ax1*temp2;
REAL temp9=6.66133814775094e-16*temp8;
REAL temp10=ay1*temp0;
REAL temp11=6.66133814775094e-16*temp10;
REAL temp12=t16+t17;
REAL temp13=az1*temp12;
REAL temp14=4.44089209850063e-16*temp13;
REAL temp15=temp6+temp7;
REAL temp16=temp15+temp9;
REAL temp17=temp11+temp16;
REAL temp18=temp14+temp17;
REAL temp19=t28+t31;
REAL temp20=t42+t43;
REAL temp21=t50+t53;
REAL temp22=t64+t65;
REAL temp23=6.66133814775094e-16*temp20;
REAL temp24=6.66133814775094e-16*temp22;
REAL temp25=6.66133814775094e-16*t99;
REAL temp26=ax3*temp24;
REAL temp27=t73+t74;
REAL temp28=az3*temp27;
REAL temp29=6.66133814775094e-16*temp28;
REAL temp30=temp26+temp29;
REAL temp31=6.66133814775094e-16*t92;
REAL temp32=ax2*temp23;
REAL temp33=t84+t85;
REAL temp34=az2*temp33;
REAL temp35=6.66133814775094e-16*temp34;
REAL temp36=temp32+temp35;
REAL temp37=2.22044604925031e-16*t155;
REAL temp38=t148*temp18;
REAL temp39=2.22044604925031e-16*t145;
REAL temp40=2.22044604925031e-16*t146;
REAL temp41=temp19*temp20;
REAL temp42=8.88178419700126e-16*temp41;
REAL temp43=temp21*temp22;
REAL temp44=8.88178419700126e-16*temp43;
REAL temp45=temp39+temp40;
REAL temp46=temp42+temp45;
REAL temp47=temp44+temp46;
REAL temp48=temp18*temp47;
REAL temp49=t136*temp21;
REAL temp50=6.66133814775094e-16*temp49;
REAL temp51=t139*temp19;
REAL temp52=6.66133814775094e-16*temp51;
REAL temp53=1.11022302462516e-16*t136;
REAL temp54=temp19*temp25;
REAL temp55=temp30+temp53;
REAL temp56=temp54+temp55;
REAL temp57=temp21*temp56;
REAL temp58=1.0*temp57;
REAL temp59=1.11022302462516e-16*t139;
REAL temp60=temp21*temp31;
REAL temp61=temp36+temp59;
REAL temp62=temp60+temp61;
REAL temp63=temp19*temp62;
REAL temp64=1.0*temp63;
REAL temp65=temp50+temp52;
REAL temp66=temp58+temp65;
REAL temp67=temp64+temp66;
REAL temp68=temp2*temp67;
REAL temp69=1.0*temp68;
REAL temp70=t123*temp23;
REAL temp71=t130*temp24;
REAL temp72=1.11022302462516e-16*t123;
REAL temp73=t118+t119;
REAL temp74=temp25*temp73;
REAL temp75=temp30+temp72;
REAL temp76=temp74+temp75;
REAL temp77=temp20*temp76;
REAL temp78=1.0*temp77;
REAL temp79=1.11022302462516e-16*t130;
REAL temp80=t125+t126;
REAL temp81=temp31*temp80;
REAL temp82=temp36+temp79;
REAL temp83=temp81+temp82;
REAL temp84=temp22*temp83;
REAL temp85=1.0*temp84;
REAL temp86=temp70+temp71;
REAL temp87=temp78+temp86;
REAL temp88=temp85+temp87;
REAL temp89=temp0*temp88;
REAL temp90=1.0*temp89;
REAL temp91=temp37+temp38;
REAL temp92=temp48+temp91;
REAL temp93=temp69+temp92;
REAL temp94=temp90+temp93;
REAL temp95=t154*temp47;
REAL temp96=temp5+temp94;
REAL temp97=temp95+temp96;
REAL temp98=t168*temp12;
REAL temp99=t177+t180;
REAL temp100=t174*temp99;
REAL temp101=temp100+temp98;
REAL temp102=7.7715611723761e-16*temp101;
REAL temp103=t102+t103;
REAL temp104=t110+t111;
REAL temp105=temp103*temp21;
REAL temp106=temp104*temp19;
REAL temp107=temp105+temp106;
REAL temp108=temp107*temp12;
REAL temp109=1.11022302462516e-15*temp108;
REAL temp110=t95+t96;
REAL temp111=temp103*temp110;
REAL temp112=t106+t107;
REAL temp113=temp104*temp112;
REAL temp114=temp111+temp113;
REAL temp115=temp114*temp99;
REAL temp116=1.11022302462516e-15*temp115;
REAL temp117=temp109+temp116;
REAL temp118=7.7715611723761e-16*t162;
REAL temp119=2.22044604925032e-16*t159;
REAL temp120=2.22044604925032e-16*t160;
REAL temp121=temp112*temp21;
REAL temp122=8.88178419700126e-16*temp121;
REAL temp123=temp110*temp19;
REAL temp124=8.88178419700126e-16*temp123;
REAL temp125=temp118+temp119;
REAL temp126=temp120+temp125;
REAL temp127=temp122+temp126;
REAL temp128=temp124+temp127;
REAL temp129=temp128*temp2;
REAL temp130=temp102+temp117;
REAL temp131=temp129+temp130;
*result=t187;
t187=fabs(t187);
REAL A=(t158*temp131+t186*temp97)*1.0000000001;
REAL B=(t186-1.0000000001*temp131)*0.9999999999999998;
REAL E=(A/t186*1.0000000000000002/B*1.0000000000000002+t187*1.1102230246251565e-16*1.0000000000000002)*1.0000000000000002;
*error_bound=t186<=temp131?1/0:E;
}
//IA Function
void IA_Intersection3D(REAL ax1,REAL ax2,REAL ax3,REAL ay1,REAL ay2,REAL ay3,REAL az1,REAL az2,REAL az3,REAL bx1,REAL bx2,REAL bx3,REAL by1,REAL by2,REAL by3,REAL bz1,REAL bz2,REAL bz3,REAL cx1,REAL cx2,REAL cx3,REAL cy1,REAL cy2,REAL cy3,REAL cz1,REAL cz2,REAL cz3,REAL* result, REAL* hbound, REAL* lbound){
REAL t1=ay1-by1;
REAL t3=az1-cz1;
REAL t5=ay1-cy1;
REAL t7=az1-bz1;
REAL t8=t1*t3;
REAL t9=t5*t7;
REAL t11=t8-t9;
REAL t13=ax1-bx1;
REAL t15=ax1-cx1;
REAL t16=t13*t5;
REAL t17=t1*t15;
REAL t19=t16-t17;
REAL t20=t13*t3;
REAL t21=t15*t7;
REAL t23=t20-t21;
REAL t25=bx3-ax3;
REAL t27=cz3-az3;
REAL t28=t25*t27;
REAL t29=cx3-ax3;
REAL t30=bz3-az3;
REAL t31=t29*t30;
REAL t33=t28-t31;
REAL t35=ay2-by2;
REAL t37=az2-cz2;
REAL t39=ay2-cy2;
REAL t41=az2-bz2;
REAL t42=t35*t37;
REAL t43=t39*t41;
REAL t45=t42-t43;
REAL t47=bx2-ax2;
REAL t49=cz2-az2;
REAL t50=t47*t49;
REAL t51=cx2-ax2;
REAL t52=bz2-az2;
REAL t53=t51*t52;
REAL t55=t50-t53;
REAL t57=ay3-by3;
REAL t59=az3-cz3;
REAL t61=ay3-cy3;
REAL t63=az3-bz3;
REAL t64=t57*t59;
REAL t65=t61*t63;
REAL t67=t64-t65;
REAL t69=ax3-bx3;
REAL t71=ax3-cx3;
REAL t72=ax3*t67;
REAL t73=t61*t69;
REAL t74=t57*t71;
REAL t76=t73-t74;
REAL t77=az3*t76;
REAL t78=t72+t77;
REAL t80=ax2-bx2;
REAL t82=ax2-cx2;
REAL t83=ax2*t45;
REAL t84=t39*t80;
REAL t85=t35*t82;
REAL t87=t84-t85;
REAL t88=az2*t87;
REAL t89=t83+t88;
REAL t92=-ay2;
REAL t93=cy2+t92;
REAL t94=by2+t92;
REAL t95=t47*t93;
REAL t96=t51*t94;
REAL t98=t95-t96;
REAL t99=-ay3;
REAL t100=by3+t99;
REAL t101=cy3+t99;
REAL t102=t100*t27;
REAL t103=t101*t30;
REAL t105=t102-t103;
REAL t106=t101*t25;
REAL t107=t100*t29;
REAL t109=t106-t107;
REAL t110=t49*t94;
REAL t111=t52*t93;
REAL t113=t110-t111;
REAL t115=t53-t50;
REAL t117=t31-t28;
REAL t118=t59*t69;
REAL t119=t63*t71;
REAL t121=t118-t119;
REAL t122=t121*t99;
REAL t123=t122+t78;
REAL t124=t123*t45;
REAL t125=t37*t80;
REAL t126=t41*t82;
REAL t128=t125-t126;
REAL t129=t128*t92;
REAL t130=t129+t89;
REAL t131=t130*t67;
REAL t133=t124-t131;
REAL t134=t133*t23;
REAL t135=t33*t99;
REAL t136=t135+t78;
REAL t137=t136*t55;
REAL t138=t55*t92;
REAL t139=t138+t89;
REAL t140=t139*t33;
REAL t142=t137-t140;
REAL t143=t11*t142;
REAL t145=t33*t45;
REAL t146=t55*t67;
REAL t148=t145-t146;
REAL t149=ax1*t11;
REAL t150=az1*t19;
REAL t151=ay1*t23;
REAL t153=t149+t150;
REAL t154=t153-t151;
REAL t155=t148*t154;
REAL t157=t134-t143;
REAL t158=t157-t155;
REAL t159=t117*t98;
REAL t160=t109*t115;
REAL t162=t159-t160;
REAL t163=t11*t162;
REAL t164=-t163;
REAL t165=t105*t115;
REAL t166=t113*t117;
REAL t168=t165-t166;
REAL t169=t168*t19;
REAL t171=t105*t98;
REAL t172=t109*t113;
REAL t174=t171-t172;
REAL t175=bx1-ax1;
REAL t176=cz1-az1;
REAL t177=t175*t176;
REAL t178=bz1-az1;
REAL t179=cx1-ax1;
REAL t180=t178*t179;
REAL t182=t177-t180;
REAL t183=t174*t182;
REAL t185=t164-t169;
REAL t186=t185-t183;
REAL t187=t158/t186;
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
REAL t7_e_1=fabs(t7);
REAL t7_e=eps*t7_e_1;
REAL t7_h=t7+t7_e;
REAL t7_l_n=t7_e-t7;
REAL t8_e_1=fabs(t8);
REAL t8_e=eps*t8_e_1;
REAL  t8_h, t8_l_n;
if(t1_h<=0){
if(t3_h<=0){
REAL t1_h_n=-t1_h;
REAL t8_l_n_1=t1_h_n*t3_h;
REAL t8_h_1=t1_l_n*t3_l_n;
t8_h=t8_h_1+t8_e;
t8_l_n=t8_l_n_1+t8_e;
}
else if(t3_l_n>=0){
REAL t8_l_n_1=t1_l_n*t3_h;
REAL t8_h_1=t1_l_n*t3_l_n;
t8_h=t8_h_1+t8_e;
t8_l_n=t8_l_n_1+t8_e;
}
else{
REAL t8_l_n_1=t1_l_n*t3_h;
REAL t1_h_n=-t1_h;
REAL t8_h_1=t1_h_n*t3_l_n;
t8_h=t8_h_1+t8_e;
t8_l_n=t8_l_n_1+t8_e;
}
}
else if(t1_l_n>=0){
if(t3_h<=0){
REAL t8_l_n_1=t1_h*t3_l_n;
REAL t8_h_1=t1_l_n*t3_l_n;
t8_h=t8_h_1+t8_e;
t8_l_n=t8_l_n_1+t8_e;
}
else if(t3_l_n>=0){
REAL t8_l_n_1_1=t1_l_n*t3_h;
REAL t8_l_n_1_2=t3_l_n*t1_h;
REAL t8_l_n_1=max(t8_l_n_1_1,t8_l_n_1_2);
REAL t8_h_1_1=t1_l_n*t3_l_n;
REAL t8_h_1_2=t1_h*t3_h;
REAL t8_h_1=max(t8_h_1_1,t8_h_1_2);
t8_h=t8_h_1+t8_e;
t8_l_n=t8_l_n_1+t8_e;
}
else{
REAL t8_l_n_1=t1_l_n*t3_h;
REAL t8_h_1=t1_h*t3_h;
t8_h=t8_h_1+t8_e;
t8_l_n=t8_l_n_1+t8_e;
}
}
else{
if(t3_h<=0){
REAL t8_l_n_1=t3_l_n*t1_h;
REAL t3_h_n=-t3_h;
REAL t8_h_1=t3_h_n*t1_l_n;
t8_h=t8_h_1+t8_e;
t8_l_n=t8_l_n_1+t8_e;
}
else if(t3_l_n>=0){
REAL t8_l_n_1=t3_l_n*t1_h;
REAL t8_h_1=t1_h*t3_h;
t8_h=t8_h_1+t8_e;
t8_l_n=t8_l_n_1+t8_e;
}
else{
REAL t1_l=-t1_l_n;
REAL t8_l_n_1=t1_l*t3_l_n;
REAL t8_h_1=t1_h*t3_h;
t8_h=t8_h_1+t8_e;
t8_l_n=t8_l_n_1+t8_e;
}
}
REAL t9_e_1=fabs(t9);
REAL t9_e=eps*t9_e_1;
REAL  t9_h, t9_l_n;
if(t5_h<=0){
if(t7_h<=0){
REAL t5_h_n=-t5_h;
REAL t9_l_n_1=t5_h_n*t7_h;
REAL t9_h_1=t5_l_n*t7_l_n;
t9_h=t9_h_1+t9_e;
t9_l_n=t9_l_n_1+t9_e;
}
else if(t7_l_n>=0){
REAL t9_l_n_1=t5_l_n*t7_h;
REAL t9_h_1=t5_l_n*t7_l_n;
t9_h=t9_h_1+t9_e;
t9_l_n=t9_l_n_1+t9_e;
}
else{
REAL t9_l_n_1=t5_l_n*t7_h;
REAL t5_h_n=-t5_h;
REAL t9_h_1=t5_h_n*t7_l_n;
t9_h=t9_h_1+t9_e;
t9_l_n=t9_l_n_1+t9_e;
}
}
else if(t5_l_n>=0){
if(t7_h<=0){
REAL t9_l_n_1=t5_h*t7_l_n;
REAL t9_h_1=t5_l_n*t7_l_n;
t9_h=t9_h_1+t9_e;
t9_l_n=t9_l_n_1+t9_e;
}
else if(t7_l_n>=0){
REAL t9_l_n_1_1=t5_l_n*t7_h;
REAL t9_l_n_1_2=t7_l_n*t5_h;
REAL t9_l_n_1=max(t9_l_n_1_1,t9_l_n_1_2);
REAL t9_h_1_1=t5_l_n*t7_l_n;
REAL t9_h_1_2=t5_h*t7_h;
REAL t9_h_1=max(t9_h_1_1,t9_h_1_2);
t9_h=t9_h_1+t9_e;
t9_l_n=t9_l_n_1+t9_e;
}
else{
REAL t9_l_n_1=t5_l_n*t7_h;
REAL t9_h_1=t5_h*t7_h;
t9_h=t9_h_1+t9_e;
t9_l_n=t9_l_n_1+t9_e;
}
}
else{
if(t7_h<=0){
REAL t9_l_n_1=t7_l_n*t5_h;
REAL t7_h_n=-t7_h;
REAL t9_h_1=t7_h_n*t5_l_n;
t9_h=t9_h_1+t9_e;
t9_l_n=t9_l_n_1+t9_e;
}
else if(t7_l_n>=0){
REAL t9_l_n_1=t7_l_n*t5_h;
REAL t9_h_1=t5_h*t7_h;
t9_h=t9_h_1+t9_e;
t9_l_n=t9_l_n_1+t9_e;
}
else{
REAL t5_l=-t5_l_n;
REAL t9_l_n_1=t5_l*t7_l_n;
REAL t9_h_1=t5_h*t7_h;
t9_h=t9_h_1+t9_e;
t9_l_n=t9_l_n_1+t9_e;
}
}
REAL t11_e_1=fabs(t11);
REAL t11_e=eps*t11_e_1;
REAL t11_h_1=t8_h+t9_l_n;
REAL t11_h=t11_h_1+t11_e;
REAL t11_l_n_1=t9_h+t8_l_n;
REAL t11_l_n=t11_l_n_1+t11_e;
REAL t13_e_1=fabs(t13);
REAL t13_e=eps*t13_e_1;
REAL t13_h=t13+t13_e;
REAL t13_l_n=t13_e-t13;
REAL t15_e_1=fabs(t15);
REAL t15_e=eps*t15_e_1;
REAL t15_h=t15+t15_e;
REAL t15_l_n=t15_e-t15;
REAL t16_e_1=fabs(t16);
REAL t16_e=eps*t16_e_1;
REAL  t16_h, t16_l_n;
if(t13_h<=0){
if(t5_h<=0){
REAL t13_h_n=-t13_h;
REAL t16_l_n_1=t13_h_n*t5_h;
REAL t16_h_1=t13_l_n*t5_l_n;
t16_h=t16_h_1+t16_e;
t16_l_n=t16_l_n_1+t16_e;
}
else if(t5_l_n>=0){
REAL t16_l_n_1=t13_l_n*t5_h;
REAL t16_h_1=t13_l_n*t5_l_n;
t16_h=t16_h_1+t16_e;
t16_l_n=t16_l_n_1+t16_e;
}
else{
REAL t16_l_n_1=t13_l_n*t5_h;
REAL t13_h_n=-t13_h;
REAL t16_h_1=t13_h_n*t5_l_n;
t16_h=t16_h_1+t16_e;
t16_l_n=t16_l_n_1+t16_e;
}
}
else if(t13_l_n>=0){
if(t5_h<=0){
REAL t16_l_n_1=t13_h*t5_l_n;
REAL t16_h_1=t13_l_n*t5_l_n;
t16_h=t16_h_1+t16_e;
t16_l_n=t16_l_n_1+t16_e;
}
else if(t5_l_n>=0){
REAL t16_l_n_1_1=t13_l_n*t5_h;
REAL t16_l_n_1_2=t5_l_n*t13_h;
REAL t16_l_n_1=max(t16_l_n_1_1,t16_l_n_1_2);
REAL t16_h_1_1=t13_l_n*t5_l_n;
REAL t16_h_1_2=t13_h*t5_h;
REAL t16_h_1=max(t16_h_1_1,t16_h_1_2);
t16_h=t16_h_1+t16_e;
t16_l_n=t16_l_n_1+t16_e;
}
else{
REAL t16_l_n_1=t13_l_n*t5_h;
REAL t16_h_1=t13_h*t5_h;
t16_h=t16_h_1+t16_e;
t16_l_n=t16_l_n_1+t16_e;
}
}
else{
if(t5_h<=0){
REAL t16_l_n_1=t5_l_n*t13_h;
REAL t5_h_n=-t5_h;
REAL t16_h_1=t5_h_n*t13_l_n;
t16_h=t16_h_1+t16_e;
t16_l_n=t16_l_n_1+t16_e;
}
else if(t5_l_n>=0){
REAL t16_l_n_1=t5_l_n*t13_h;
REAL t16_h_1=t13_h*t5_h;
t16_h=t16_h_1+t16_e;
t16_l_n=t16_l_n_1+t16_e;
}
else{
REAL t13_l=-t13_l_n;
REAL t16_l_n_1=t13_l*t5_l_n;
REAL t16_h_1=t13_h*t5_h;
t16_h=t16_h_1+t16_e;
t16_l_n=t16_l_n_1+t16_e;
}
}
REAL t17_e_1=fabs(t17);
REAL t17_e=eps*t17_e_1;
REAL  t17_h, t17_l_n;
if(t1_h<=0){
if(t15_h<=0){
REAL t1_h_n=-t1_h;
REAL t17_l_n_1=t1_h_n*t15_h;
REAL t17_h_1=t1_l_n*t15_l_n;
t17_h=t17_h_1+t17_e;
t17_l_n=t17_l_n_1+t17_e;
}
else if(t15_l_n>=0){
REAL t17_l_n_1=t1_l_n*t15_h;
REAL t17_h_1=t1_l_n*t15_l_n;
t17_h=t17_h_1+t17_e;
t17_l_n=t17_l_n_1+t17_e;
}
else{
REAL t17_l_n_1=t1_l_n*t15_h;
REAL t1_h_n=-t1_h;
REAL t17_h_1=t1_h_n*t15_l_n;
t17_h=t17_h_1+t17_e;
t17_l_n=t17_l_n_1+t17_e;
}
}
else if(t1_l_n>=0){
if(t15_h<=0){
REAL t17_l_n_1=t1_h*t15_l_n;
REAL t17_h_1=t1_l_n*t15_l_n;
t17_h=t17_h_1+t17_e;
t17_l_n=t17_l_n_1+t17_e;
}
else if(t15_l_n>=0){
REAL t17_l_n_1_1=t1_l_n*t15_h;
REAL t17_l_n_1_2=t15_l_n*t1_h;
REAL t17_l_n_1=max(t17_l_n_1_1,t17_l_n_1_2);
REAL t17_h_1_1=t1_l_n*t15_l_n;
REAL t17_h_1_2=t1_h*t15_h;
REAL t17_h_1=max(t17_h_1_1,t17_h_1_2);
t17_h=t17_h_1+t17_e;
t17_l_n=t17_l_n_1+t17_e;
}
else{
REAL t17_l_n_1=t1_l_n*t15_h;
REAL t17_h_1=t1_h*t15_h;
t17_h=t17_h_1+t17_e;
t17_l_n=t17_l_n_1+t17_e;
}
}
else{
if(t15_h<=0){
REAL t17_l_n_1=t15_l_n*t1_h;
REAL t15_h_n=-t15_h;
REAL t17_h_1=t15_h_n*t1_l_n;
t17_h=t17_h_1+t17_e;
t17_l_n=t17_l_n_1+t17_e;
}
else if(t15_l_n>=0){
REAL t17_l_n_1=t15_l_n*t1_h;
REAL t17_h_1=t1_h*t15_h;
t17_h=t17_h_1+t17_e;
t17_l_n=t17_l_n_1+t17_e;
}
else{
REAL t1_l=-t1_l_n;
REAL t17_l_n_1=t1_l*t15_l_n;
REAL t17_h_1=t1_h*t15_h;
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
if(t13_h<=0){
if(t3_h<=0){
REAL t13_h_n=-t13_h;
REAL t20_l_n_1=t13_h_n*t3_h;
REAL t20_h_1=t13_l_n*t3_l_n;
t20_h=t20_h_1+t20_e;
t20_l_n=t20_l_n_1+t20_e;
}
else if(t3_l_n>=0){
REAL t20_l_n_1=t13_l_n*t3_h;
REAL t20_h_1=t13_l_n*t3_l_n;
t20_h=t20_h_1+t20_e;
t20_l_n=t20_l_n_1+t20_e;
}
else{
REAL t20_l_n_1=t13_l_n*t3_h;
REAL t13_h_n=-t13_h;
REAL t20_h_1=t13_h_n*t3_l_n;
t20_h=t20_h_1+t20_e;
t20_l_n=t20_l_n_1+t20_e;
}
}
else if(t13_l_n>=0){
if(t3_h<=0){
REAL t20_l_n_1=t13_h*t3_l_n;
REAL t20_h_1=t13_l_n*t3_l_n;
t20_h=t20_h_1+t20_e;
t20_l_n=t20_l_n_1+t20_e;
}
else if(t3_l_n>=0){
REAL t20_l_n_1_1=t13_l_n*t3_h;
REAL t20_l_n_1_2=t3_l_n*t13_h;
REAL t20_l_n_1=max(t20_l_n_1_1,t20_l_n_1_2);
REAL t20_h_1_1=t13_l_n*t3_l_n;
REAL t20_h_1_2=t13_h*t3_h;
REAL t20_h_1=max(t20_h_1_1,t20_h_1_2);
t20_h=t20_h_1+t20_e;
t20_l_n=t20_l_n_1+t20_e;
}
else{
REAL t20_l_n_1=t13_l_n*t3_h;
REAL t20_h_1=t13_h*t3_h;
t20_h=t20_h_1+t20_e;
t20_l_n=t20_l_n_1+t20_e;
}
}
else{
if(t3_h<=0){
REAL t20_l_n_1=t3_l_n*t13_h;
REAL t3_h_n=-t3_h;
REAL t20_h_1=t3_h_n*t13_l_n;
t20_h=t20_h_1+t20_e;
t20_l_n=t20_l_n_1+t20_e;
}
else if(t3_l_n>=0){
REAL t20_l_n_1=t3_l_n*t13_h;
REAL t20_h_1=t13_h*t3_h;
t20_h=t20_h_1+t20_e;
t20_l_n=t20_l_n_1+t20_e;
}
else{
REAL t13_l=-t13_l_n;
REAL t20_l_n_1=t13_l*t3_l_n;
REAL t20_h_1=t13_h*t3_h;
t20_h=t20_h_1+t20_e;
t20_l_n=t20_l_n_1+t20_e;
}
}
REAL t21_e_1=fabs(t21);
REAL t21_e=eps*t21_e_1;
REAL  t21_h, t21_l_n;
if(t15_h<=0){
if(t7_h<=0){
REAL t15_h_n=-t15_h;
REAL t21_l_n_1=t15_h_n*t7_h;
REAL t21_h_1=t15_l_n*t7_l_n;
t21_h=t21_h_1+t21_e;
t21_l_n=t21_l_n_1+t21_e;
}
else if(t7_l_n>=0){
REAL t21_l_n_1=t15_l_n*t7_h;
REAL t21_h_1=t15_l_n*t7_l_n;
t21_h=t21_h_1+t21_e;
t21_l_n=t21_l_n_1+t21_e;
}
else{
REAL t21_l_n_1=t15_l_n*t7_h;
REAL t15_h_n=-t15_h;
REAL t21_h_1=t15_h_n*t7_l_n;
t21_h=t21_h_1+t21_e;
t21_l_n=t21_l_n_1+t21_e;
}
}
else if(t15_l_n>=0){
if(t7_h<=0){
REAL t21_l_n_1=t15_h*t7_l_n;
REAL t21_h_1=t15_l_n*t7_l_n;
t21_h=t21_h_1+t21_e;
t21_l_n=t21_l_n_1+t21_e;
}
else if(t7_l_n>=0){
REAL t21_l_n_1_1=t15_l_n*t7_h;
REAL t21_l_n_1_2=t7_l_n*t15_h;
REAL t21_l_n_1=max(t21_l_n_1_1,t21_l_n_1_2);
REAL t21_h_1_1=t15_l_n*t7_l_n;
REAL t21_h_1_2=t15_h*t7_h;
REAL t21_h_1=max(t21_h_1_1,t21_h_1_2);
t21_h=t21_h_1+t21_e;
t21_l_n=t21_l_n_1+t21_e;
}
else{
REAL t21_l_n_1=t15_l_n*t7_h;
REAL t21_h_1=t15_h*t7_h;
t21_h=t21_h_1+t21_e;
t21_l_n=t21_l_n_1+t21_e;
}
}
else{
if(t7_h<=0){
REAL t21_l_n_1=t7_l_n*t15_h;
REAL t7_h_n=-t7_h;
REAL t21_h_1=t7_h_n*t15_l_n;
t21_h=t21_h_1+t21_e;
t21_l_n=t21_l_n_1+t21_e;
}
else if(t7_l_n>=0){
REAL t21_l_n_1=t7_l_n*t15_h;
REAL t21_h_1=t15_h*t7_h;
t21_h=t21_h_1+t21_e;
t21_l_n=t21_l_n_1+t21_e;
}
else{
REAL t15_l=-t15_l_n;
REAL t21_l_n_1=t15_l*t7_l_n;
REAL t21_h_1=t15_h*t7_h;
t21_h=t21_h_1+t21_e;
t21_l_n=t21_l_n_1+t21_e;
}
}
REAL t23_e_1=fabs(t23);
REAL t23_e=eps*t23_e_1;
REAL t23_h_1=t20_h+t21_l_n;
REAL t23_h=t23_h_1+t23_e;
REAL t23_l_n_1=t21_h+t20_l_n;
REAL t23_l_n=t23_l_n_1+t23_e;
REAL t25_e_1=fabs(t25);
REAL t25_e=eps*t25_e_1;
REAL t25_h=t25+t25_e;
REAL t25_l_n=t25_e-t25;
REAL t27_e_1=fabs(t27);
REAL t27_e=eps*t27_e_1;
REAL t27_h=t27+t27_e;
REAL t27_l_n=t27_e-t27;
REAL t28_e_1=fabs(t28);
REAL t28_e=eps*t28_e_1;
REAL  t28_h, t28_l_n;
if(t25_h<=0){
if(t27_h<=0){
REAL t25_h_n=-t25_h;
REAL t28_l_n_1=t25_h_n*t27_h;
REAL t28_h_1=t25_l_n*t27_l_n;
t28_h=t28_h_1+t28_e;
t28_l_n=t28_l_n_1+t28_e;
}
else if(t27_l_n>=0){
REAL t28_l_n_1=t25_l_n*t27_h;
REAL t28_h_1=t25_l_n*t27_l_n;
t28_h=t28_h_1+t28_e;
t28_l_n=t28_l_n_1+t28_e;
}
else{
REAL t28_l_n_1=t25_l_n*t27_h;
REAL t25_h_n=-t25_h;
REAL t28_h_1=t25_h_n*t27_l_n;
t28_h=t28_h_1+t28_e;
t28_l_n=t28_l_n_1+t28_e;
}
}
else if(t25_l_n>=0){
if(t27_h<=0){
REAL t28_l_n_1=t25_h*t27_l_n;
REAL t28_h_1=t25_l_n*t27_l_n;
t28_h=t28_h_1+t28_e;
t28_l_n=t28_l_n_1+t28_e;
}
else if(t27_l_n>=0){
REAL t28_l_n_1_1=t25_l_n*t27_h;
REAL t28_l_n_1_2=t27_l_n*t25_h;
REAL t28_l_n_1=max(t28_l_n_1_1,t28_l_n_1_2);
REAL t28_h_1_1=t25_l_n*t27_l_n;
REAL t28_h_1_2=t25_h*t27_h;
REAL t28_h_1=max(t28_h_1_1,t28_h_1_2);
t28_h=t28_h_1+t28_e;
t28_l_n=t28_l_n_1+t28_e;
}
else{
REAL t28_l_n_1=t25_l_n*t27_h;
REAL t28_h_1=t25_h*t27_h;
t28_h=t28_h_1+t28_e;
t28_l_n=t28_l_n_1+t28_e;
}
}
else{
if(t27_h<=0){
REAL t28_l_n_1=t27_l_n*t25_h;
REAL t27_h_n=-t27_h;
REAL t28_h_1=t27_h_n*t25_l_n;
t28_h=t28_h_1+t28_e;
t28_l_n=t28_l_n_1+t28_e;
}
else if(t27_l_n>=0){
REAL t28_l_n_1=t27_l_n*t25_h;
REAL t28_h_1=t25_h*t27_h;
t28_h=t28_h_1+t28_e;
t28_l_n=t28_l_n_1+t28_e;
}
else{
REAL t25_l=-t25_l_n;
REAL t28_l_n_1=t25_l*t27_l_n;
REAL t28_h_1=t25_h*t27_h;
t28_h=t28_h_1+t28_e;
t28_l_n=t28_l_n_1+t28_e;
}
}
REAL t29_e_1=fabs(t29);
REAL t29_e=eps*t29_e_1;
REAL t29_h=t29+t29_e;
REAL t29_l_n=t29_e-t29;
REAL t30_e_1=fabs(t30);
REAL t30_e=eps*t30_e_1;
REAL t30_h=t30+t30_e;
REAL t30_l_n=t30_e-t30;
REAL t31_e_1=fabs(t31);
REAL t31_e=eps*t31_e_1;
REAL  t31_h, t31_l_n;
if(t29_h<=0){
if(t30_h<=0){
REAL t29_h_n=-t29_h;
REAL t31_l_n_1=t29_h_n*t30_h;
REAL t31_h_1=t29_l_n*t30_l_n;
t31_h=t31_h_1+t31_e;
t31_l_n=t31_l_n_1+t31_e;
}
else if(t30_l_n>=0){
REAL t31_l_n_1=t29_l_n*t30_h;
REAL t31_h_1=t29_l_n*t30_l_n;
t31_h=t31_h_1+t31_e;
t31_l_n=t31_l_n_1+t31_e;
}
else{
REAL t31_l_n_1=t29_l_n*t30_h;
REAL t29_h_n=-t29_h;
REAL t31_h_1=t29_h_n*t30_l_n;
t31_h=t31_h_1+t31_e;
t31_l_n=t31_l_n_1+t31_e;
}
}
else if(t29_l_n>=0){
if(t30_h<=0){
REAL t31_l_n_1=t29_h*t30_l_n;
REAL t31_h_1=t29_l_n*t30_l_n;
t31_h=t31_h_1+t31_e;
t31_l_n=t31_l_n_1+t31_e;
}
else if(t30_l_n>=0){
REAL t31_l_n_1_1=t29_l_n*t30_h;
REAL t31_l_n_1_2=t30_l_n*t29_h;
REAL t31_l_n_1=max(t31_l_n_1_1,t31_l_n_1_2);
REAL t31_h_1_1=t29_l_n*t30_l_n;
REAL t31_h_1_2=t29_h*t30_h;
REAL t31_h_1=max(t31_h_1_1,t31_h_1_2);
t31_h=t31_h_1+t31_e;
t31_l_n=t31_l_n_1+t31_e;
}
else{
REAL t31_l_n_1=t29_l_n*t30_h;
REAL t31_h_1=t29_h*t30_h;
t31_h=t31_h_1+t31_e;
t31_l_n=t31_l_n_1+t31_e;
}
}
else{
if(t30_h<=0){
REAL t31_l_n_1=t30_l_n*t29_h;
REAL t30_h_n=-t30_h;
REAL t31_h_1=t30_h_n*t29_l_n;
t31_h=t31_h_1+t31_e;
t31_l_n=t31_l_n_1+t31_e;
}
else if(t30_l_n>=0){
REAL t31_l_n_1=t30_l_n*t29_h;
REAL t31_h_1=t29_h*t30_h;
t31_h=t31_h_1+t31_e;
t31_l_n=t31_l_n_1+t31_e;
}
else{
REAL t29_l=-t29_l_n;
REAL t31_l_n_1=t29_l*t30_l_n;
REAL t31_h_1=t29_h*t30_h;
t31_h=t31_h_1+t31_e;
t31_l_n=t31_l_n_1+t31_e;
}
}
REAL t33_e_1=fabs(t33);
REAL t33_e=eps*t33_e_1;
REAL t33_h_1=t28_h+t31_l_n;
REAL t33_h=t33_h_1+t33_e;
REAL t33_l_n_1=t31_h+t28_l_n;
REAL t33_l_n=t33_l_n_1+t33_e;
REAL t35_e_1=fabs(t35);
REAL t35_e=eps*t35_e_1;
REAL t35_h=t35+t35_e;
REAL t35_l_n=t35_e-t35;
REAL t37_e_1=fabs(t37);
REAL t37_e=eps*t37_e_1;
REAL t37_h=t37+t37_e;
REAL t37_l_n=t37_e-t37;
REAL t39_e_1=fabs(t39);
REAL t39_e=eps*t39_e_1;
REAL t39_h=t39+t39_e;
REAL t39_l_n=t39_e-t39;
REAL t41_e_1=fabs(t41);
REAL t41_e=eps*t41_e_1;
REAL t41_h=t41+t41_e;
REAL t41_l_n=t41_e-t41;
REAL t42_e_1=fabs(t42);
REAL t42_e=eps*t42_e_1;
REAL  t42_h, t42_l_n;
if(t35_h<=0){
if(t37_h<=0){
REAL t35_h_n=-t35_h;
REAL t42_l_n_1=t35_h_n*t37_h;
REAL t42_h_1=t35_l_n*t37_l_n;
t42_h=t42_h_1+t42_e;
t42_l_n=t42_l_n_1+t42_e;
}
else if(t37_l_n>=0){
REAL t42_l_n_1=t35_l_n*t37_h;
REAL t42_h_1=t35_l_n*t37_l_n;
t42_h=t42_h_1+t42_e;
t42_l_n=t42_l_n_1+t42_e;
}
else{
REAL t42_l_n_1=t35_l_n*t37_h;
REAL t35_h_n=-t35_h;
REAL t42_h_1=t35_h_n*t37_l_n;
t42_h=t42_h_1+t42_e;
t42_l_n=t42_l_n_1+t42_e;
}
}
else if(t35_l_n>=0){
if(t37_h<=0){
REAL t42_l_n_1=t35_h*t37_l_n;
REAL t42_h_1=t35_l_n*t37_l_n;
t42_h=t42_h_1+t42_e;
t42_l_n=t42_l_n_1+t42_e;
}
else if(t37_l_n>=0){
REAL t42_l_n_1_1=t35_l_n*t37_h;
REAL t42_l_n_1_2=t37_l_n*t35_h;
REAL t42_l_n_1=max(t42_l_n_1_1,t42_l_n_1_2);
REAL t42_h_1_1=t35_l_n*t37_l_n;
REAL t42_h_1_2=t35_h*t37_h;
REAL t42_h_1=max(t42_h_1_1,t42_h_1_2);
t42_h=t42_h_1+t42_e;
t42_l_n=t42_l_n_1+t42_e;
}
else{
REAL t42_l_n_1=t35_l_n*t37_h;
REAL t42_h_1=t35_h*t37_h;
t42_h=t42_h_1+t42_e;
t42_l_n=t42_l_n_1+t42_e;
}
}
else{
if(t37_h<=0){
REAL t42_l_n_1=t37_l_n*t35_h;
REAL t37_h_n=-t37_h;
REAL t42_h_1=t37_h_n*t35_l_n;
t42_h=t42_h_1+t42_e;
t42_l_n=t42_l_n_1+t42_e;
}
else if(t37_l_n>=0){
REAL t42_l_n_1=t37_l_n*t35_h;
REAL t42_h_1=t35_h*t37_h;
t42_h=t42_h_1+t42_e;
t42_l_n=t42_l_n_1+t42_e;
}
else{
REAL t35_l=-t35_l_n;
REAL t42_l_n_1=t35_l*t37_l_n;
REAL t42_h_1=t35_h*t37_h;
t42_h=t42_h_1+t42_e;
t42_l_n=t42_l_n_1+t42_e;
}
}
REAL t43_e_1=fabs(t43);
REAL t43_e=eps*t43_e_1;
REAL  t43_h, t43_l_n;
if(t39_h<=0){
if(t41_h<=0){
REAL t39_h_n=-t39_h;
REAL t43_l_n_1=t39_h_n*t41_h;
REAL t43_h_1=t39_l_n*t41_l_n;
t43_h=t43_h_1+t43_e;
t43_l_n=t43_l_n_1+t43_e;
}
else if(t41_l_n>=0){
REAL t43_l_n_1=t39_l_n*t41_h;
REAL t43_h_1=t39_l_n*t41_l_n;
t43_h=t43_h_1+t43_e;
t43_l_n=t43_l_n_1+t43_e;
}
else{
REAL t43_l_n_1=t39_l_n*t41_h;
REAL t39_h_n=-t39_h;
REAL t43_h_1=t39_h_n*t41_l_n;
t43_h=t43_h_1+t43_e;
t43_l_n=t43_l_n_1+t43_e;
}
}
else if(t39_l_n>=0){
if(t41_h<=0){
REAL t43_l_n_1=t39_h*t41_l_n;
REAL t43_h_1=t39_l_n*t41_l_n;
t43_h=t43_h_1+t43_e;
t43_l_n=t43_l_n_1+t43_e;
}
else if(t41_l_n>=0){
REAL t43_l_n_1_1=t39_l_n*t41_h;
REAL t43_l_n_1_2=t41_l_n*t39_h;
REAL t43_l_n_1=max(t43_l_n_1_1,t43_l_n_1_2);
REAL t43_h_1_1=t39_l_n*t41_l_n;
REAL t43_h_1_2=t39_h*t41_h;
REAL t43_h_1=max(t43_h_1_1,t43_h_1_2);
t43_h=t43_h_1+t43_e;
t43_l_n=t43_l_n_1+t43_e;
}
else{
REAL t43_l_n_1=t39_l_n*t41_h;
REAL t43_h_1=t39_h*t41_h;
t43_h=t43_h_1+t43_e;
t43_l_n=t43_l_n_1+t43_e;
}
}
else{
if(t41_h<=0){
REAL t43_l_n_1=t41_l_n*t39_h;
REAL t41_h_n=-t41_h;
REAL t43_h_1=t41_h_n*t39_l_n;
t43_h=t43_h_1+t43_e;
t43_l_n=t43_l_n_1+t43_e;
}
else if(t41_l_n>=0){
REAL t43_l_n_1=t41_l_n*t39_h;
REAL t43_h_1=t39_h*t41_h;
t43_h=t43_h_1+t43_e;
t43_l_n=t43_l_n_1+t43_e;
}
else{
REAL t39_l=-t39_l_n;
REAL t43_l_n_1=t39_l*t41_l_n;
REAL t43_h_1=t39_h*t41_h;
t43_h=t43_h_1+t43_e;
t43_l_n=t43_l_n_1+t43_e;
}
}
REAL t45_e_1=fabs(t45);
REAL t45_e=eps*t45_e_1;
REAL t45_h_1=t42_h+t43_l_n;
REAL t45_h=t45_h_1+t45_e;
REAL t45_l_n_1=t43_h+t42_l_n;
REAL t45_l_n=t45_l_n_1+t45_e;
REAL t47_e_1=fabs(t47);
REAL t47_e=eps*t47_e_1;
REAL t47_h=t47+t47_e;
REAL t47_l_n=t47_e-t47;
REAL t49_e_1=fabs(t49);
REAL t49_e=eps*t49_e_1;
REAL t49_h=t49+t49_e;
REAL t49_l_n=t49_e-t49;
REAL t50_e_1=fabs(t50);
REAL t50_e=eps*t50_e_1;
REAL  t50_h, t50_l_n;
if(t47_h<=0){
if(t49_h<=0){
REAL t47_h_n=-t47_h;
REAL t50_l_n_1=t47_h_n*t49_h;
REAL t50_h_1=t47_l_n*t49_l_n;
t50_h=t50_h_1+t50_e;
t50_l_n=t50_l_n_1+t50_e;
}
else if(t49_l_n>=0){
REAL t50_l_n_1=t47_l_n*t49_h;
REAL t50_h_1=t47_l_n*t49_l_n;
t50_h=t50_h_1+t50_e;
t50_l_n=t50_l_n_1+t50_e;
}
else{
REAL t50_l_n_1=t47_l_n*t49_h;
REAL t47_h_n=-t47_h;
REAL t50_h_1=t47_h_n*t49_l_n;
t50_h=t50_h_1+t50_e;
t50_l_n=t50_l_n_1+t50_e;
}
}
else if(t47_l_n>=0){
if(t49_h<=0){
REAL t50_l_n_1=t47_h*t49_l_n;
REAL t50_h_1=t47_l_n*t49_l_n;
t50_h=t50_h_1+t50_e;
t50_l_n=t50_l_n_1+t50_e;
}
else if(t49_l_n>=0){
REAL t50_l_n_1_1=t47_l_n*t49_h;
REAL t50_l_n_1_2=t49_l_n*t47_h;
REAL t50_l_n_1=max(t50_l_n_1_1,t50_l_n_1_2);
REAL t50_h_1_1=t47_l_n*t49_l_n;
REAL t50_h_1_2=t47_h*t49_h;
REAL t50_h_1=max(t50_h_1_1,t50_h_1_2);
t50_h=t50_h_1+t50_e;
t50_l_n=t50_l_n_1+t50_e;
}
else{
REAL t50_l_n_1=t47_l_n*t49_h;
REAL t50_h_1=t47_h*t49_h;
t50_h=t50_h_1+t50_e;
t50_l_n=t50_l_n_1+t50_e;
}
}
else{
if(t49_h<=0){
REAL t50_l_n_1=t49_l_n*t47_h;
REAL t49_h_n=-t49_h;
REAL t50_h_1=t49_h_n*t47_l_n;
t50_h=t50_h_1+t50_e;
t50_l_n=t50_l_n_1+t50_e;
}
else if(t49_l_n>=0){
REAL t50_l_n_1=t49_l_n*t47_h;
REAL t50_h_1=t47_h*t49_h;
t50_h=t50_h_1+t50_e;
t50_l_n=t50_l_n_1+t50_e;
}
else{
REAL t47_l=-t47_l_n;
REAL t50_l_n_1=t47_l*t49_l_n;
REAL t50_h_1=t47_h*t49_h;
t50_h=t50_h_1+t50_e;
t50_l_n=t50_l_n_1+t50_e;
}
}
REAL t51_e_1=fabs(t51);
REAL t51_e=eps*t51_e_1;
REAL t51_h=t51+t51_e;
REAL t51_l_n=t51_e-t51;
REAL t52_e_1=fabs(t52);
REAL t52_e=eps*t52_e_1;
REAL t52_h=t52+t52_e;
REAL t52_l_n=t52_e-t52;
REAL t53_e_1=fabs(t53);
REAL t53_e=eps*t53_e_1;
REAL  t53_h, t53_l_n;
if(t51_h<=0){
if(t52_h<=0){
REAL t51_h_n=-t51_h;
REAL t53_l_n_1=t51_h_n*t52_h;
REAL t53_h_1=t51_l_n*t52_l_n;
t53_h=t53_h_1+t53_e;
t53_l_n=t53_l_n_1+t53_e;
}
else if(t52_l_n>=0){
REAL t53_l_n_1=t51_l_n*t52_h;
REAL t53_h_1=t51_l_n*t52_l_n;
t53_h=t53_h_1+t53_e;
t53_l_n=t53_l_n_1+t53_e;
}
else{
REAL t53_l_n_1=t51_l_n*t52_h;
REAL t51_h_n=-t51_h;
REAL t53_h_1=t51_h_n*t52_l_n;
t53_h=t53_h_1+t53_e;
t53_l_n=t53_l_n_1+t53_e;
}
}
else if(t51_l_n>=0){
if(t52_h<=0){
REAL t53_l_n_1=t51_h*t52_l_n;
REAL t53_h_1=t51_l_n*t52_l_n;
t53_h=t53_h_1+t53_e;
t53_l_n=t53_l_n_1+t53_e;
}
else if(t52_l_n>=0){
REAL t53_l_n_1_1=t51_l_n*t52_h;
REAL t53_l_n_1_2=t52_l_n*t51_h;
REAL t53_l_n_1=max(t53_l_n_1_1,t53_l_n_1_2);
REAL t53_h_1_1=t51_l_n*t52_l_n;
REAL t53_h_1_2=t51_h*t52_h;
REAL t53_h_1=max(t53_h_1_1,t53_h_1_2);
t53_h=t53_h_1+t53_e;
t53_l_n=t53_l_n_1+t53_e;
}
else{
REAL t53_l_n_1=t51_l_n*t52_h;
REAL t53_h_1=t51_h*t52_h;
t53_h=t53_h_1+t53_e;
t53_l_n=t53_l_n_1+t53_e;
}
}
else{
if(t52_h<=0){
REAL t53_l_n_1=t52_l_n*t51_h;
REAL t52_h_n=-t52_h;
REAL t53_h_1=t52_h_n*t51_l_n;
t53_h=t53_h_1+t53_e;
t53_l_n=t53_l_n_1+t53_e;
}
else if(t52_l_n>=0){
REAL t53_l_n_1=t52_l_n*t51_h;
REAL t53_h_1=t51_h*t52_h;
t53_h=t53_h_1+t53_e;
t53_l_n=t53_l_n_1+t53_e;
}
else{
REAL t51_l=-t51_l_n;
REAL t53_l_n_1=t51_l*t52_l_n;
REAL t53_h_1=t51_h*t52_h;
t53_h=t53_h_1+t53_e;
t53_l_n=t53_l_n_1+t53_e;
}
}
REAL t55_e_1=fabs(t55);
REAL t55_e=eps*t55_e_1;
REAL t55_h_1=t50_h+t53_l_n;
REAL t55_h=t55_h_1+t55_e;
REAL t55_l_n_1=t53_h+t50_l_n;
REAL t55_l_n=t55_l_n_1+t55_e;
REAL t57_e_1=fabs(t57);
REAL t57_e=eps*t57_e_1;
REAL t57_h=t57+t57_e;
REAL t57_l_n=t57_e-t57;
REAL t59_e_1=fabs(t59);
REAL t59_e=eps*t59_e_1;
REAL t59_h=t59+t59_e;
REAL t59_l_n=t59_e-t59;
REAL t61_e_1=fabs(t61);
REAL t61_e=eps*t61_e_1;
REAL t61_h=t61+t61_e;
REAL t61_l_n=t61_e-t61;
REAL t63_e_1=fabs(t63);
REAL t63_e=eps*t63_e_1;
REAL t63_h=t63+t63_e;
REAL t63_l_n=t63_e-t63;
REAL t64_e_1=fabs(t64);
REAL t64_e=eps*t64_e_1;
REAL  t64_h, t64_l_n;
if(t57_h<=0){
if(t59_h<=0){
REAL t57_h_n=-t57_h;
REAL t64_l_n_1=t57_h_n*t59_h;
REAL t64_h_1=t57_l_n*t59_l_n;
t64_h=t64_h_1+t64_e;
t64_l_n=t64_l_n_1+t64_e;
}
else if(t59_l_n>=0){
REAL t64_l_n_1=t57_l_n*t59_h;
REAL t64_h_1=t57_l_n*t59_l_n;
t64_h=t64_h_1+t64_e;
t64_l_n=t64_l_n_1+t64_e;
}
else{
REAL t64_l_n_1=t57_l_n*t59_h;
REAL t57_h_n=-t57_h;
REAL t64_h_1=t57_h_n*t59_l_n;
t64_h=t64_h_1+t64_e;
t64_l_n=t64_l_n_1+t64_e;
}
}
else if(t57_l_n>=0){
if(t59_h<=0){
REAL t64_l_n_1=t57_h*t59_l_n;
REAL t64_h_1=t57_l_n*t59_l_n;
t64_h=t64_h_1+t64_e;
t64_l_n=t64_l_n_1+t64_e;
}
else if(t59_l_n>=0){
REAL t64_l_n_1_1=t57_l_n*t59_h;
REAL t64_l_n_1_2=t59_l_n*t57_h;
REAL t64_l_n_1=max(t64_l_n_1_1,t64_l_n_1_2);
REAL t64_h_1_1=t57_l_n*t59_l_n;
REAL t64_h_1_2=t57_h*t59_h;
REAL t64_h_1=max(t64_h_1_1,t64_h_1_2);
t64_h=t64_h_1+t64_e;
t64_l_n=t64_l_n_1+t64_e;
}
else{
REAL t64_l_n_1=t57_l_n*t59_h;
REAL t64_h_1=t57_h*t59_h;
t64_h=t64_h_1+t64_e;
t64_l_n=t64_l_n_1+t64_e;
}
}
else{
if(t59_h<=0){
REAL t64_l_n_1=t59_l_n*t57_h;
REAL t59_h_n=-t59_h;
REAL t64_h_1=t59_h_n*t57_l_n;
t64_h=t64_h_1+t64_e;
t64_l_n=t64_l_n_1+t64_e;
}
else if(t59_l_n>=0){
REAL t64_l_n_1=t59_l_n*t57_h;
REAL t64_h_1=t57_h*t59_h;
t64_h=t64_h_1+t64_e;
t64_l_n=t64_l_n_1+t64_e;
}
else{
REAL t57_l=-t57_l_n;
REAL t64_l_n_1=t57_l*t59_l_n;
REAL t64_h_1=t57_h*t59_h;
t64_h=t64_h_1+t64_e;
t64_l_n=t64_l_n_1+t64_e;
}
}
REAL t65_e_1=fabs(t65);
REAL t65_e=eps*t65_e_1;
REAL  t65_h, t65_l_n;
if(t61_h<=0){
if(t63_h<=0){
REAL t61_h_n=-t61_h;
REAL t65_l_n_1=t61_h_n*t63_h;
REAL t65_h_1=t61_l_n*t63_l_n;
t65_h=t65_h_1+t65_e;
t65_l_n=t65_l_n_1+t65_e;
}
else if(t63_l_n>=0){
REAL t65_l_n_1=t61_l_n*t63_h;
REAL t65_h_1=t61_l_n*t63_l_n;
t65_h=t65_h_1+t65_e;
t65_l_n=t65_l_n_1+t65_e;
}
else{
REAL t65_l_n_1=t61_l_n*t63_h;
REAL t61_h_n=-t61_h;
REAL t65_h_1=t61_h_n*t63_l_n;
t65_h=t65_h_1+t65_e;
t65_l_n=t65_l_n_1+t65_e;
}
}
else if(t61_l_n>=0){
if(t63_h<=0){
REAL t65_l_n_1=t61_h*t63_l_n;
REAL t65_h_1=t61_l_n*t63_l_n;
t65_h=t65_h_1+t65_e;
t65_l_n=t65_l_n_1+t65_e;
}
else if(t63_l_n>=0){
REAL t65_l_n_1_1=t61_l_n*t63_h;
REAL t65_l_n_1_2=t63_l_n*t61_h;
REAL t65_l_n_1=max(t65_l_n_1_1,t65_l_n_1_2);
REAL t65_h_1_1=t61_l_n*t63_l_n;
REAL t65_h_1_2=t61_h*t63_h;
REAL t65_h_1=max(t65_h_1_1,t65_h_1_2);
t65_h=t65_h_1+t65_e;
t65_l_n=t65_l_n_1+t65_e;
}
else{
REAL t65_l_n_1=t61_l_n*t63_h;
REAL t65_h_1=t61_h*t63_h;
t65_h=t65_h_1+t65_e;
t65_l_n=t65_l_n_1+t65_e;
}
}
else{
if(t63_h<=0){
REAL t65_l_n_1=t63_l_n*t61_h;
REAL t63_h_n=-t63_h;
REAL t65_h_1=t63_h_n*t61_l_n;
t65_h=t65_h_1+t65_e;
t65_l_n=t65_l_n_1+t65_e;
}
else if(t63_l_n>=0){
REAL t65_l_n_1=t63_l_n*t61_h;
REAL t65_h_1=t61_h*t63_h;
t65_h=t65_h_1+t65_e;
t65_l_n=t65_l_n_1+t65_e;
}
else{
REAL t61_l=-t61_l_n;
REAL t65_l_n_1=t61_l*t63_l_n;
REAL t65_h_1=t61_h*t63_h;
t65_h=t65_h_1+t65_e;
t65_l_n=t65_l_n_1+t65_e;
}
}
REAL t67_e_1=fabs(t67);
REAL t67_e=eps*t67_e_1;
REAL t67_h_1=t64_h+t65_l_n;
REAL t67_h=t67_h_1+t67_e;
REAL t67_l_n_1=t65_h+t64_l_n;
REAL t67_l_n=t67_l_n_1+t67_e;
REAL t69_e_1=fabs(t69);
REAL t69_e=eps*t69_e_1;
REAL t69_h=t69+t69_e;
REAL t69_l_n=t69_e-t69;
REAL t71_e_1=fabs(t71);
REAL t71_e=eps*t71_e_1;
REAL t71_h=t71+t71_e;
REAL t71_l_n=t71_e-t71;
REAL t72_e_1=fabs(t72);
REAL t72_e=eps*t72_e_1;
REAL ax3_n=-ax3;
REAL  t72_h, t72_l_n;
if(ax3<=0){
if(t67_h<=0){
REAL ax3_n=-ax3;
REAL t72_l_n_1=ax3_n*t67_h;
REAL t72_h_1=ax3_n*t67_l_n;
t72_h=t72_h_1+t72_e;
t72_l_n=t72_l_n_1+t72_e;
}
else if(t67_l_n>=0){
REAL t72_l_n_1=ax3_n*t67_h;
REAL t72_h_1=ax3_n*t67_l_n;
t72_h=t72_h_1+t72_e;
t72_l_n=t72_l_n_1+t72_e;
}
else{
REAL t72_l_n_1=ax3_n*t67_h;
REAL ax3_n=-ax3;
REAL t72_h_1=ax3_n*t67_l_n;
t72_h=t72_h_1+t72_e;
t72_l_n=t72_l_n_1+t72_e;
}
}
else if(ax3_n>=0){
if(t67_h<=0){
REAL t72_l_n_1=ax3*t67_l_n;
REAL t72_h_1=ax3_n*t67_l_n;
t72_h=t72_h_1+t72_e;
t72_l_n=t72_l_n_1+t72_e;
}
else if(t67_l_n>=0){
REAL t72_l_n_1_1=ax3_n*t67_h;
REAL t72_l_n_1_2=t67_l_n*ax3;
REAL t72_l_n_1=max(t72_l_n_1_1,t72_l_n_1_2);
REAL t72_h_1_1=ax3_n*t67_l_n;
REAL t72_h_1_2=ax3*t67_h;
REAL t72_h_1=max(t72_h_1_1,t72_h_1_2);
t72_h=t72_h_1+t72_e;
t72_l_n=t72_l_n_1+t72_e;
}
else{
REAL t72_l_n_1=ax3_n*t67_h;
REAL t72_h_1=ax3*t67_h;
t72_h=t72_h_1+t72_e;
t72_l_n=t72_l_n_1+t72_e;
}
}
else{
if(t67_h<=0){
REAL t72_l_n_1=t67_l_n*ax3;
REAL t67_h_n=-t67_h;
REAL t72_h_1=t67_h_n*ax3_n;
t72_h=t72_h_1+t72_e;
t72_l_n=t72_l_n_1+t72_e;
}
else if(t67_l_n>=0){
REAL t72_l_n_1=t67_l_n*ax3;
REAL t72_h_1=ax3*t67_h;
t72_h=t72_h_1+t72_e;
t72_l_n=t72_l_n_1+t72_e;
}
else{
REAL ax3=-ax3_n;
REAL t72_l_n_1=ax3*t67_l_n;
REAL t72_h_1=ax3*t67_h;
t72_h=t72_h_1+t72_e;
t72_l_n=t72_l_n_1+t72_e;
}
}
REAL t73_e_1=fabs(t73);
REAL t73_e=eps*t73_e_1;
REAL  t73_h, t73_l_n;
if(t61_h<=0){
if(t69_h<=0){
REAL t61_h_n=-t61_h;
REAL t73_l_n_1=t61_h_n*t69_h;
REAL t73_h_1=t61_l_n*t69_l_n;
t73_h=t73_h_1+t73_e;
t73_l_n=t73_l_n_1+t73_e;
}
else if(t69_l_n>=0){
REAL t73_l_n_1=t61_l_n*t69_h;
REAL t73_h_1=t61_l_n*t69_l_n;
t73_h=t73_h_1+t73_e;
t73_l_n=t73_l_n_1+t73_e;
}
else{
REAL t73_l_n_1=t61_l_n*t69_h;
REAL t61_h_n=-t61_h;
REAL t73_h_1=t61_h_n*t69_l_n;
t73_h=t73_h_1+t73_e;
t73_l_n=t73_l_n_1+t73_e;
}
}
else if(t61_l_n>=0){
if(t69_h<=0){
REAL t73_l_n_1=t61_h*t69_l_n;
REAL t73_h_1=t61_l_n*t69_l_n;
t73_h=t73_h_1+t73_e;
t73_l_n=t73_l_n_1+t73_e;
}
else if(t69_l_n>=0){
REAL t73_l_n_1_1=t61_l_n*t69_h;
REAL t73_l_n_1_2=t69_l_n*t61_h;
REAL t73_l_n_1=max(t73_l_n_1_1,t73_l_n_1_2);
REAL t73_h_1_1=t61_l_n*t69_l_n;
REAL t73_h_1_2=t61_h*t69_h;
REAL t73_h_1=max(t73_h_1_1,t73_h_1_2);
t73_h=t73_h_1+t73_e;
t73_l_n=t73_l_n_1+t73_e;
}
else{
REAL t73_l_n_1=t61_l_n*t69_h;
REAL t73_h_1=t61_h*t69_h;
t73_h=t73_h_1+t73_e;
t73_l_n=t73_l_n_1+t73_e;
}
}
else{
if(t69_h<=0){
REAL t73_l_n_1=t69_l_n*t61_h;
REAL t69_h_n=-t69_h;
REAL t73_h_1=t69_h_n*t61_l_n;
t73_h=t73_h_1+t73_e;
t73_l_n=t73_l_n_1+t73_e;
}
else if(t69_l_n>=0){
REAL t73_l_n_1=t69_l_n*t61_h;
REAL t73_h_1=t61_h*t69_h;
t73_h=t73_h_1+t73_e;
t73_l_n=t73_l_n_1+t73_e;
}
else{
REAL t61_l=-t61_l_n;
REAL t73_l_n_1=t61_l*t69_l_n;
REAL t73_h_1=t61_h*t69_h;
t73_h=t73_h_1+t73_e;
t73_l_n=t73_l_n_1+t73_e;
}
}
REAL t74_e_1=fabs(t74);
REAL t74_e=eps*t74_e_1;
REAL  t74_h, t74_l_n;
if(t57_h<=0){
if(t71_h<=0){
REAL t57_h_n=-t57_h;
REAL t74_l_n_1=t57_h_n*t71_h;
REAL t74_h_1=t57_l_n*t71_l_n;
t74_h=t74_h_1+t74_e;
t74_l_n=t74_l_n_1+t74_e;
}
else if(t71_l_n>=0){
REAL t74_l_n_1=t57_l_n*t71_h;
REAL t74_h_1=t57_l_n*t71_l_n;
t74_h=t74_h_1+t74_e;
t74_l_n=t74_l_n_1+t74_e;
}
else{
REAL t74_l_n_1=t57_l_n*t71_h;
REAL t57_h_n=-t57_h;
REAL t74_h_1=t57_h_n*t71_l_n;
t74_h=t74_h_1+t74_e;
t74_l_n=t74_l_n_1+t74_e;
}
}
else if(t57_l_n>=0){
if(t71_h<=0){
REAL t74_l_n_1=t57_h*t71_l_n;
REAL t74_h_1=t57_l_n*t71_l_n;
t74_h=t74_h_1+t74_e;
t74_l_n=t74_l_n_1+t74_e;
}
else if(t71_l_n>=0){
REAL t74_l_n_1_1=t57_l_n*t71_h;
REAL t74_l_n_1_2=t71_l_n*t57_h;
REAL t74_l_n_1=max(t74_l_n_1_1,t74_l_n_1_2);
REAL t74_h_1_1=t57_l_n*t71_l_n;
REAL t74_h_1_2=t57_h*t71_h;
REAL t74_h_1=max(t74_h_1_1,t74_h_1_2);
t74_h=t74_h_1+t74_e;
t74_l_n=t74_l_n_1+t74_e;
}
else{
REAL t74_l_n_1=t57_l_n*t71_h;
REAL t74_h_1=t57_h*t71_h;
t74_h=t74_h_1+t74_e;
t74_l_n=t74_l_n_1+t74_e;
}
}
else{
if(t71_h<=0){
REAL t74_l_n_1=t71_l_n*t57_h;
REAL t71_h_n=-t71_h;
REAL t74_h_1=t71_h_n*t57_l_n;
t74_h=t74_h_1+t74_e;
t74_l_n=t74_l_n_1+t74_e;
}
else if(t71_l_n>=0){
REAL t74_l_n_1=t71_l_n*t57_h;
REAL t74_h_1=t57_h*t71_h;
t74_h=t74_h_1+t74_e;
t74_l_n=t74_l_n_1+t74_e;
}
else{
REAL t57_l=-t57_l_n;
REAL t74_l_n_1=t57_l*t71_l_n;
REAL t74_h_1=t57_h*t71_h;
t74_h=t74_h_1+t74_e;
t74_l_n=t74_l_n_1+t74_e;
}
}
REAL t76_e_1=fabs(t76);
REAL t76_e=eps*t76_e_1;
REAL t76_h_1=t73_h+t74_l_n;
REAL t76_h=t76_h_1+t76_e;
REAL t76_l_n_1=t74_h+t73_l_n;
REAL t76_l_n=t76_l_n_1+t76_e;
REAL t77_e_1=fabs(t77);
REAL t77_e=eps*t77_e_1;
REAL az3_n=-az3;
REAL  t77_h, t77_l_n;
if(az3<=0){
if(t76_h<=0){
REAL az3_n=-az3;
REAL t77_l_n_1=az3_n*t76_h;
REAL t77_h_1=az3_n*t76_l_n;
t77_h=t77_h_1+t77_e;
t77_l_n=t77_l_n_1+t77_e;
}
else if(t76_l_n>=0){
REAL t77_l_n_1=az3_n*t76_h;
REAL t77_h_1=az3_n*t76_l_n;
t77_h=t77_h_1+t77_e;
t77_l_n=t77_l_n_1+t77_e;
}
else{
REAL t77_l_n_1=az3_n*t76_h;
REAL az3_n=-az3;
REAL t77_h_1=az3_n*t76_l_n;
t77_h=t77_h_1+t77_e;
t77_l_n=t77_l_n_1+t77_e;
}
}
else if(az3_n>=0){
if(t76_h<=0){
REAL t77_l_n_1=az3*t76_l_n;
REAL t77_h_1=az3_n*t76_l_n;
t77_h=t77_h_1+t77_e;
t77_l_n=t77_l_n_1+t77_e;
}
else if(t76_l_n>=0){
REAL t77_l_n_1_1=az3_n*t76_h;
REAL t77_l_n_1_2=t76_l_n*az3;
REAL t77_l_n_1=max(t77_l_n_1_1,t77_l_n_1_2);
REAL t77_h_1_1=az3_n*t76_l_n;
REAL t77_h_1_2=az3*t76_h;
REAL t77_h_1=max(t77_h_1_1,t77_h_1_2);
t77_h=t77_h_1+t77_e;
t77_l_n=t77_l_n_1+t77_e;
}
else{
REAL t77_l_n_1=az3_n*t76_h;
REAL t77_h_1=az3*t76_h;
t77_h=t77_h_1+t77_e;
t77_l_n=t77_l_n_1+t77_e;
}
}
else{
if(t76_h<=0){
REAL t77_l_n_1=t76_l_n*az3;
REAL t76_h_n=-t76_h;
REAL t77_h_1=t76_h_n*az3_n;
t77_h=t77_h_1+t77_e;
t77_l_n=t77_l_n_1+t77_e;
}
else if(t76_l_n>=0){
REAL t77_l_n_1=t76_l_n*az3;
REAL t77_h_1=az3*t76_h;
t77_h=t77_h_1+t77_e;
t77_l_n=t77_l_n_1+t77_e;
}
else{
REAL az3=-az3_n;
REAL t77_l_n_1=az3*t76_l_n;
REAL t77_h_1=az3*t76_h;
t77_h=t77_h_1+t77_e;
t77_l_n=t77_l_n_1+t77_e;
}
}
REAL t78_e_1=fabs(t78);
REAL t78_e=eps*t78_e_1;
REAL t78_h_1=t72_h+t77_h;
REAL t78_h=t78_h_1+t78_e;
REAL t78_l_n_1=t72_l_n+t77_l_n;
REAL t78_l_n=t78_l_n_1+t78_e;
REAL t80_e_1=fabs(t80);
REAL t80_e=eps*t80_e_1;
REAL t80_h=t80+t80_e;
REAL t80_l_n=t80_e-t80;
REAL t82_e_1=fabs(t82);
REAL t82_e=eps*t82_e_1;
REAL t82_h=t82+t82_e;
REAL t82_l_n=t82_e-t82;
REAL t83_e_1=fabs(t83);
REAL t83_e=eps*t83_e_1;
REAL ax2_n=-ax2;
REAL  t83_h, t83_l_n;
if(ax2<=0){
if(t45_h<=0){
REAL ax2_n=-ax2;
REAL t83_l_n_1=ax2_n*t45_h;
REAL t83_h_1=ax2_n*t45_l_n;
t83_h=t83_h_1+t83_e;
t83_l_n=t83_l_n_1+t83_e;
}
else if(t45_l_n>=0){
REAL t83_l_n_1=ax2_n*t45_h;
REAL t83_h_1=ax2_n*t45_l_n;
t83_h=t83_h_1+t83_e;
t83_l_n=t83_l_n_1+t83_e;
}
else{
REAL t83_l_n_1=ax2_n*t45_h;
REAL ax2_n=-ax2;
REAL t83_h_1=ax2_n*t45_l_n;
t83_h=t83_h_1+t83_e;
t83_l_n=t83_l_n_1+t83_e;
}
}
else if(ax2_n>=0){
if(t45_h<=0){
REAL t83_l_n_1=ax2*t45_l_n;
REAL t83_h_1=ax2_n*t45_l_n;
t83_h=t83_h_1+t83_e;
t83_l_n=t83_l_n_1+t83_e;
}
else if(t45_l_n>=0){
REAL t83_l_n_1_1=ax2_n*t45_h;
REAL t83_l_n_1_2=t45_l_n*ax2;
REAL t83_l_n_1=max(t83_l_n_1_1,t83_l_n_1_2);
REAL t83_h_1_1=ax2_n*t45_l_n;
REAL t83_h_1_2=ax2*t45_h;
REAL t83_h_1=max(t83_h_1_1,t83_h_1_2);
t83_h=t83_h_1+t83_e;
t83_l_n=t83_l_n_1+t83_e;
}
else{
REAL t83_l_n_1=ax2_n*t45_h;
REAL t83_h_1=ax2*t45_h;
t83_h=t83_h_1+t83_e;
t83_l_n=t83_l_n_1+t83_e;
}
}
else{
if(t45_h<=0){
REAL t83_l_n_1=t45_l_n*ax2;
REAL t45_h_n=-t45_h;
REAL t83_h_1=t45_h_n*ax2_n;
t83_h=t83_h_1+t83_e;
t83_l_n=t83_l_n_1+t83_e;
}
else if(t45_l_n>=0){
REAL t83_l_n_1=t45_l_n*ax2;
REAL t83_h_1=ax2*t45_h;
t83_h=t83_h_1+t83_e;
t83_l_n=t83_l_n_1+t83_e;
}
else{
REAL ax2=-ax2_n;
REAL t83_l_n_1=ax2*t45_l_n;
REAL t83_h_1=ax2*t45_h;
t83_h=t83_h_1+t83_e;
t83_l_n=t83_l_n_1+t83_e;
}
}
REAL t84_e_1=fabs(t84);
REAL t84_e=eps*t84_e_1;
REAL  t84_h, t84_l_n;
if(t39_h<=0){
if(t80_h<=0){
REAL t39_h_n=-t39_h;
REAL t84_l_n_1=t39_h_n*t80_h;
REAL t84_h_1=t39_l_n*t80_l_n;
t84_h=t84_h_1+t84_e;
t84_l_n=t84_l_n_1+t84_e;
}
else if(t80_l_n>=0){
REAL t84_l_n_1=t39_l_n*t80_h;
REAL t84_h_1=t39_l_n*t80_l_n;
t84_h=t84_h_1+t84_e;
t84_l_n=t84_l_n_1+t84_e;
}
else{
REAL t84_l_n_1=t39_l_n*t80_h;
REAL t39_h_n=-t39_h;
REAL t84_h_1=t39_h_n*t80_l_n;
t84_h=t84_h_1+t84_e;
t84_l_n=t84_l_n_1+t84_e;
}
}
else if(t39_l_n>=0){
if(t80_h<=0){
REAL t84_l_n_1=t39_h*t80_l_n;
REAL t84_h_1=t39_l_n*t80_l_n;
t84_h=t84_h_1+t84_e;
t84_l_n=t84_l_n_1+t84_e;
}
else if(t80_l_n>=0){
REAL t84_l_n_1_1=t39_l_n*t80_h;
REAL t84_l_n_1_2=t80_l_n*t39_h;
REAL t84_l_n_1=max(t84_l_n_1_1,t84_l_n_1_2);
REAL t84_h_1_1=t39_l_n*t80_l_n;
REAL t84_h_1_2=t39_h*t80_h;
REAL t84_h_1=max(t84_h_1_1,t84_h_1_2);
t84_h=t84_h_1+t84_e;
t84_l_n=t84_l_n_1+t84_e;
}
else{
REAL t84_l_n_1=t39_l_n*t80_h;
REAL t84_h_1=t39_h*t80_h;
t84_h=t84_h_1+t84_e;
t84_l_n=t84_l_n_1+t84_e;
}
}
else{
if(t80_h<=0){
REAL t84_l_n_1=t80_l_n*t39_h;
REAL t80_h_n=-t80_h;
REAL t84_h_1=t80_h_n*t39_l_n;
t84_h=t84_h_1+t84_e;
t84_l_n=t84_l_n_1+t84_e;
}
else if(t80_l_n>=0){
REAL t84_l_n_1=t80_l_n*t39_h;
REAL t84_h_1=t39_h*t80_h;
t84_h=t84_h_1+t84_e;
t84_l_n=t84_l_n_1+t84_e;
}
else{
REAL t39_l=-t39_l_n;
REAL t84_l_n_1=t39_l*t80_l_n;
REAL t84_h_1=t39_h*t80_h;
t84_h=t84_h_1+t84_e;
t84_l_n=t84_l_n_1+t84_e;
}
}
REAL t85_e_1=fabs(t85);
REAL t85_e=eps*t85_e_1;
REAL  t85_h, t85_l_n;
if(t35_h<=0){
if(t82_h<=0){
REAL t35_h_n=-t35_h;
REAL t85_l_n_1=t35_h_n*t82_h;
REAL t85_h_1=t35_l_n*t82_l_n;
t85_h=t85_h_1+t85_e;
t85_l_n=t85_l_n_1+t85_e;
}
else if(t82_l_n>=0){
REAL t85_l_n_1=t35_l_n*t82_h;
REAL t85_h_1=t35_l_n*t82_l_n;
t85_h=t85_h_1+t85_e;
t85_l_n=t85_l_n_1+t85_e;
}
else{
REAL t85_l_n_1=t35_l_n*t82_h;
REAL t35_h_n=-t35_h;
REAL t85_h_1=t35_h_n*t82_l_n;
t85_h=t85_h_1+t85_e;
t85_l_n=t85_l_n_1+t85_e;
}
}
else if(t35_l_n>=0){
if(t82_h<=0){
REAL t85_l_n_1=t35_h*t82_l_n;
REAL t85_h_1=t35_l_n*t82_l_n;
t85_h=t85_h_1+t85_e;
t85_l_n=t85_l_n_1+t85_e;
}
else if(t82_l_n>=0){
REAL t85_l_n_1_1=t35_l_n*t82_h;
REAL t85_l_n_1_2=t82_l_n*t35_h;
REAL t85_l_n_1=max(t85_l_n_1_1,t85_l_n_1_2);
REAL t85_h_1_1=t35_l_n*t82_l_n;
REAL t85_h_1_2=t35_h*t82_h;
REAL t85_h_1=max(t85_h_1_1,t85_h_1_2);
t85_h=t85_h_1+t85_e;
t85_l_n=t85_l_n_1+t85_e;
}
else{
REAL t85_l_n_1=t35_l_n*t82_h;
REAL t85_h_1=t35_h*t82_h;
t85_h=t85_h_1+t85_e;
t85_l_n=t85_l_n_1+t85_e;
}
}
else{
if(t82_h<=0){
REAL t85_l_n_1=t82_l_n*t35_h;
REAL t82_h_n=-t82_h;
REAL t85_h_1=t82_h_n*t35_l_n;
t85_h=t85_h_1+t85_e;
t85_l_n=t85_l_n_1+t85_e;
}
else if(t82_l_n>=0){
REAL t85_l_n_1=t82_l_n*t35_h;
REAL t85_h_1=t35_h*t82_h;
t85_h=t85_h_1+t85_e;
t85_l_n=t85_l_n_1+t85_e;
}
else{
REAL t35_l=-t35_l_n;
REAL t85_l_n_1=t35_l*t82_l_n;
REAL t85_h_1=t35_h*t82_h;
t85_h=t85_h_1+t85_e;
t85_l_n=t85_l_n_1+t85_e;
}
}
REAL t87_e_1=fabs(t87);
REAL t87_e=eps*t87_e_1;
REAL t87_h_1=t84_h+t85_l_n;
REAL t87_h=t87_h_1+t87_e;
REAL t87_l_n_1=t85_h+t84_l_n;
REAL t87_l_n=t87_l_n_1+t87_e;
REAL t88_e_1=fabs(t88);
REAL t88_e=eps*t88_e_1;
REAL az2_n=-az2;
REAL  t88_h, t88_l_n;
if(az2<=0){
if(t87_h<=0){
REAL az2_n=-az2;
REAL t88_l_n_1=az2_n*t87_h;
REAL t88_h_1=az2_n*t87_l_n;
t88_h=t88_h_1+t88_e;
t88_l_n=t88_l_n_1+t88_e;
}
else if(t87_l_n>=0){
REAL t88_l_n_1=az2_n*t87_h;
REAL t88_h_1=az2_n*t87_l_n;
t88_h=t88_h_1+t88_e;
t88_l_n=t88_l_n_1+t88_e;
}
else{
REAL t88_l_n_1=az2_n*t87_h;
REAL az2_n=-az2;
REAL t88_h_1=az2_n*t87_l_n;
t88_h=t88_h_1+t88_e;
t88_l_n=t88_l_n_1+t88_e;
}
}
else if(az2_n>=0){
if(t87_h<=0){
REAL t88_l_n_1=az2*t87_l_n;
REAL t88_h_1=az2_n*t87_l_n;
t88_h=t88_h_1+t88_e;
t88_l_n=t88_l_n_1+t88_e;
}
else if(t87_l_n>=0){
REAL t88_l_n_1_1=az2_n*t87_h;
REAL t88_l_n_1_2=t87_l_n*az2;
REAL t88_l_n_1=max(t88_l_n_1_1,t88_l_n_1_2);
REAL t88_h_1_1=az2_n*t87_l_n;
REAL t88_h_1_2=az2*t87_h;
REAL t88_h_1=max(t88_h_1_1,t88_h_1_2);
t88_h=t88_h_1+t88_e;
t88_l_n=t88_l_n_1+t88_e;
}
else{
REAL t88_l_n_1=az2_n*t87_h;
REAL t88_h_1=az2*t87_h;
t88_h=t88_h_1+t88_e;
t88_l_n=t88_l_n_1+t88_e;
}
}
else{
if(t87_h<=0){
REAL t88_l_n_1=t87_l_n*az2;
REAL t87_h_n=-t87_h;
REAL t88_h_1=t87_h_n*az2_n;
t88_h=t88_h_1+t88_e;
t88_l_n=t88_l_n_1+t88_e;
}
else if(t87_l_n>=0){
REAL t88_l_n_1=t87_l_n*az2;
REAL t88_h_1=az2*t87_h;
t88_h=t88_h_1+t88_e;
t88_l_n=t88_l_n_1+t88_e;
}
else{
REAL az2=-az2_n;
REAL t88_l_n_1=az2*t87_l_n;
REAL t88_h_1=az2*t87_h;
t88_h=t88_h_1+t88_e;
t88_l_n=t88_l_n_1+t88_e;
}
}
REAL t89_e_1=fabs(t89);
REAL t89_e=eps*t89_e_1;
REAL t89_h_1=t83_h+t88_h;
REAL t89_h=t89_h_1+t89_e;
REAL t89_l_n_1=t83_l_n+t88_l_n;
REAL t89_l_n=t89_l_n_1+t89_e;
REAL ay2_n=-ay2;
REAL t92_h=ay2_n;
REAL t92_l_n=ay2;
REAL t93_e_1=fabs(t93);
REAL t93_e=eps*t93_e_1;
REAL cy2_n=-cy2;
REAL t93_h_1=cy2+t92_h;
REAL t93_h=t93_h_1+t93_e;
REAL t93_l_n_1=cy2_n+t92_l_n;
REAL t93_l_n=t93_l_n_1+t93_e;
REAL t94_e_1=fabs(t94);
REAL t94_e=eps*t94_e_1;
REAL by2_n=-by2;
REAL t94_h_1=by2+t92_h;
REAL t94_h=t94_h_1+t94_e;
REAL t94_l_n_1=by2_n+t92_l_n;
REAL t94_l_n=t94_l_n_1+t94_e;
REAL t95_e_1=fabs(t95);
REAL t95_e=eps*t95_e_1;
REAL  t95_h, t95_l_n;
if(t47_h<=0){
if(t93_h<=0){
REAL t47_h_n=-t47_h;
REAL t95_l_n_1=t47_h_n*t93_h;
REAL t95_h_1=t47_l_n*t93_l_n;
t95_h=t95_h_1+t95_e;
t95_l_n=t95_l_n_1+t95_e;
}
else if(t93_l_n>=0){
REAL t95_l_n_1=t47_l_n*t93_h;
REAL t95_h_1=t47_l_n*t93_l_n;
t95_h=t95_h_1+t95_e;
t95_l_n=t95_l_n_1+t95_e;
}
else{
REAL t95_l_n_1=t47_l_n*t93_h;
REAL t47_h_n=-t47_h;
REAL t95_h_1=t47_h_n*t93_l_n;
t95_h=t95_h_1+t95_e;
t95_l_n=t95_l_n_1+t95_e;
}
}
else if(t47_l_n>=0){
if(t93_h<=0){
REAL t95_l_n_1=t47_h*t93_l_n;
REAL t95_h_1=t47_l_n*t93_l_n;
t95_h=t95_h_1+t95_e;
t95_l_n=t95_l_n_1+t95_e;
}
else if(t93_l_n>=0){
REAL t95_l_n_1_1=t47_l_n*t93_h;
REAL t95_l_n_1_2=t93_l_n*t47_h;
REAL t95_l_n_1=max(t95_l_n_1_1,t95_l_n_1_2);
REAL t95_h_1_1=t47_l_n*t93_l_n;
REAL t95_h_1_2=t47_h*t93_h;
REAL t95_h_1=max(t95_h_1_1,t95_h_1_2);
t95_h=t95_h_1+t95_e;
t95_l_n=t95_l_n_1+t95_e;
}
else{
REAL t95_l_n_1=t47_l_n*t93_h;
REAL t95_h_1=t47_h*t93_h;
t95_h=t95_h_1+t95_e;
t95_l_n=t95_l_n_1+t95_e;
}
}
else{
if(t93_h<=0){
REAL t95_l_n_1=t93_l_n*t47_h;
REAL t93_h_n=-t93_h;
REAL t95_h_1=t93_h_n*t47_l_n;
t95_h=t95_h_1+t95_e;
t95_l_n=t95_l_n_1+t95_e;
}
else if(t93_l_n>=0){
REAL t95_l_n_1=t93_l_n*t47_h;
REAL t95_h_1=t47_h*t93_h;
t95_h=t95_h_1+t95_e;
t95_l_n=t95_l_n_1+t95_e;
}
else{
REAL t47_l=-t47_l_n;
REAL t95_l_n_1=t47_l*t93_l_n;
REAL t95_h_1=t47_h*t93_h;
t95_h=t95_h_1+t95_e;
t95_l_n=t95_l_n_1+t95_e;
}
}
REAL t96_e_1=fabs(t96);
REAL t96_e=eps*t96_e_1;
REAL  t96_h, t96_l_n;
if(t51_h<=0){
if(t94_h<=0){
REAL t51_h_n=-t51_h;
REAL t96_l_n_1=t51_h_n*t94_h;
REAL t96_h_1=t51_l_n*t94_l_n;
t96_h=t96_h_1+t96_e;
t96_l_n=t96_l_n_1+t96_e;
}
else if(t94_l_n>=0){
REAL t96_l_n_1=t51_l_n*t94_h;
REAL t96_h_1=t51_l_n*t94_l_n;
t96_h=t96_h_1+t96_e;
t96_l_n=t96_l_n_1+t96_e;
}
else{
REAL t96_l_n_1=t51_l_n*t94_h;
REAL t51_h_n=-t51_h;
REAL t96_h_1=t51_h_n*t94_l_n;
t96_h=t96_h_1+t96_e;
t96_l_n=t96_l_n_1+t96_e;
}
}
else if(t51_l_n>=0){
if(t94_h<=0){
REAL t96_l_n_1=t51_h*t94_l_n;
REAL t96_h_1=t51_l_n*t94_l_n;
t96_h=t96_h_1+t96_e;
t96_l_n=t96_l_n_1+t96_e;
}
else if(t94_l_n>=0){
REAL t96_l_n_1_1=t51_l_n*t94_h;
REAL t96_l_n_1_2=t94_l_n*t51_h;
REAL t96_l_n_1=max(t96_l_n_1_1,t96_l_n_1_2);
REAL t96_h_1_1=t51_l_n*t94_l_n;
REAL t96_h_1_2=t51_h*t94_h;
REAL t96_h_1=max(t96_h_1_1,t96_h_1_2);
t96_h=t96_h_1+t96_e;
t96_l_n=t96_l_n_1+t96_e;
}
else{
REAL t96_l_n_1=t51_l_n*t94_h;
REAL t96_h_1=t51_h*t94_h;
t96_h=t96_h_1+t96_e;
t96_l_n=t96_l_n_1+t96_e;
}
}
else{
if(t94_h<=0){
REAL t96_l_n_1=t94_l_n*t51_h;
REAL t94_h_n=-t94_h;
REAL t96_h_1=t94_h_n*t51_l_n;
t96_h=t96_h_1+t96_e;
t96_l_n=t96_l_n_1+t96_e;
}
else if(t94_l_n>=0){
REAL t96_l_n_1=t94_l_n*t51_h;
REAL t96_h_1=t51_h*t94_h;
t96_h=t96_h_1+t96_e;
t96_l_n=t96_l_n_1+t96_e;
}
else{
REAL t51_l=-t51_l_n;
REAL t96_l_n_1=t51_l*t94_l_n;
REAL t96_h_1=t51_h*t94_h;
t96_h=t96_h_1+t96_e;
t96_l_n=t96_l_n_1+t96_e;
}
}
REAL t98_e_1=fabs(t98);
REAL t98_e=eps*t98_e_1;
REAL t98_h_1=t95_h+t96_l_n;
REAL t98_h=t98_h_1+t98_e;
REAL t98_l_n_1=t96_h+t95_l_n;
REAL t98_l_n=t98_l_n_1+t98_e;
REAL ay3_n=-ay3;
REAL t99_h=ay3_n;
REAL t99_l_n=ay3;
REAL t100_e_1=fabs(t100);
REAL t100_e=eps*t100_e_1;
REAL by3_n=-by3;
REAL t100_h_1=by3+t99_h;
REAL t100_h=t100_h_1+t100_e;
REAL t100_l_n_1=by3_n+t99_l_n;
REAL t100_l_n=t100_l_n_1+t100_e;
REAL t101_e_1=fabs(t101);
REAL t101_e=eps*t101_e_1;
REAL cy3_n=-cy3;
REAL t101_h_1=cy3+t99_h;
REAL t101_h=t101_h_1+t101_e;
REAL t101_l_n_1=cy3_n+t99_l_n;
REAL t101_l_n=t101_l_n_1+t101_e;
REAL t102_e_1=fabs(t102);
REAL t102_e=eps*t102_e_1;
REAL  t102_h, t102_l_n;
if(t100_h<=0){
if(t27_h<=0){
REAL t100_h_n=-t100_h;
REAL t102_l_n_1=t100_h_n*t27_h;
REAL t102_h_1=t100_l_n*t27_l_n;
t102_h=t102_h_1+t102_e;
t102_l_n=t102_l_n_1+t102_e;
}
else if(t27_l_n>=0){
REAL t102_l_n_1=t100_l_n*t27_h;
REAL t102_h_1=t100_l_n*t27_l_n;
t102_h=t102_h_1+t102_e;
t102_l_n=t102_l_n_1+t102_e;
}
else{
REAL t102_l_n_1=t100_l_n*t27_h;
REAL t100_h_n=-t100_h;
REAL t102_h_1=t100_h_n*t27_l_n;
t102_h=t102_h_1+t102_e;
t102_l_n=t102_l_n_1+t102_e;
}
}
else if(t100_l_n>=0){
if(t27_h<=0){
REAL t102_l_n_1=t100_h*t27_l_n;
REAL t102_h_1=t100_l_n*t27_l_n;
t102_h=t102_h_1+t102_e;
t102_l_n=t102_l_n_1+t102_e;
}
else if(t27_l_n>=0){
REAL t102_l_n_1_1=t100_l_n*t27_h;
REAL t102_l_n_1_2=t27_l_n*t100_h;
REAL t102_l_n_1=max(t102_l_n_1_1,t102_l_n_1_2);
REAL t102_h_1_1=t100_l_n*t27_l_n;
REAL t102_h_1_2=t100_h*t27_h;
REAL t102_h_1=max(t102_h_1_1,t102_h_1_2);
t102_h=t102_h_1+t102_e;
t102_l_n=t102_l_n_1+t102_e;
}
else{
REAL t102_l_n_1=t100_l_n*t27_h;
REAL t102_h_1=t100_h*t27_h;
t102_h=t102_h_1+t102_e;
t102_l_n=t102_l_n_1+t102_e;
}
}
else{
if(t27_h<=0){
REAL t102_l_n_1=t27_l_n*t100_h;
REAL t27_h_n=-t27_h;
REAL t102_h_1=t27_h_n*t100_l_n;
t102_h=t102_h_1+t102_e;
t102_l_n=t102_l_n_1+t102_e;
}
else if(t27_l_n>=0){
REAL t102_l_n_1=t27_l_n*t100_h;
REAL t102_h_1=t100_h*t27_h;
t102_h=t102_h_1+t102_e;
t102_l_n=t102_l_n_1+t102_e;
}
else{
REAL t100_l=-t100_l_n;
REAL t102_l_n_1=t100_l*t27_l_n;
REAL t102_h_1=t100_h*t27_h;
t102_h=t102_h_1+t102_e;
t102_l_n=t102_l_n_1+t102_e;
}
}
REAL t103_e_1=fabs(t103);
REAL t103_e=eps*t103_e_1;
REAL  t103_h, t103_l_n;
if(t101_h<=0){
if(t30_h<=0){
REAL t101_h_n=-t101_h;
REAL t103_l_n_1=t101_h_n*t30_h;
REAL t103_h_1=t101_l_n*t30_l_n;
t103_h=t103_h_1+t103_e;
t103_l_n=t103_l_n_1+t103_e;
}
else if(t30_l_n>=0){
REAL t103_l_n_1=t101_l_n*t30_h;
REAL t103_h_1=t101_l_n*t30_l_n;
t103_h=t103_h_1+t103_e;
t103_l_n=t103_l_n_1+t103_e;
}
else{
REAL t103_l_n_1=t101_l_n*t30_h;
REAL t101_h_n=-t101_h;
REAL t103_h_1=t101_h_n*t30_l_n;
t103_h=t103_h_1+t103_e;
t103_l_n=t103_l_n_1+t103_e;
}
}
else if(t101_l_n>=0){
if(t30_h<=0){
REAL t103_l_n_1=t101_h*t30_l_n;
REAL t103_h_1=t101_l_n*t30_l_n;
t103_h=t103_h_1+t103_e;
t103_l_n=t103_l_n_1+t103_e;
}
else if(t30_l_n>=0){
REAL t103_l_n_1_1=t101_l_n*t30_h;
REAL t103_l_n_1_2=t30_l_n*t101_h;
REAL t103_l_n_1=max(t103_l_n_1_1,t103_l_n_1_2);
REAL t103_h_1_1=t101_l_n*t30_l_n;
REAL t103_h_1_2=t101_h*t30_h;
REAL t103_h_1=max(t103_h_1_1,t103_h_1_2);
t103_h=t103_h_1+t103_e;
t103_l_n=t103_l_n_1+t103_e;
}
else{
REAL t103_l_n_1=t101_l_n*t30_h;
REAL t103_h_1=t101_h*t30_h;
t103_h=t103_h_1+t103_e;
t103_l_n=t103_l_n_1+t103_e;
}
}
else{
if(t30_h<=0){
REAL t103_l_n_1=t30_l_n*t101_h;
REAL t30_h_n=-t30_h;
REAL t103_h_1=t30_h_n*t101_l_n;
t103_h=t103_h_1+t103_e;
t103_l_n=t103_l_n_1+t103_e;
}
else if(t30_l_n>=0){
REAL t103_l_n_1=t30_l_n*t101_h;
REAL t103_h_1=t101_h*t30_h;
t103_h=t103_h_1+t103_e;
t103_l_n=t103_l_n_1+t103_e;
}
else{
REAL t101_l=-t101_l_n;
REAL t103_l_n_1=t101_l*t30_l_n;
REAL t103_h_1=t101_h*t30_h;
t103_h=t103_h_1+t103_e;
t103_l_n=t103_l_n_1+t103_e;
}
}
REAL t105_e_1=fabs(t105);
REAL t105_e=eps*t105_e_1;
REAL t105_h_1=t102_h+t103_l_n;
REAL t105_h=t105_h_1+t105_e;
REAL t105_l_n_1=t103_h+t102_l_n;
REAL t105_l_n=t105_l_n_1+t105_e;
REAL t106_e_1=fabs(t106);
REAL t106_e=eps*t106_e_1;
REAL  t106_h, t106_l_n;
if(t101_h<=0){
if(t25_h<=0){
REAL t101_h_n=-t101_h;
REAL t106_l_n_1=t101_h_n*t25_h;
REAL t106_h_1=t101_l_n*t25_l_n;
t106_h=t106_h_1+t106_e;
t106_l_n=t106_l_n_1+t106_e;
}
else if(t25_l_n>=0){
REAL t106_l_n_1=t101_l_n*t25_h;
REAL t106_h_1=t101_l_n*t25_l_n;
t106_h=t106_h_1+t106_e;
t106_l_n=t106_l_n_1+t106_e;
}
else{
REAL t106_l_n_1=t101_l_n*t25_h;
REAL t101_h_n=-t101_h;
REAL t106_h_1=t101_h_n*t25_l_n;
t106_h=t106_h_1+t106_e;
t106_l_n=t106_l_n_1+t106_e;
}
}
else if(t101_l_n>=0){
if(t25_h<=0){
REAL t106_l_n_1=t101_h*t25_l_n;
REAL t106_h_1=t101_l_n*t25_l_n;
t106_h=t106_h_1+t106_e;
t106_l_n=t106_l_n_1+t106_e;
}
else if(t25_l_n>=0){
REAL t106_l_n_1_1=t101_l_n*t25_h;
REAL t106_l_n_1_2=t25_l_n*t101_h;
REAL t106_l_n_1=max(t106_l_n_1_1,t106_l_n_1_2);
REAL t106_h_1_1=t101_l_n*t25_l_n;
REAL t106_h_1_2=t101_h*t25_h;
REAL t106_h_1=max(t106_h_1_1,t106_h_1_2);
t106_h=t106_h_1+t106_e;
t106_l_n=t106_l_n_1+t106_e;
}
else{
REAL t106_l_n_1=t101_l_n*t25_h;
REAL t106_h_1=t101_h*t25_h;
t106_h=t106_h_1+t106_e;
t106_l_n=t106_l_n_1+t106_e;
}
}
else{
if(t25_h<=0){
REAL t106_l_n_1=t25_l_n*t101_h;
REAL t25_h_n=-t25_h;
REAL t106_h_1=t25_h_n*t101_l_n;
t106_h=t106_h_1+t106_e;
t106_l_n=t106_l_n_1+t106_e;
}
else if(t25_l_n>=0){
REAL t106_l_n_1=t25_l_n*t101_h;
REAL t106_h_1=t101_h*t25_h;
t106_h=t106_h_1+t106_e;
t106_l_n=t106_l_n_1+t106_e;
}
else{
REAL t101_l=-t101_l_n;
REAL t106_l_n_1=t101_l*t25_l_n;
REAL t106_h_1=t101_h*t25_h;
t106_h=t106_h_1+t106_e;
t106_l_n=t106_l_n_1+t106_e;
}
}
REAL t107_e_1=fabs(t107);
REAL t107_e=eps*t107_e_1;
REAL  t107_h, t107_l_n;
if(t100_h<=0){
if(t29_h<=0){
REAL t100_h_n=-t100_h;
REAL t107_l_n_1=t100_h_n*t29_h;
REAL t107_h_1=t100_l_n*t29_l_n;
t107_h=t107_h_1+t107_e;
t107_l_n=t107_l_n_1+t107_e;
}
else if(t29_l_n>=0){
REAL t107_l_n_1=t100_l_n*t29_h;
REAL t107_h_1=t100_l_n*t29_l_n;
t107_h=t107_h_1+t107_e;
t107_l_n=t107_l_n_1+t107_e;
}
else{
REAL t107_l_n_1=t100_l_n*t29_h;
REAL t100_h_n=-t100_h;
REAL t107_h_1=t100_h_n*t29_l_n;
t107_h=t107_h_1+t107_e;
t107_l_n=t107_l_n_1+t107_e;
}
}
else if(t100_l_n>=0){
if(t29_h<=0){
REAL t107_l_n_1=t100_h*t29_l_n;
REAL t107_h_1=t100_l_n*t29_l_n;
t107_h=t107_h_1+t107_e;
t107_l_n=t107_l_n_1+t107_e;
}
else if(t29_l_n>=0){
REAL t107_l_n_1_1=t100_l_n*t29_h;
REAL t107_l_n_1_2=t29_l_n*t100_h;
REAL t107_l_n_1=max(t107_l_n_1_1,t107_l_n_1_2);
REAL t107_h_1_1=t100_l_n*t29_l_n;
REAL t107_h_1_2=t100_h*t29_h;
REAL t107_h_1=max(t107_h_1_1,t107_h_1_2);
t107_h=t107_h_1+t107_e;
t107_l_n=t107_l_n_1+t107_e;
}
else{
REAL t107_l_n_1=t100_l_n*t29_h;
REAL t107_h_1=t100_h*t29_h;
t107_h=t107_h_1+t107_e;
t107_l_n=t107_l_n_1+t107_e;
}
}
else{
if(t29_h<=0){
REAL t107_l_n_1=t29_l_n*t100_h;
REAL t29_h_n=-t29_h;
REAL t107_h_1=t29_h_n*t100_l_n;
t107_h=t107_h_1+t107_e;
t107_l_n=t107_l_n_1+t107_e;
}
else if(t29_l_n>=0){
REAL t107_l_n_1=t29_l_n*t100_h;
REAL t107_h_1=t100_h*t29_h;
t107_h=t107_h_1+t107_e;
t107_l_n=t107_l_n_1+t107_e;
}
else{
REAL t100_l=-t100_l_n;
REAL t107_l_n_1=t100_l*t29_l_n;
REAL t107_h_1=t100_h*t29_h;
t107_h=t107_h_1+t107_e;
t107_l_n=t107_l_n_1+t107_e;
}
}
REAL t109_e_1=fabs(t109);
REAL t109_e=eps*t109_e_1;
REAL t109_h_1=t106_h+t107_l_n;
REAL t109_h=t109_h_1+t109_e;
REAL t109_l_n_1=t107_h+t106_l_n;
REAL t109_l_n=t109_l_n_1+t109_e;
REAL t110_e_1=fabs(t110);
REAL t110_e=eps*t110_e_1;
REAL  t110_h, t110_l_n;
if(t49_h<=0){
if(t94_h<=0){
REAL t49_h_n=-t49_h;
REAL t110_l_n_1=t49_h_n*t94_h;
REAL t110_h_1=t49_l_n*t94_l_n;
t110_h=t110_h_1+t110_e;
t110_l_n=t110_l_n_1+t110_e;
}
else if(t94_l_n>=0){
REAL t110_l_n_1=t49_l_n*t94_h;
REAL t110_h_1=t49_l_n*t94_l_n;
t110_h=t110_h_1+t110_e;
t110_l_n=t110_l_n_1+t110_e;
}
else{
REAL t110_l_n_1=t49_l_n*t94_h;
REAL t49_h_n=-t49_h;
REAL t110_h_1=t49_h_n*t94_l_n;
t110_h=t110_h_1+t110_e;
t110_l_n=t110_l_n_1+t110_e;
}
}
else if(t49_l_n>=0){
if(t94_h<=0){
REAL t110_l_n_1=t49_h*t94_l_n;
REAL t110_h_1=t49_l_n*t94_l_n;
t110_h=t110_h_1+t110_e;
t110_l_n=t110_l_n_1+t110_e;
}
else if(t94_l_n>=0){
REAL t110_l_n_1_1=t49_l_n*t94_h;
REAL t110_l_n_1_2=t94_l_n*t49_h;
REAL t110_l_n_1=max(t110_l_n_1_1,t110_l_n_1_2);
REAL t110_h_1_1=t49_l_n*t94_l_n;
REAL t110_h_1_2=t49_h*t94_h;
REAL t110_h_1=max(t110_h_1_1,t110_h_1_2);
t110_h=t110_h_1+t110_e;
t110_l_n=t110_l_n_1+t110_e;
}
else{
REAL t110_l_n_1=t49_l_n*t94_h;
REAL t110_h_1=t49_h*t94_h;
t110_h=t110_h_1+t110_e;
t110_l_n=t110_l_n_1+t110_e;
}
}
else{
if(t94_h<=0){
REAL t110_l_n_1=t94_l_n*t49_h;
REAL t94_h_n=-t94_h;
REAL t110_h_1=t94_h_n*t49_l_n;
t110_h=t110_h_1+t110_e;
t110_l_n=t110_l_n_1+t110_e;
}
else if(t94_l_n>=0){
REAL t110_l_n_1=t94_l_n*t49_h;
REAL t110_h_1=t49_h*t94_h;
t110_h=t110_h_1+t110_e;
t110_l_n=t110_l_n_1+t110_e;
}
else{
REAL t49_l=-t49_l_n;
REAL t110_l_n_1=t49_l*t94_l_n;
REAL t110_h_1=t49_h*t94_h;
t110_h=t110_h_1+t110_e;
t110_l_n=t110_l_n_1+t110_e;
}
}
REAL t111_e_1=fabs(t111);
REAL t111_e=eps*t111_e_1;
REAL  t111_h, t111_l_n;
if(t52_h<=0){
if(t93_h<=0){
REAL t52_h_n=-t52_h;
REAL t111_l_n_1=t52_h_n*t93_h;
REAL t111_h_1=t52_l_n*t93_l_n;
t111_h=t111_h_1+t111_e;
t111_l_n=t111_l_n_1+t111_e;
}
else if(t93_l_n>=0){
REAL t111_l_n_1=t52_l_n*t93_h;
REAL t111_h_1=t52_l_n*t93_l_n;
t111_h=t111_h_1+t111_e;
t111_l_n=t111_l_n_1+t111_e;
}
else{
REAL t111_l_n_1=t52_l_n*t93_h;
REAL t52_h_n=-t52_h;
REAL t111_h_1=t52_h_n*t93_l_n;
t111_h=t111_h_1+t111_e;
t111_l_n=t111_l_n_1+t111_e;
}
}
else if(t52_l_n>=0){
if(t93_h<=0){
REAL t111_l_n_1=t52_h*t93_l_n;
REAL t111_h_1=t52_l_n*t93_l_n;
t111_h=t111_h_1+t111_e;
t111_l_n=t111_l_n_1+t111_e;
}
else if(t93_l_n>=0){
REAL t111_l_n_1_1=t52_l_n*t93_h;
REAL t111_l_n_1_2=t93_l_n*t52_h;
REAL t111_l_n_1=max(t111_l_n_1_1,t111_l_n_1_2);
REAL t111_h_1_1=t52_l_n*t93_l_n;
REAL t111_h_1_2=t52_h*t93_h;
REAL t111_h_1=max(t111_h_1_1,t111_h_1_2);
t111_h=t111_h_1+t111_e;
t111_l_n=t111_l_n_1+t111_e;
}
else{
REAL t111_l_n_1=t52_l_n*t93_h;
REAL t111_h_1=t52_h*t93_h;
t111_h=t111_h_1+t111_e;
t111_l_n=t111_l_n_1+t111_e;
}
}
else{
if(t93_h<=0){
REAL t111_l_n_1=t93_l_n*t52_h;
REAL t93_h_n=-t93_h;
REAL t111_h_1=t93_h_n*t52_l_n;
t111_h=t111_h_1+t111_e;
t111_l_n=t111_l_n_1+t111_e;
}
else if(t93_l_n>=0){
REAL t111_l_n_1=t93_l_n*t52_h;
REAL t111_h_1=t52_h*t93_h;
t111_h=t111_h_1+t111_e;
t111_l_n=t111_l_n_1+t111_e;
}
else{
REAL t52_l=-t52_l_n;
REAL t111_l_n_1=t52_l*t93_l_n;
REAL t111_h_1=t52_h*t93_h;
t111_h=t111_h_1+t111_e;
t111_l_n=t111_l_n_1+t111_e;
}
}
REAL t113_e_1=fabs(t113);
REAL t113_e=eps*t113_e_1;
REAL t113_h_1=t110_h+t111_l_n;
REAL t113_h=t113_h_1+t113_e;
REAL t113_l_n_1=t111_h+t110_l_n;
REAL t113_l_n=t113_l_n_1+t113_e;
REAL t115_e_1=fabs(t115);
REAL t115_e=eps*t115_e_1;
REAL t115_h_1=t53_h+t50_l_n;
REAL t115_h=t115_h_1+t115_e;
REAL t115_l_n_1=t50_h+t53_l_n;
REAL t115_l_n=t115_l_n_1+t115_e;
REAL t117_e_1=fabs(t117);
REAL t117_e=eps*t117_e_1;
REAL t117_h_1=t31_h+t28_l_n;
REAL t117_h=t117_h_1+t117_e;
REAL t117_l_n_1=t28_h+t31_l_n;
REAL t117_l_n=t117_l_n_1+t117_e;
REAL t118_e_1=fabs(t118);
REAL t118_e=eps*t118_e_1;
REAL  t118_h, t118_l_n;
if(t59_h<=0){
if(t69_h<=0){
REAL t59_h_n=-t59_h;
REAL t118_l_n_1=t59_h_n*t69_h;
REAL t118_h_1=t59_l_n*t69_l_n;
t118_h=t118_h_1+t118_e;
t118_l_n=t118_l_n_1+t118_e;
}
else if(t69_l_n>=0){
REAL t118_l_n_1=t59_l_n*t69_h;
REAL t118_h_1=t59_l_n*t69_l_n;
t118_h=t118_h_1+t118_e;
t118_l_n=t118_l_n_1+t118_e;
}
else{
REAL t118_l_n_1=t59_l_n*t69_h;
REAL t59_h_n=-t59_h;
REAL t118_h_1=t59_h_n*t69_l_n;
t118_h=t118_h_1+t118_e;
t118_l_n=t118_l_n_1+t118_e;
}
}
else if(t59_l_n>=0){
if(t69_h<=0){
REAL t118_l_n_1=t59_h*t69_l_n;
REAL t118_h_1=t59_l_n*t69_l_n;
t118_h=t118_h_1+t118_e;
t118_l_n=t118_l_n_1+t118_e;
}
else if(t69_l_n>=0){
REAL t118_l_n_1_1=t59_l_n*t69_h;
REAL t118_l_n_1_2=t69_l_n*t59_h;
REAL t118_l_n_1=max(t118_l_n_1_1,t118_l_n_1_2);
REAL t118_h_1_1=t59_l_n*t69_l_n;
REAL t118_h_1_2=t59_h*t69_h;
REAL t118_h_1=max(t118_h_1_1,t118_h_1_2);
t118_h=t118_h_1+t118_e;
t118_l_n=t118_l_n_1+t118_e;
}
else{
REAL t118_l_n_1=t59_l_n*t69_h;
REAL t118_h_1=t59_h*t69_h;
t118_h=t118_h_1+t118_e;
t118_l_n=t118_l_n_1+t118_e;
}
}
else{
if(t69_h<=0){
REAL t118_l_n_1=t69_l_n*t59_h;
REAL t69_h_n=-t69_h;
REAL t118_h_1=t69_h_n*t59_l_n;
t118_h=t118_h_1+t118_e;
t118_l_n=t118_l_n_1+t118_e;
}
else if(t69_l_n>=0){
REAL t118_l_n_1=t69_l_n*t59_h;
REAL t118_h_1=t59_h*t69_h;
t118_h=t118_h_1+t118_e;
t118_l_n=t118_l_n_1+t118_e;
}
else{
REAL t59_l=-t59_l_n;
REAL t118_l_n_1=t59_l*t69_l_n;
REAL t118_h_1=t59_h*t69_h;
t118_h=t118_h_1+t118_e;
t118_l_n=t118_l_n_1+t118_e;
}
}
REAL t119_e_1=fabs(t119);
REAL t119_e=eps*t119_e_1;
REAL  t119_h, t119_l_n;
if(t63_h<=0){
if(t71_h<=0){
REAL t63_h_n=-t63_h;
REAL t119_l_n_1=t63_h_n*t71_h;
REAL t119_h_1=t63_l_n*t71_l_n;
t119_h=t119_h_1+t119_e;
t119_l_n=t119_l_n_1+t119_e;
}
else if(t71_l_n>=0){
REAL t119_l_n_1=t63_l_n*t71_h;
REAL t119_h_1=t63_l_n*t71_l_n;
t119_h=t119_h_1+t119_e;
t119_l_n=t119_l_n_1+t119_e;
}
else{
REAL t119_l_n_1=t63_l_n*t71_h;
REAL t63_h_n=-t63_h;
REAL t119_h_1=t63_h_n*t71_l_n;
t119_h=t119_h_1+t119_e;
t119_l_n=t119_l_n_1+t119_e;
}
}
else if(t63_l_n>=0){
if(t71_h<=0){
REAL t119_l_n_1=t63_h*t71_l_n;
REAL t119_h_1=t63_l_n*t71_l_n;
t119_h=t119_h_1+t119_e;
t119_l_n=t119_l_n_1+t119_e;
}
else if(t71_l_n>=0){
REAL t119_l_n_1_1=t63_l_n*t71_h;
REAL t119_l_n_1_2=t71_l_n*t63_h;
REAL t119_l_n_1=max(t119_l_n_1_1,t119_l_n_1_2);
REAL t119_h_1_1=t63_l_n*t71_l_n;
REAL t119_h_1_2=t63_h*t71_h;
REAL t119_h_1=max(t119_h_1_1,t119_h_1_2);
t119_h=t119_h_1+t119_e;
t119_l_n=t119_l_n_1+t119_e;
}
else{
REAL t119_l_n_1=t63_l_n*t71_h;
REAL t119_h_1=t63_h*t71_h;
t119_h=t119_h_1+t119_e;
t119_l_n=t119_l_n_1+t119_e;
}
}
else{
if(t71_h<=0){
REAL t119_l_n_1=t71_l_n*t63_h;
REAL t71_h_n=-t71_h;
REAL t119_h_1=t71_h_n*t63_l_n;
t119_h=t119_h_1+t119_e;
t119_l_n=t119_l_n_1+t119_e;
}
else if(t71_l_n>=0){
REAL t119_l_n_1=t71_l_n*t63_h;
REAL t119_h_1=t63_h*t71_h;
t119_h=t119_h_1+t119_e;
t119_l_n=t119_l_n_1+t119_e;
}
else{
REAL t63_l=-t63_l_n;
REAL t119_l_n_1=t63_l*t71_l_n;
REAL t119_h_1=t63_h*t71_h;
t119_h=t119_h_1+t119_e;
t119_l_n=t119_l_n_1+t119_e;
}
}
REAL t121_e_1=fabs(t121);
REAL t121_e=eps*t121_e_1;
REAL t121_h_1=t118_h+t119_l_n;
REAL t121_h=t121_h_1+t121_e;
REAL t121_l_n_1=t119_h+t118_l_n;
REAL t121_l_n=t121_l_n_1+t121_e;
REAL t122_e_1=fabs(t122);
REAL t122_e=eps*t122_e_1;
REAL  t122_h, t122_l_n;
if(t121_h<=0){
if(t99_h<=0){
REAL t121_h_n=-t121_h;
REAL t122_l_n_1=t121_h_n*t99_h;
REAL t122_h_1=t121_l_n*t99_l_n;
t122_h=t122_h_1+t122_e;
t122_l_n=t122_l_n_1+t122_e;
}
else if(t99_l_n>=0){
REAL t122_l_n_1=t121_l_n*t99_h;
REAL t122_h_1=t121_l_n*t99_l_n;
t122_h=t122_h_1+t122_e;
t122_l_n=t122_l_n_1+t122_e;
}
else{
REAL t122_l_n_1=t121_l_n*t99_h;
REAL t121_h_n=-t121_h;
REAL t122_h_1=t121_h_n*t99_l_n;
t122_h=t122_h_1+t122_e;
t122_l_n=t122_l_n_1+t122_e;
}
}
else if(t121_l_n>=0){
if(t99_h<=0){
REAL t122_l_n_1=t121_h*t99_l_n;
REAL t122_h_1=t121_l_n*t99_l_n;
t122_h=t122_h_1+t122_e;
t122_l_n=t122_l_n_1+t122_e;
}
else if(t99_l_n>=0){
REAL t122_l_n_1_1=t121_l_n*t99_h;
REAL t122_l_n_1_2=t99_l_n*t121_h;
REAL t122_l_n_1=max(t122_l_n_1_1,t122_l_n_1_2);
REAL t122_h_1_1=t121_l_n*t99_l_n;
REAL t122_h_1_2=t121_h*t99_h;
REAL t122_h_1=max(t122_h_1_1,t122_h_1_2);
t122_h=t122_h_1+t122_e;
t122_l_n=t122_l_n_1+t122_e;
}
else{
REAL t122_l_n_1=t121_l_n*t99_h;
REAL t122_h_1=t121_h*t99_h;
t122_h=t122_h_1+t122_e;
t122_l_n=t122_l_n_1+t122_e;
}
}
else{
if(t99_h<=0){
REAL t122_l_n_1=t99_l_n*t121_h;
REAL t99_h_n=-t99_h;
REAL t122_h_1=t99_h_n*t121_l_n;
t122_h=t122_h_1+t122_e;
t122_l_n=t122_l_n_1+t122_e;
}
else if(t99_l_n>=0){
REAL t122_l_n_1=t99_l_n*t121_h;
REAL t122_h_1=t121_h*t99_h;
t122_h=t122_h_1+t122_e;
t122_l_n=t122_l_n_1+t122_e;
}
else{
REAL t121_l=-t121_l_n;
REAL t122_l_n_1=t121_l*t99_l_n;
REAL t122_h_1=t121_h*t99_h;
t122_h=t122_h_1+t122_e;
t122_l_n=t122_l_n_1+t122_e;
}
}
REAL t123_e_1=fabs(t123);
REAL t123_e=eps*t123_e_1;
REAL t123_h_1=t122_h+t78_h;
REAL t123_h=t123_h_1+t123_e;
REAL t123_l_n_1=t122_l_n+t78_l_n;
REAL t123_l_n=t123_l_n_1+t123_e;
REAL t124_e_1=fabs(t124);
REAL t124_e=eps*t124_e_1;
REAL  t124_h, t124_l_n;
if(t123_h<=0){
if(t45_h<=0){
REAL t123_h_n=-t123_h;
REAL t124_l_n_1=t123_h_n*t45_h;
REAL t124_h_1=t123_l_n*t45_l_n;
t124_h=t124_h_1+t124_e;
t124_l_n=t124_l_n_1+t124_e;
}
else if(t45_l_n>=0){
REAL t124_l_n_1=t123_l_n*t45_h;
REAL t124_h_1=t123_l_n*t45_l_n;
t124_h=t124_h_1+t124_e;
t124_l_n=t124_l_n_1+t124_e;
}
else{
REAL t124_l_n_1=t123_l_n*t45_h;
REAL t123_h_n=-t123_h;
REAL t124_h_1=t123_h_n*t45_l_n;
t124_h=t124_h_1+t124_e;
t124_l_n=t124_l_n_1+t124_e;
}
}
else if(t123_l_n>=0){
if(t45_h<=0){
REAL t124_l_n_1=t123_h*t45_l_n;
REAL t124_h_1=t123_l_n*t45_l_n;
t124_h=t124_h_1+t124_e;
t124_l_n=t124_l_n_1+t124_e;
}
else if(t45_l_n>=0){
REAL t124_l_n_1_1=t123_l_n*t45_h;
REAL t124_l_n_1_2=t45_l_n*t123_h;
REAL t124_l_n_1=max(t124_l_n_1_1,t124_l_n_1_2);
REAL t124_h_1_1=t123_l_n*t45_l_n;
REAL t124_h_1_2=t123_h*t45_h;
REAL t124_h_1=max(t124_h_1_1,t124_h_1_2);
t124_h=t124_h_1+t124_e;
t124_l_n=t124_l_n_1+t124_e;
}
else{
REAL t124_l_n_1=t123_l_n*t45_h;
REAL t124_h_1=t123_h*t45_h;
t124_h=t124_h_1+t124_e;
t124_l_n=t124_l_n_1+t124_e;
}
}
else{
if(t45_h<=0){
REAL t124_l_n_1=t45_l_n*t123_h;
REAL t45_h_n=-t45_h;
REAL t124_h_1=t45_h_n*t123_l_n;
t124_h=t124_h_1+t124_e;
t124_l_n=t124_l_n_1+t124_e;
}
else if(t45_l_n>=0){
REAL t124_l_n_1=t45_l_n*t123_h;
REAL t124_h_1=t123_h*t45_h;
t124_h=t124_h_1+t124_e;
t124_l_n=t124_l_n_1+t124_e;
}
else{
REAL t123_l=-t123_l_n;
REAL t124_l_n_1=t123_l*t45_l_n;
REAL t124_h_1=t123_h*t45_h;
t124_h=t124_h_1+t124_e;
t124_l_n=t124_l_n_1+t124_e;
}
}
REAL t125_e_1=fabs(t125);
REAL t125_e=eps*t125_e_1;
REAL  t125_h, t125_l_n;
if(t37_h<=0){
if(t80_h<=0){
REAL t37_h_n=-t37_h;
REAL t125_l_n_1=t37_h_n*t80_h;
REAL t125_h_1=t37_l_n*t80_l_n;
t125_h=t125_h_1+t125_e;
t125_l_n=t125_l_n_1+t125_e;
}
else if(t80_l_n>=0){
REAL t125_l_n_1=t37_l_n*t80_h;
REAL t125_h_1=t37_l_n*t80_l_n;
t125_h=t125_h_1+t125_e;
t125_l_n=t125_l_n_1+t125_e;
}
else{
REAL t125_l_n_1=t37_l_n*t80_h;
REAL t37_h_n=-t37_h;
REAL t125_h_1=t37_h_n*t80_l_n;
t125_h=t125_h_1+t125_e;
t125_l_n=t125_l_n_1+t125_e;
}
}
else if(t37_l_n>=0){
if(t80_h<=0){
REAL t125_l_n_1=t37_h*t80_l_n;
REAL t125_h_1=t37_l_n*t80_l_n;
t125_h=t125_h_1+t125_e;
t125_l_n=t125_l_n_1+t125_e;
}
else if(t80_l_n>=0){
REAL t125_l_n_1_1=t37_l_n*t80_h;
REAL t125_l_n_1_2=t80_l_n*t37_h;
REAL t125_l_n_1=max(t125_l_n_1_1,t125_l_n_1_2);
REAL t125_h_1_1=t37_l_n*t80_l_n;
REAL t125_h_1_2=t37_h*t80_h;
REAL t125_h_1=max(t125_h_1_1,t125_h_1_2);
t125_h=t125_h_1+t125_e;
t125_l_n=t125_l_n_1+t125_e;
}
else{
REAL t125_l_n_1=t37_l_n*t80_h;
REAL t125_h_1=t37_h*t80_h;
t125_h=t125_h_1+t125_e;
t125_l_n=t125_l_n_1+t125_e;
}
}
else{
if(t80_h<=0){
REAL t125_l_n_1=t80_l_n*t37_h;
REAL t80_h_n=-t80_h;
REAL t125_h_1=t80_h_n*t37_l_n;
t125_h=t125_h_1+t125_e;
t125_l_n=t125_l_n_1+t125_e;
}
else if(t80_l_n>=0){
REAL t125_l_n_1=t80_l_n*t37_h;
REAL t125_h_1=t37_h*t80_h;
t125_h=t125_h_1+t125_e;
t125_l_n=t125_l_n_1+t125_e;
}
else{
REAL t37_l=-t37_l_n;
REAL t125_l_n_1=t37_l*t80_l_n;
REAL t125_h_1=t37_h*t80_h;
t125_h=t125_h_1+t125_e;
t125_l_n=t125_l_n_1+t125_e;
}
}
REAL t126_e_1=fabs(t126);
REAL t126_e=eps*t126_e_1;
REAL  t126_h, t126_l_n;
if(t41_h<=0){
if(t82_h<=0){
REAL t41_h_n=-t41_h;
REAL t126_l_n_1=t41_h_n*t82_h;
REAL t126_h_1=t41_l_n*t82_l_n;
t126_h=t126_h_1+t126_e;
t126_l_n=t126_l_n_1+t126_e;
}
else if(t82_l_n>=0){
REAL t126_l_n_1=t41_l_n*t82_h;
REAL t126_h_1=t41_l_n*t82_l_n;
t126_h=t126_h_1+t126_e;
t126_l_n=t126_l_n_1+t126_e;
}
else{
REAL t126_l_n_1=t41_l_n*t82_h;
REAL t41_h_n=-t41_h;
REAL t126_h_1=t41_h_n*t82_l_n;
t126_h=t126_h_1+t126_e;
t126_l_n=t126_l_n_1+t126_e;
}
}
else if(t41_l_n>=0){
if(t82_h<=0){
REAL t126_l_n_1=t41_h*t82_l_n;
REAL t126_h_1=t41_l_n*t82_l_n;
t126_h=t126_h_1+t126_e;
t126_l_n=t126_l_n_1+t126_e;
}
else if(t82_l_n>=0){
REAL t126_l_n_1_1=t41_l_n*t82_h;
REAL t126_l_n_1_2=t82_l_n*t41_h;
REAL t126_l_n_1=max(t126_l_n_1_1,t126_l_n_1_2);
REAL t126_h_1_1=t41_l_n*t82_l_n;
REAL t126_h_1_2=t41_h*t82_h;
REAL t126_h_1=max(t126_h_1_1,t126_h_1_2);
t126_h=t126_h_1+t126_e;
t126_l_n=t126_l_n_1+t126_e;
}
else{
REAL t126_l_n_1=t41_l_n*t82_h;
REAL t126_h_1=t41_h*t82_h;
t126_h=t126_h_1+t126_e;
t126_l_n=t126_l_n_1+t126_e;
}
}
else{
if(t82_h<=0){
REAL t126_l_n_1=t82_l_n*t41_h;
REAL t82_h_n=-t82_h;
REAL t126_h_1=t82_h_n*t41_l_n;
t126_h=t126_h_1+t126_e;
t126_l_n=t126_l_n_1+t126_e;
}
else if(t82_l_n>=0){
REAL t126_l_n_1=t82_l_n*t41_h;
REAL t126_h_1=t41_h*t82_h;
t126_h=t126_h_1+t126_e;
t126_l_n=t126_l_n_1+t126_e;
}
else{
REAL t41_l=-t41_l_n;
REAL t126_l_n_1=t41_l*t82_l_n;
REAL t126_h_1=t41_h*t82_h;
t126_h=t126_h_1+t126_e;
t126_l_n=t126_l_n_1+t126_e;
}
}
REAL t128_e_1=fabs(t128);
REAL t128_e=eps*t128_e_1;
REAL t128_h_1=t125_h+t126_l_n;
REAL t128_h=t128_h_1+t128_e;
REAL t128_l_n_1=t126_h+t125_l_n;
REAL t128_l_n=t128_l_n_1+t128_e;
REAL t129_e_1=fabs(t129);
REAL t129_e=eps*t129_e_1;
REAL  t129_h, t129_l_n;
if(t128_h<=0){
if(t92_h<=0){
REAL t128_h_n=-t128_h;
REAL t129_l_n_1=t128_h_n*t92_h;
REAL t129_h_1=t128_l_n*t92_l_n;
t129_h=t129_h_1+t129_e;
t129_l_n=t129_l_n_1+t129_e;
}
else if(t92_l_n>=0){
REAL t129_l_n_1=t128_l_n*t92_h;
REAL t129_h_1=t128_l_n*t92_l_n;
t129_h=t129_h_1+t129_e;
t129_l_n=t129_l_n_1+t129_e;
}
else{
REAL t129_l_n_1=t128_l_n*t92_h;
REAL t128_h_n=-t128_h;
REAL t129_h_1=t128_h_n*t92_l_n;
t129_h=t129_h_1+t129_e;
t129_l_n=t129_l_n_1+t129_e;
}
}
else if(t128_l_n>=0){
if(t92_h<=0){
REAL t129_l_n_1=t128_h*t92_l_n;
REAL t129_h_1=t128_l_n*t92_l_n;
t129_h=t129_h_1+t129_e;
t129_l_n=t129_l_n_1+t129_e;
}
else if(t92_l_n>=0){
REAL t129_l_n_1_1=t128_l_n*t92_h;
REAL t129_l_n_1_2=t92_l_n*t128_h;
REAL t129_l_n_1=max(t129_l_n_1_1,t129_l_n_1_2);
REAL t129_h_1_1=t128_l_n*t92_l_n;
REAL t129_h_1_2=t128_h*t92_h;
REAL t129_h_1=max(t129_h_1_1,t129_h_1_2);
t129_h=t129_h_1+t129_e;
t129_l_n=t129_l_n_1+t129_e;
}
else{
REAL t129_l_n_1=t128_l_n*t92_h;
REAL t129_h_1=t128_h*t92_h;
t129_h=t129_h_1+t129_e;
t129_l_n=t129_l_n_1+t129_e;
}
}
else{
if(t92_h<=0){
REAL t129_l_n_1=t92_l_n*t128_h;
REAL t92_h_n=-t92_h;
REAL t129_h_1=t92_h_n*t128_l_n;
t129_h=t129_h_1+t129_e;
t129_l_n=t129_l_n_1+t129_e;
}
else if(t92_l_n>=0){
REAL t129_l_n_1=t92_l_n*t128_h;
REAL t129_h_1=t128_h*t92_h;
t129_h=t129_h_1+t129_e;
t129_l_n=t129_l_n_1+t129_e;
}
else{
REAL t128_l=-t128_l_n;
REAL t129_l_n_1=t128_l*t92_l_n;
REAL t129_h_1=t128_h*t92_h;
t129_h=t129_h_1+t129_e;
t129_l_n=t129_l_n_1+t129_e;
}
}
REAL t130_e_1=fabs(t130);
REAL t130_e=eps*t130_e_1;
REAL t130_h_1=t129_h+t89_h;
REAL t130_h=t130_h_1+t130_e;
REAL t130_l_n_1=t129_l_n+t89_l_n;
REAL t130_l_n=t130_l_n_1+t130_e;
REAL t131_e_1=fabs(t131);
REAL t131_e=eps*t131_e_1;
REAL  t131_h, t131_l_n;
if(t130_h<=0){
if(t67_h<=0){
REAL t130_h_n=-t130_h;
REAL t131_l_n_1=t130_h_n*t67_h;
REAL t131_h_1=t130_l_n*t67_l_n;
t131_h=t131_h_1+t131_e;
t131_l_n=t131_l_n_1+t131_e;
}
else if(t67_l_n>=0){
REAL t131_l_n_1=t130_l_n*t67_h;
REAL t131_h_1=t130_l_n*t67_l_n;
t131_h=t131_h_1+t131_e;
t131_l_n=t131_l_n_1+t131_e;
}
else{
REAL t131_l_n_1=t130_l_n*t67_h;
REAL t130_h_n=-t130_h;
REAL t131_h_1=t130_h_n*t67_l_n;
t131_h=t131_h_1+t131_e;
t131_l_n=t131_l_n_1+t131_e;
}
}
else if(t130_l_n>=0){
if(t67_h<=0){
REAL t131_l_n_1=t130_h*t67_l_n;
REAL t131_h_1=t130_l_n*t67_l_n;
t131_h=t131_h_1+t131_e;
t131_l_n=t131_l_n_1+t131_e;
}
else if(t67_l_n>=0){
REAL t131_l_n_1_1=t130_l_n*t67_h;
REAL t131_l_n_1_2=t67_l_n*t130_h;
REAL t131_l_n_1=max(t131_l_n_1_1,t131_l_n_1_2);
REAL t131_h_1_1=t130_l_n*t67_l_n;
REAL t131_h_1_2=t130_h*t67_h;
REAL t131_h_1=max(t131_h_1_1,t131_h_1_2);
t131_h=t131_h_1+t131_e;
t131_l_n=t131_l_n_1+t131_e;
}
else{
REAL t131_l_n_1=t130_l_n*t67_h;
REAL t131_h_1=t130_h*t67_h;
t131_h=t131_h_1+t131_e;
t131_l_n=t131_l_n_1+t131_e;
}
}
else{
if(t67_h<=0){
REAL t131_l_n_1=t67_l_n*t130_h;
REAL t67_h_n=-t67_h;
REAL t131_h_1=t67_h_n*t130_l_n;
t131_h=t131_h_1+t131_e;
t131_l_n=t131_l_n_1+t131_e;
}
else if(t67_l_n>=0){
REAL t131_l_n_1=t67_l_n*t130_h;
REAL t131_h_1=t130_h*t67_h;
t131_h=t131_h_1+t131_e;
t131_l_n=t131_l_n_1+t131_e;
}
else{
REAL t130_l=-t130_l_n;
REAL t131_l_n_1=t130_l*t67_l_n;
REAL t131_h_1=t130_h*t67_h;
t131_h=t131_h_1+t131_e;
t131_l_n=t131_l_n_1+t131_e;
}
}
REAL t133_e_1=fabs(t133);
REAL t133_e=eps*t133_e_1;
REAL t133_h_1=t124_h+t131_l_n;
REAL t133_h=t133_h_1+t133_e;
REAL t133_l_n_1=t131_h+t124_l_n;
REAL t133_l_n=t133_l_n_1+t133_e;
REAL t134_e_1=fabs(t134);
REAL t134_e=eps*t134_e_1;
REAL  t134_h, t134_l_n;
if(t133_h<=0){
if(t23_h<=0){
REAL t133_h_n=-t133_h;
REAL t134_l_n_1=t133_h_n*t23_h;
REAL t134_h_1=t133_l_n*t23_l_n;
t134_h=t134_h_1+t134_e;
t134_l_n=t134_l_n_1+t134_e;
}
else if(t23_l_n>=0){
REAL t134_l_n_1=t133_l_n*t23_h;
REAL t134_h_1=t133_l_n*t23_l_n;
t134_h=t134_h_1+t134_e;
t134_l_n=t134_l_n_1+t134_e;
}
else{
REAL t134_l_n_1=t133_l_n*t23_h;
REAL t133_h_n=-t133_h;
REAL t134_h_1=t133_h_n*t23_l_n;
t134_h=t134_h_1+t134_e;
t134_l_n=t134_l_n_1+t134_e;
}
}
else if(t133_l_n>=0){
if(t23_h<=0){
REAL t134_l_n_1=t133_h*t23_l_n;
REAL t134_h_1=t133_l_n*t23_l_n;
t134_h=t134_h_1+t134_e;
t134_l_n=t134_l_n_1+t134_e;
}
else if(t23_l_n>=0){
REAL t134_l_n_1_1=t133_l_n*t23_h;
REAL t134_l_n_1_2=t23_l_n*t133_h;
REAL t134_l_n_1=max(t134_l_n_1_1,t134_l_n_1_2);
REAL t134_h_1_1=t133_l_n*t23_l_n;
REAL t134_h_1_2=t133_h*t23_h;
REAL t134_h_1=max(t134_h_1_1,t134_h_1_2);
t134_h=t134_h_1+t134_e;
t134_l_n=t134_l_n_1+t134_e;
}
else{
REAL t134_l_n_1=t133_l_n*t23_h;
REAL t134_h_1=t133_h*t23_h;
t134_h=t134_h_1+t134_e;
t134_l_n=t134_l_n_1+t134_e;
}
}
else{
if(t23_h<=0){
REAL t134_l_n_1=t23_l_n*t133_h;
REAL t23_h_n=-t23_h;
REAL t134_h_1=t23_h_n*t133_l_n;
t134_h=t134_h_1+t134_e;
t134_l_n=t134_l_n_1+t134_e;
}
else if(t23_l_n>=0){
REAL t134_l_n_1=t23_l_n*t133_h;
REAL t134_h_1=t133_h*t23_h;
t134_h=t134_h_1+t134_e;
t134_l_n=t134_l_n_1+t134_e;
}
else{
REAL t133_l=-t133_l_n;
REAL t134_l_n_1=t133_l*t23_l_n;
REAL t134_h_1=t133_h*t23_h;
t134_h=t134_h_1+t134_e;
t134_l_n=t134_l_n_1+t134_e;
}
}
REAL t135_e_1=fabs(t135);
REAL t135_e=eps*t135_e_1;
REAL  t135_h, t135_l_n;
if(t33_h<=0){
if(t99_h<=0){
REAL t33_h_n=-t33_h;
REAL t135_l_n_1=t33_h_n*t99_h;
REAL t135_h_1=t33_l_n*t99_l_n;
t135_h=t135_h_1+t135_e;
t135_l_n=t135_l_n_1+t135_e;
}
else if(t99_l_n>=0){
REAL t135_l_n_1=t33_l_n*t99_h;
REAL t135_h_1=t33_l_n*t99_l_n;
t135_h=t135_h_1+t135_e;
t135_l_n=t135_l_n_1+t135_e;
}
else{
REAL t135_l_n_1=t33_l_n*t99_h;
REAL t33_h_n=-t33_h;
REAL t135_h_1=t33_h_n*t99_l_n;
t135_h=t135_h_1+t135_e;
t135_l_n=t135_l_n_1+t135_e;
}
}
else if(t33_l_n>=0){
if(t99_h<=0){
REAL t135_l_n_1=t33_h*t99_l_n;
REAL t135_h_1=t33_l_n*t99_l_n;
t135_h=t135_h_1+t135_e;
t135_l_n=t135_l_n_1+t135_e;
}
else if(t99_l_n>=0){
REAL t135_l_n_1_1=t33_l_n*t99_h;
REAL t135_l_n_1_2=t99_l_n*t33_h;
REAL t135_l_n_1=max(t135_l_n_1_1,t135_l_n_1_2);
REAL t135_h_1_1=t33_l_n*t99_l_n;
REAL t135_h_1_2=t33_h*t99_h;
REAL t135_h_1=max(t135_h_1_1,t135_h_1_2);
t135_h=t135_h_1+t135_e;
t135_l_n=t135_l_n_1+t135_e;
}
else{
REAL t135_l_n_1=t33_l_n*t99_h;
REAL t135_h_1=t33_h*t99_h;
t135_h=t135_h_1+t135_e;
t135_l_n=t135_l_n_1+t135_e;
}
}
else{
if(t99_h<=0){
REAL t135_l_n_1=t99_l_n*t33_h;
REAL t99_h_n=-t99_h;
REAL t135_h_1=t99_h_n*t33_l_n;
t135_h=t135_h_1+t135_e;
t135_l_n=t135_l_n_1+t135_e;
}
else if(t99_l_n>=0){
REAL t135_l_n_1=t99_l_n*t33_h;
REAL t135_h_1=t33_h*t99_h;
t135_h=t135_h_1+t135_e;
t135_l_n=t135_l_n_1+t135_e;
}
else{
REAL t33_l=-t33_l_n;
REAL t135_l_n_1=t33_l*t99_l_n;
REAL t135_h_1=t33_h*t99_h;
t135_h=t135_h_1+t135_e;
t135_l_n=t135_l_n_1+t135_e;
}
}
REAL t136_e_1=fabs(t136);
REAL t136_e=eps*t136_e_1;
REAL t136_h_1=t135_h+t78_h;
REAL t136_h=t136_h_1+t136_e;
REAL t136_l_n_1=t135_l_n+t78_l_n;
REAL t136_l_n=t136_l_n_1+t136_e;
REAL t137_e_1=fabs(t137);
REAL t137_e=eps*t137_e_1;
REAL  t137_h, t137_l_n;
if(t136_h<=0){
if(t55_h<=0){
REAL t136_h_n=-t136_h;
REAL t137_l_n_1=t136_h_n*t55_h;
REAL t137_h_1=t136_l_n*t55_l_n;
t137_h=t137_h_1+t137_e;
t137_l_n=t137_l_n_1+t137_e;
}
else if(t55_l_n>=0){
REAL t137_l_n_1=t136_l_n*t55_h;
REAL t137_h_1=t136_l_n*t55_l_n;
t137_h=t137_h_1+t137_e;
t137_l_n=t137_l_n_1+t137_e;
}
else{
REAL t137_l_n_1=t136_l_n*t55_h;
REAL t136_h_n=-t136_h;
REAL t137_h_1=t136_h_n*t55_l_n;
t137_h=t137_h_1+t137_e;
t137_l_n=t137_l_n_1+t137_e;
}
}
else if(t136_l_n>=0){
if(t55_h<=0){
REAL t137_l_n_1=t136_h*t55_l_n;
REAL t137_h_1=t136_l_n*t55_l_n;
t137_h=t137_h_1+t137_e;
t137_l_n=t137_l_n_1+t137_e;
}
else if(t55_l_n>=0){
REAL t137_l_n_1_1=t136_l_n*t55_h;
REAL t137_l_n_1_2=t55_l_n*t136_h;
REAL t137_l_n_1=max(t137_l_n_1_1,t137_l_n_1_2);
REAL t137_h_1_1=t136_l_n*t55_l_n;
REAL t137_h_1_2=t136_h*t55_h;
REAL t137_h_1=max(t137_h_1_1,t137_h_1_2);
t137_h=t137_h_1+t137_e;
t137_l_n=t137_l_n_1+t137_e;
}
else{
REAL t137_l_n_1=t136_l_n*t55_h;
REAL t137_h_1=t136_h*t55_h;
t137_h=t137_h_1+t137_e;
t137_l_n=t137_l_n_1+t137_e;
}
}
else{
if(t55_h<=0){
REAL t137_l_n_1=t55_l_n*t136_h;
REAL t55_h_n=-t55_h;
REAL t137_h_1=t55_h_n*t136_l_n;
t137_h=t137_h_1+t137_e;
t137_l_n=t137_l_n_1+t137_e;
}
else if(t55_l_n>=0){
REAL t137_l_n_1=t55_l_n*t136_h;
REAL t137_h_1=t136_h*t55_h;
t137_h=t137_h_1+t137_e;
t137_l_n=t137_l_n_1+t137_e;
}
else{
REAL t136_l=-t136_l_n;
REAL t137_l_n_1=t136_l*t55_l_n;
REAL t137_h_1=t136_h*t55_h;
t137_h=t137_h_1+t137_e;
t137_l_n=t137_l_n_1+t137_e;
}
}
REAL t138_e_1=fabs(t138);
REAL t138_e=eps*t138_e_1;
REAL  t138_h, t138_l_n;
if(t55_h<=0){
if(t92_h<=0){
REAL t55_h_n=-t55_h;
REAL t138_l_n_1=t55_h_n*t92_h;
REAL t138_h_1=t55_l_n*t92_l_n;
t138_h=t138_h_1+t138_e;
t138_l_n=t138_l_n_1+t138_e;
}
else if(t92_l_n>=0){
REAL t138_l_n_1=t55_l_n*t92_h;
REAL t138_h_1=t55_l_n*t92_l_n;
t138_h=t138_h_1+t138_e;
t138_l_n=t138_l_n_1+t138_e;
}
else{
REAL t138_l_n_1=t55_l_n*t92_h;
REAL t55_h_n=-t55_h;
REAL t138_h_1=t55_h_n*t92_l_n;
t138_h=t138_h_1+t138_e;
t138_l_n=t138_l_n_1+t138_e;
}
}
else if(t55_l_n>=0){
if(t92_h<=0){
REAL t138_l_n_1=t55_h*t92_l_n;
REAL t138_h_1=t55_l_n*t92_l_n;
t138_h=t138_h_1+t138_e;
t138_l_n=t138_l_n_1+t138_e;
}
else if(t92_l_n>=0){
REAL t138_l_n_1_1=t55_l_n*t92_h;
REAL t138_l_n_1_2=t92_l_n*t55_h;
REAL t138_l_n_1=max(t138_l_n_1_1,t138_l_n_1_2);
REAL t138_h_1_1=t55_l_n*t92_l_n;
REAL t138_h_1_2=t55_h*t92_h;
REAL t138_h_1=max(t138_h_1_1,t138_h_1_2);
t138_h=t138_h_1+t138_e;
t138_l_n=t138_l_n_1+t138_e;
}
else{
REAL t138_l_n_1=t55_l_n*t92_h;
REAL t138_h_1=t55_h*t92_h;
t138_h=t138_h_1+t138_e;
t138_l_n=t138_l_n_1+t138_e;
}
}
else{
if(t92_h<=0){
REAL t138_l_n_1=t92_l_n*t55_h;
REAL t92_h_n=-t92_h;
REAL t138_h_1=t92_h_n*t55_l_n;
t138_h=t138_h_1+t138_e;
t138_l_n=t138_l_n_1+t138_e;
}
else if(t92_l_n>=0){
REAL t138_l_n_1=t92_l_n*t55_h;
REAL t138_h_1=t55_h*t92_h;
t138_h=t138_h_1+t138_e;
t138_l_n=t138_l_n_1+t138_e;
}
else{
REAL t55_l=-t55_l_n;
REAL t138_l_n_1=t55_l*t92_l_n;
REAL t138_h_1=t55_h*t92_h;
t138_h=t138_h_1+t138_e;
t138_l_n=t138_l_n_1+t138_e;
}
}
REAL t139_e_1=fabs(t139);
REAL t139_e=eps*t139_e_1;
REAL t139_h_1=t138_h+t89_h;
REAL t139_h=t139_h_1+t139_e;
REAL t139_l_n_1=t138_l_n+t89_l_n;
REAL t139_l_n=t139_l_n_1+t139_e;
REAL t140_e_1=fabs(t140);
REAL t140_e=eps*t140_e_1;
REAL  t140_h, t140_l_n;
if(t139_h<=0){
if(t33_h<=0){
REAL t139_h_n=-t139_h;
REAL t140_l_n_1=t139_h_n*t33_h;
REAL t140_h_1=t139_l_n*t33_l_n;
t140_h=t140_h_1+t140_e;
t140_l_n=t140_l_n_1+t140_e;
}
else if(t33_l_n>=0){
REAL t140_l_n_1=t139_l_n*t33_h;
REAL t140_h_1=t139_l_n*t33_l_n;
t140_h=t140_h_1+t140_e;
t140_l_n=t140_l_n_1+t140_e;
}
else{
REAL t140_l_n_1=t139_l_n*t33_h;
REAL t139_h_n=-t139_h;
REAL t140_h_1=t139_h_n*t33_l_n;
t140_h=t140_h_1+t140_e;
t140_l_n=t140_l_n_1+t140_e;
}
}
else if(t139_l_n>=0){
if(t33_h<=0){
REAL t140_l_n_1=t139_h*t33_l_n;
REAL t140_h_1=t139_l_n*t33_l_n;
t140_h=t140_h_1+t140_e;
t140_l_n=t140_l_n_1+t140_e;
}
else if(t33_l_n>=0){
REAL t140_l_n_1_1=t139_l_n*t33_h;
REAL t140_l_n_1_2=t33_l_n*t139_h;
REAL t140_l_n_1=max(t140_l_n_1_1,t140_l_n_1_2);
REAL t140_h_1_1=t139_l_n*t33_l_n;
REAL t140_h_1_2=t139_h*t33_h;
REAL t140_h_1=max(t140_h_1_1,t140_h_1_2);
t140_h=t140_h_1+t140_e;
t140_l_n=t140_l_n_1+t140_e;
}
else{
REAL t140_l_n_1=t139_l_n*t33_h;
REAL t140_h_1=t139_h*t33_h;
t140_h=t140_h_1+t140_e;
t140_l_n=t140_l_n_1+t140_e;
}
}
else{
if(t33_h<=0){
REAL t140_l_n_1=t33_l_n*t139_h;
REAL t33_h_n=-t33_h;
REAL t140_h_1=t33_h_n*t139_l_n;
t140_h=t140_h_1+t140_e;
t140_l_n=t140_l_n_1+t140_e;
}
else if(t33_l_n>=0){
REAL t140_l_n_1=t33_l_n*t139_h;
REAL t140_h_1=t139_h*t33_h;
t140_h=t140_h_1+t140_e;
t140_l_n=t140_l_n_1+t140_e;
}
else{
REAL t139_l=-t139_l_n;
REAL t140_l_n_1=t139_l*t33_l_n;
REAL t140_h_1=t139_h*t33_h;
t140_h=t140_h_1+t140_e;
t140_l_n=t140_l_n_1+t140_e;
}
}
REAL t142_e_1=fabs(t142);
REAL t142_e=eps*t142_e_1;
REAL t142_h_1=t137_h+t140_l_n;
REAL t142_h=t142_h_1+t142_e;
REAL t142_l_n_1=t140_h+t137_l_n;
REAL t142_l_n=t142_l_n_1+t142_e;
REAL t143_e_1=fabs(t143);
REAL t143_e=eps*t143_e_1;
REAL  t143_h, t143_l_n;
if(t11_h<=0){
if(t142_h<=0){
REAL t11_h_n=-t11_h;
REAL t143_l_n_1=t11_h_n*t142_h;
REAL t143_h_1=t11_l_n*t142_l_n;
t143_h=t143_h_1+t143_e;
t143_l_n=t143_l_n_1+t143_e;
}
else if(t142_l_n>=0){
REAL t143_l_n_1=t11_l_n*t142_h;
REAL t143_h_1=t11_l_n*t142_l_n;
t143_h=t143_h_1+t143_e;
t143_l_n=t143_l_n_1+t143_e;
}
else{
REAL t143_l_n_1=t11_l_n*t142_h;
REAL t11_h_n=-t11_h;
REAL t143_h_1=t11_h_n*t142_l_n;
t143_h=t143_h_1+t143_e;
t143_l_n=t143_l_n_1+t143_e;
}
}
else if(t11_l_n>=0){
if(t142_h<=0){
REAL t143_l_n_1=t11_h*t142_l_n;
REAL t143_h_1=t11_l_n*t142_l_n;
t143_h=t143_h_1+t143_e;
t143_l_n=t143_l_n_1+t143_e;
}
else if(t142_l_n>=0){
REAL t143_l_n_1_1=t11_l_n*t142_h;
REAL t143_l_n_1_2=t142_l_n*t11_h;
REAL t143_l_n_1=max(t143_l_n_1_1,t143_l_n_1_2);
REAL t143_h_1_1=t11_l_n*t142_l_n;
REAL t143_h_1_2=t11_h*t142_h;
REAL t143_h_1=max(t143_h_1_1,t143_h_1_2);
t143_h=t143_h_1+t143_e;
t143_l_n=t143_l_n_1+t143_e;
}
else{
REAL t143_l_n_1=t11_l_n*t142_h;
REAL t143_h_1=t11_h*t142_h;
t143_h=t143_h_1+t143_e;
t143_l_n=t143_l_n_1+t143_e;
}
}
else{
if(t142_h<=0){
REAL t143_l_n_1=t142_l_n*t11_h;
REAL t142_h_n=-t142_h;
REAL t143_h_1=t142_h_n*t11_l_n;
t143_h=t143_h_1+t143_e;
t143_l_n=t143_l_n_1+t143_e;
}
else if(t142_l_n>=0){
REAL t143_l_n_1=t142_l_n*t11_h;
REAL t143_h_1=t11_h*t142_h;
t143_h=t143_h_1+t143_e;
t143_l_n=t143_l_n_1+t143_e;
}
else{
REAL t11_l=-t11_l_n;
REAL t143_l_n_1=t11_l*t142_l_n;
REAL t143_h_1=t11_h*t142_h;
t143_h=t143_h_1+t143_e;
t143_l_n=t143_l_n_1+t143_e;
}
}
REAL t145_e_1=fabs(t145);
REAL t145_e=eps*t145_e_1;
REAL  t145_h, t145_l_n;
if(t33_h<=0){
if(t45_h<=0){
REAL t33_h_n=-t33_h;
REAL t145_l_n_1=t33_h_n*t45_h;
REAL t145_h_1=t33_l_n*t45_l_n;
t145_h=t145_h_1+t145_e;
t145_l_n=t145_l_n_1+t145_e;
}
else if(t45_l_n>=0){
REAL t145_l_n_1=t33_l_n*t45_h;
REAL t145_h_1=t33_l_n*t45_l_n;
t145_h=t145_h_1+t145_e;
t145_l_n=t145_l_n_1+t145_e;
}
else{
REAL t145_l_n_1=t33_l_n*t45_h;
REAL t33_h_n=-t33_h;
REAL t145_h_1=t33_h_n*t45_l_n;
t145_h=t145_h_1+t145_e;
t145_l_n=t145_l_n_1+t145_e;
}
}
else if(t33_l_n>=0){
if(t45_h<=0){
REAL t145_l_n_1=t33_h*t45_l_n;
REAL t145_h_1=t33_l_n*t45_l_n;
t145_h=t145_h_1+t145_e;
t145_l_n=t145_l_n_1+t145_e;
}
else if(t45_l_n>=0){
REAL t145_l_n_1_1=t33_l_n*t45_h;
REAL t145_l_n_1_2=t45_l_n*t33_h;
REAL t145_l_n_1=max(t145_l_n_1_1,t145_l_n_1_2);
REAL t145_h_1_1=t33_l_n*t45_l_n;
REAL t145_h_1_2=t33_h*t45_h;
REAL t145_h_1=max(t145_h_1_1,t145_h_1_2);
t145_h=t145_h_1+t145_e;
t145_l_n=t145_l_n_1+t145_e;
}
else{
REAL t145_l_n_1=t33_l_n*t45_h;
REAL t145_h_1=t33_h*t45_h;
t145_h=t145_h_1+t145_e;
t145_l_n=t145_l_n_1+t145_e;
}
}
else{
if(t45_h<=0){
REAL t145_l_n_1=t45_l_n*t33_h;
REAL t45_h_n=-t45_h;
REAL t145_h_1=t45_h_n*t33_l_n;
t145_h=t145_h_1+t145_e;
t145_l_n=t145_l_n_1+t145_e;
}
else if(t45_l_n>=0){
REAL t145_l_n_1=t45_l_n*t33_h;
REAL t145_h_1=t33_h*t45_h;
t145_h=t145_h_1+t145_e;
t145_l_n=t145_l_n_1+t145_e;
}
else{
REAL t33_l=-t33_l_n;
REAL t145_l_n_1=t33_l*t45_l_n;
REAL t145_h_1=t33_h*t45_h;
t145_h=t145_h_1+t145_e;
t145_l_n=t145_l_n_1+t145_e;
}
}
REAL t146_e_1=fabs(t146);
REAL t146_e=eps*t146_e_1;
REAL  t146_h, t146_l_n;
if(t55_h<=0){
if(t67_h<=0){
REAL t55_h_n=-t55_h;
REAL t146_l_n_1=t55_h_n*t67_h;
REAL t146_h_1=t55_l_n*t67_l_n;
t146_h=t146_h_1+t146_e;
t146_l_n=t146_l_n_1+t146_e;
}
else if(t67_l_n>=0){
REAL t146_l_n_1=t55_l_n*t67_h;
REAL t146_h_1=t55_l_n*t67_l_n;
t146_h=t146_h_1+t146_e;
t146_l_n=t146_l_n_1+t146_e;
}
else{
REAL t146_l_n_1=t55_l_n*t67_h;
REAL t55_h_n=-t55_h;
REAL t146_h_1=t55_h_n*t67_l_n;
t146_h=t146_h_1+t146_e;
t146_l_n=t146_l_n_1+t146_e;
}
}
else if(t55_l_n>=0){
if(t67_h<=0){
REAL t146_l_n_1=t55_h*t67_l_n;
REAL t146_h_1=t55_l_n*t67_l_n;
t146_h=t146_h_1+t146_e;
t146_l_n=t146_l_n_1+t146_e;
}
else if(t67_l_n>=0){
REAL t146_l_n_1_1=t55_l_n*t67_h;
REAL t146_l_n_1_2=t67_l_n*t55_h;
REAL t146_l_n_1=max(t146_l_n_1_1,t146_l_n_1_2);
REAL t146_h_1_1=t55_l_n*t67_l_n;
REAL t146_h_1_2=t55_h*t67_h;
REAL t146_h_1=max(t146_h_1_1,t146_h_1_2);
t146_h=t146_h_1+t146_e;
t146_l_n=t146_l_n_1+t146_e;
}
else{
REAL t146_l_n_1=t55_l_n*t67_h;
REAL t146_h_1=t55_h*t67_h;
t146_h=t146_h_1+t146_e;
t146_l_n=t146_l_n_1+t146_e;
}
}
else{
if(t67_h<=0){
REAL t146_l_n_1=t67_l_n*t55_h;
REAL t67_h_n=-t67_h;
REAL t146_h_1=t67_h_n*t55_l_n;
t146_h=t146_h_1+t146_e;
t146_l_n=t146_l_n_1+t146_e;
}
else if(t67_l_n>=0){
REAL t146_l_n_1=t67_l_n*t55_h;
REAL t146_h_1=t55_h*t67_h;
t146_h=t146_h_1+t146_e;
t146_l_n=t146_l_n_1+t146_e;
}
else{
REAL t55_l=-t55_l_n;
REAL t146_l_n_1=t55_l*t67_l_n;
REAL t146_h_1=t55_h*t67_h;
t146_h=t146_h_1+t146_e;
t146_l_n=t146_l_n_1+t146_e;
}
}
REAL t148_e_1=fabs(t148);
REAL t148_e=eps*t148_e_1;
REAL t148_h_1=t145_h+t146_l_n;
REAL t148_h=t148_h_1+t148_e;
REAL t148_l_n_1=t146_h+t145_l_n;
REAL t148_l_n=t148_l_n_1+t148_e;
REAL t149_e_1=fabs(t149);
REAL t149_e=eps*t149_e_1;
REAL ax1_n=-ax1;
REAL  t149_h, t149_l_n;
if(ax1<=0){
if(t11_h<=0){
REAL ax1_n=-ax1;
REAL t149_l_n_1=ax1_n*t11_h;
REAL t149_h_1=ax1_n*t11_l_n;
t149_h=t149_h_1+t149_e;
t149_l_n=t149_l_n_1+t149_e;
}
else if(t11_l_n>=0){
REAL t149_l_n_1=ax1_n*t11_h;
REAL t149_h_1=ax1_n*t11_l_n;
t149_h=t149_h_1+t149_e;
t149_l_n=t149_l_n_1+t149_e;
}
else{
REAL t149_l_n_1=ax1_n*t11_h;
REAL ax1_n=-ax1;
REAL t149_h_1=ax1_n*t11_l_n;
t149_h=t149_h_1+t149_e;
t149_l_n=t149_l_n_1+t149_e;
}
}
else if(ax1_n>=0){
if(t11_h<=0){
REAL t149_l_n_1=ax1*t11_l_n;
REAL t149_h_1=ax1_n*t11_l_n;
t149_h=t149_h_1+t149_e;
t149_l_n=t149_l_n_1+t149_e;
}
else if(t11_l_n>=0){
REAL t149_l_n_1_1=ax1_n*t11_h;
REAL t149_l_n_1_2=t11_l_n*ax1;
REAL t149_l_n_1=max(t149_l_n_1_1,t149_l_n_1_2);
REAL t149_h_1_1=ax1_n*t11_l_n;
REAL t149_h_1_2=ax1*t11_h;
REAL t149_h_1=max(t149_h_1_1,t149_h_1_2);
t149_h=t149_h_1+t149_e;
t149_l_n=t149_l_n_1+t149_e;
}
else{
REAL t149_l_n_1=ax1_n*t11_h;
REAL t149_h_1=ax1*t11_h;
t149_h=t149_h_1+t149_e;
t149_l_n=t149_l_n_1+t149_e;
}
}
else{
if(t11_h<=0){
REAL t149_l_n_1=t11_l_n*ax1;
REAL t11_h_n=-t11_h;
REAL t149_h_1=t11_h_n*ax1_n;
t149_h=t149_h_1+t149_e;
t149_l_n=t149_l_n_1+t149_e;
}
else if(t11_l_n>=0){
REAL t149_l_n_1=t11_l_n*ax1;
REAL t149_h_1=ax1*t11_h;
t149_h=t149_h_1+t149_e;
t149_l_n=t149_l_n_1+t149_e;
}
else{
REAL ax1=-ax1_n;
REAL t149_l_n_1=ax1*t11_l_n;
REAL t149_h_1=ax1*t11_h;
t149_h=t149_h_1+t149_e;
t149_l_n=t149_l_n_1+t149_e;
}
}
REAL t150_e_1=fabs(t150);
REAL t150_e=eps*t150_e_1;
REAL az1_n=-az1;
REAL  t150_h, t150_l_n;
if(az1<=0){
if(t19_h<=0){
REAL az1_n=-az1;
REAL t150_l_n_1=az1_n*t19_h;
REAL t150_h_1=az1_n*t19_l_n;
t150_h=t150_h_1+t150_e;
t150_l_n=t150_l_n_1+t150_e;
}
else if(t19_l_n>=0){
REAL t150_l_n_1=az1_n*t19_h;
REAL t150_h_1=az1_n*t19_l_n;
t150_h=t150_h_1+t150_e;
t150_l_n=t150_l_n_1+t150_e;
}
else{
REAL t150_l_n_1=az1_n*t19_h;
REAL az1_n=-az1;
REAL t150_h_1=az1_n*t19_l_n;
t150_h=t150_h_1+t150_e;
t150_l_n=t150_l_n_1+t150_e;
}
}
else if(az1_n>=0){
if(t19_h<=0){
REAL t150_l_n_1=az1*t19_l_n;
REAL t150_h_1=az1_n*t19_l_n;
t150_h=t150_h_1+t150_e;
t150_l_n=t150_l_n_1+t150_e;
}
else if(t19_l_n>=0){
REAL t150_l_n_1_1=az1_n*t19_h;
REAL t150_l_n_1_2=t19_l_n*az1;
REAL t150_l_n_1=max(t150_l_n_1_1,t150_l_n_1_2);
REAL t150_h_1_1=az1_n*t19_l_n;
REAL t150_h_1_2=az1*t19_h;
REAL t150_h_1=max(t150_h_1_1,t150_h_1_2);
t150_h=t150_h_1+t150_e;
t150_l_n=t150_l_n_1+t150_e;
}
else{
REAL t150_l_n_1=az1_n*t19_h;
REAL t150_h_1=az1*t19_h;
t150_h=t150_h_1+t150_e;
t150_l_n=t150_l_n_1+t150_e;
}
}
else{
if(t19_h<=0){
REAL t150_l_n_1=t19_l_n*az1;
REAL t19_h_n=-t19_h;
REAL t150_h_1=t19_h_n*az1_n;
t150_h=t150_h_1+t150_e;
t150_l_n=t150_l_n_1+t150_e;
}
else if(t19_l_n>=0){
REAL t150_l_n_1=t19_l_n*az1;
REAL t150_h_1=az1*t19_h;
t150_h=t150_h_1+t150_e;
t150_l_n=t150_l_n_1+t150_e;
}
else{
REAL az1=-az1_n;
REAL t150_l_n_1=az1*t19_l_n;
REAL t150_h_1=az1*t19_h;
t150_h=t150_h_1+t150_e;
t150_l_n=t150_l_n_1+t150_e;
}
}
REAL t151_e_1=fabs(t151);
REAL t151_e=eps*t151_e_1;
REAL ay1_n=-ay1;
REAL  t151_h, t151_l_n;
if(ay1<=0){
if(t23_h<=0){
REAL ay1_n=-ay1;
REAL t151_l_n_1=ay1_n*t23_h;
REAL t151_h_1=ay1_n*t23_l_n;
t151_h=t151_h_1+t151_e;
t151_l_n=t151_l_n_1+t151_e;
}
else if(t23_l_n>=0){
REAL t151_l_n_1=ay1_n*t23_h;
REAL t151_h_1=ay1_n*t23_l_n;
t151_h=t151_h_1+t151_e;
t151_l_n=t151_l_n_1+t151_e;
}
else{
REAL t151_l_n_1=ay1_n*t23_h;
REAL ay1_n=-ay1;
REAL t151_h_1=ay1_n*t23_l_n;
t151_h=t151_h_1+t151_e;
t151_l_n=t151_l_n_1+t151_e;
}
}
else if(ay1_n>=0){
if(t23_h<=0){
REAL t151_l_n_1=ay1*t23_l_n;
REAL t151_h_1=ay1_n*t23_l_n;
t151_h=t151_h_1+t151_e;
t151_l_n=t151_l_n_1+t151_e;
}
else if(t23_l_n>=0){
REAL t151_l_n_1_1=ay1_n*t23_h;
REAL t151_l_n_1_2=t23_l_n*ay1;
REAL t151_l_n_1=max(t151_l_n_1_1,t151_l_n_1_2);
REAL t151_h_1_1=ay1_n*t23_l_n;
REAL t151_h_1_2=ay1*t23_h;
REAL t151_h_1=max(t151_h_1_1,t151_h_1_2);
t151_h=t151_h_1+t151_e;
t151_l_n=t151_l_n_1+t151_e;
}
else{
REAL t151_l_n_1=ay1_n*t23_h;
REAL t151_h_1=ay1*t23_h;
t151_h=t151_h_1+t151_e;
t151_l_n=t151_l_n_1+t151_e;
}
}
else{
if(t23_h<=0){
REAL t151_l_n_1=t23_l_n*ay1;
REAL t23_h_n=-t23_h;
REAL t151_h_1=t23_h_n*ay1_n;
t151_h=t151_h_1+t151_e;
t151_l_n=t151_l_n_1+t151_e;
}
else if(t23_l_n>=0){
REAL t151_l_n_1=t23_l_n*ay1;
REAL t151_h_1=ay1*t23_h;
t151_h=t151_h_1+t151_e;
t151_l_n=t151_l_n_1+t151_e;
}
else{
REAL ay1=-ay1_n;
REAL t151_l_n_1=ay1*t23_l_n;
REAL t151_h_1=ay1*t23_h;
t151_h=t151_h_1+t151_e;
t151_l_n=t151_l_n_1+t151_e;
}
}
REAL t153_e_1=fabs(t153);
REAL t153_e=eps*t153_e_1;
REAL t153_h_1=t149_h+t150_h;
REAL t153_h=t153_h_1+t153_e;
REAL t153_l_n_1=t149_l_n+t150_l_n;
REAL t153_l_n=t153_l_n_1+t153_e;
REAL t154_e_1=fabs(t154);
REAL t154_e=eps*t154_e_1;
REAL t154_h_1=t153_h+t151_l_n;
REAL t154_h=t154_h_1+t154_e;
REAL t154_l_n_1=t151_h+t153_l_n;
REAL t154_l_n=t154_l_n_1+t154_e;
REAL t155_e_1=fabs(t155);
REAL t155_e=eps*t155_e_1;
REAL  t155_h, t155_l_n;
if(t148_h<=0){
if(t154_h<=0){
REAL t148_h_n=-t148_h;
REAL t155_l_n_1=t148_h_n*t154_h;
REAL t155_h_1=t148_l_n*t154_l_n;
t155_h=t155_h_1+t155_e;
t155_l_n=t155_l_n_1+t155_e;
}
else if(t154_l_n>=0){
REAL t155_l_n_1=t148_l_n*t154_h;
REAL t155_h_1=t148_l_n*t154_l_n;
t155_h=t155_h_1+t155_e;
t155_l_n=t155_l_n_1+t155_e;
}
else{
REAL t155_l_n_1=t148_l_n*t154_h;
REAL t148_h_n=-t148_h;
REAL t155_h_1=t148_h_n*t154_l_n;
t155_h=t155_h_1+t155_e;
t155_l_n=t155_l_n_1+t155_e;
}
}
else if(t148_l_n>=0){
if(t154_h<=0){
REAL t155_l_n_1=t148_h*t154_l_n;
REAL t155_h_1=t148_l_n*t154_l_n;
t155_h=t155_h_1+t155_e;
t155_l_n=t155_l_n_1+t155_e;
}
else if(t154_l_n>=0){
REAL t155_l_n_1_1=t148_l_n*t154_h;
REAL t155_l_n_1_2=t154_l_n*t148_h;
REAL t155_l_n_1=max(t155_l_n_1_1,t155_l_n_1_2);
REAL t155_h_1_1=t148_l_n*t154_l_n;
REAL t155_h_1_2=t148_h*t154_h;
REAL t155_h_1=max(t155_h_1_1,t155_h_1_2);
t155_h=t155_h_1+t155_e;
t155_l_n=t155_l_n_1+t155_e;
}
else{
REAL t155_l_n_1=t148_l_n*t154_h;
REAL t155_h_1=t148_h*t154_h;
t155_h=t155_h_1+t155_e;
t155_l_n=t155_l_n_1+t155_e;
}
}
else{
if(t154_h<=0){
REAL t155_l_n_1=t154_l_n*t148_h;
REAL t154_h_n=-t154_h;
REAL t155_h_1=t154_h_n*t148_l_n;
t155_h=t155_h_1+t155_e;
t155_l_n=t155_l_n_1+t155_e;
}
else if(t154_l_n>=0){
REAL t155_l_n_1=t154_l_n*t148_h;
REAL t155_h_1=t148_h*t154_h;
t155_h=t155_h_1+t155_e;
t155_l_n=t155_l_n_1+t155_e;
}
else{
REAL t148_l=-t148_l_n;
REAL t155_l_n_1=t148_l*t154_l_n;
REAL t155_h_1=t148_h*t154_h;
t155_h=t155_h_1+t155_e;
t155_l_n=t155_l_n_1+t155_e;
}
}
REAL t157_e_1=fabs(t157);
REAL t157_e=eps*t157_e_1;
REAL t157_h_1=t134_h+t143_l_n;
REAL t157_h=t157_h_1+t157_e;
REAL t157_l_n_1=t143_h+t134_l_n;
REAL t157_l_n=t157_l_n_1+t157_e;
REAL t158_e_1=fabs(t158);
REAL t158_e=eps*t158_e_1;
REAL t158_h_1=t157_h+t155_l_n;
REAL t158_h=t158_h_1+t158_e;
REAL t158_l_n_1=t155_h+t157_l_n;
REAL t158_l_n=t158_l_n_1+t158_e;
REAL t159_e_1=fabs(t159);
REAL t159_e=eps*t159_e_1;
REAL  t159_h, t159_l_n;
if(t117_h<=0){
if(t98_h<=0){
REAL t117_h_n=-t117_h;
REAL t159_l_n_1=t117_h_n*t98_h;
REAL t159_h_1=t117_l_n*t98_l_n;
t159_h=t159_h_1+t159_e;
t159_l_n=t159_l_n_1+t159_e;
}
else if(t98_l_n>=0){
REAL t159_l_n_1=t117_l_n*t98_h;
REAL t159_h_1=t117_l_n*t98_l_n;
t159_h=t159_h_1+t159_e;
t159_l_n=t159_l_n_1+t159_e;
}
else{
REAL t159_l_n_1=t117_l_n*t98_h;
REAL t117_h_n=-t117_h;
REAL t159_h_1=t117_h_n*t98_l_n;
t159_h=t159_h_1+t159_e;
t159_l_n=t159_l_n_1+t159_e;
}
}
else if(t117_l_n>=0){
if(t98_h<=0){
REAL t159_l_n_1=t117_h*t98_l_n;
REAL t159_h_1=t117_l_n*t98_l_n;
t159_h=t159_h_1+t159_e;
t159_l_n=t159_l_n_1+t159_e;
}
else if(t98_l_n>=0){
REAL t159_l_n_1_1=t117_l_n*t98_h;
REAL t159_l_n_1_2=t98_l_n*t117_h;
REAL t159_l_n_1=max(t159_l_n_1_1,t159_l_n_1_2);
REAL t159_h_1_1=t117_l_n*t98_l_n;
REAL t159_h_1_2=t117_h*t98_h;
REAL t159_h_1=max(t159_h_1_1,t159_h_1_2);
t159_h=t159_h_1+t159_e;
t159_l_n=t159_l_n_1+t159_e;
}
else{
REAL t159_l_n_1=t117_l_n*t98_h;
REAL t159_h_1=t117_h*t98_h;
t159_h=t159_h_1+t159_e;
t159_l_n=t159_l_n_1+t159_e;
}
}
else{
if(t98_h<=0){
REAL t159_l_n_1=t98_l_n*t117_h;
REAL t98_h_n=-t98_h;
REAL t159_h_1=t98_h_n*t117_l_n;
t159_h=t159_h_1+t159_e;
t159_l_n=t159_l_n_1+t159_e;
}
else if(t98_l_n>=0){
REAL t159_l_n_1=t98_l_n*t117_h;
REAL t159_h_1=t117_h*t98_h;
t159_h=t159_h_1+t159_e;
t159_l_n=t159_l_n_1+t159_e;
}
else{
REAL t117_l=-t117_l_n;
REAL t159_l_n_1=t117_l*t98_l_n;
REAL t159_h_1=t117_h*t98_h;
t159_h=t159_h_1+t159_e;
t159_l_n=t159_l_n_1+t159_e;
}
}
REAL t160_e_1=fabs(t160);
REAL t160_e=eps*t160_e_1;
REAL  t160_h, t160_l_n;
if(t109_h<=0){
if(t115_h<=0){
REAL t109_h_n=-t109_h;
REAL t160_l_n_1=t109_h_n*t115_h;
REAL t160_h_1=t109_l_n*t115_l_n;
t160_h=t160_h_1+t160_e;
t160_l_n=t160_l_n_1+t160_e;
}
else if(t115_l_n>=0){
REAL t160_l_n_1=t109_l_n*t115_h;
REAL t160_h_1=t109_l_n*t115_l_n;
t160_h=t160_h_1+t160_e;
t160_l_n=t160_l_n_1+t160_e;
}
else{
REAL t160_l_n_1=t109_l_n*t115_h;
REAL t109_h_n=-t109_h;
REAL t160_h_1=t109_h_n*t115_l_n;
t160_h=t160_h_1+t160_e;
t160_l_n=t160_l_n_1+t160_e;
}
}
else if(t109_l_n>=0){
if(t115_h<=0){
REAL t160_l_n_1=t109_h*t115_l_n;
REAL t160_h_1=t109_l_n*t115_l_n;
t160_h=t160_h_1+t160_e;
t160_l_n=t160_l_n_1+t160_e;
}
else if(t115_l_n>=0){
REAL t160_l_n_1_1=t109_l_n*t115_h;
REAL t160_l_n_1_2=t115_l_n*t109_h;
REAL t160_l_n_1=max(t160_l_n_1_1,t160_l_n_1_2);
REAL t160_h_1_1=t109_l_n*t115_l_n;
REAL t160_h_1_2=t109_h*t115_h;
REAL t160_h_1=max(t160_h_1_1,t160_h_1_2);
t160_h=t160_h_1+t160_e;
t160_l_n=t160_l_n_1+t160_e;
}
else{
REAL t160_l_n_1=t109_l_n*t115_h;
REAL t160_h_1=t109_h*t115_h;
t160_h=t160_h_1+t160_e;
t160_l_n=t160_l_n_1+t160_e;
}
}
else{
if(t115_h<=0){
REAL t160_l_n_1=t115_l_n*t109_h;
REAL t115_h_n=-t115_h;
REAL t160_h_1=t115_h_n*t109_l_n;
t160_h=t160_h_1+t160_e;
t160_l_n=t160_l_n_1+t160_e;
}
else if(t115_l_n>=0){
REAL t160_l_n_1=t115_l_n*t109_h;
REAL t160_h_1=t109_h*t115_h;
t160_h=t160_h_1+t160_e;
t160_l_n=t160_l_n_1+t160_e;
}
else{
REAL t109_l=-t109_l_n;
REAL t160_l_n_1=t109_l*t115_l_n;
REAL t160_h_1=t109_h*t115_h;
t160_h=t160_h_1+t160_e;
t160_l_n=t160_l_n_1+t160_e;
}
}
REAL t162_e_1=fabs(t162);
REAL t162_e=eps*t162_e_1;
REAL t162_h_1=t159_h+t160_l_n;
REAL t162_h=t162_h_1+t162_e;
REAL t162_l_n_1=t160_h+t159_l_n;
REAL t162_l_n=t162_l_n_1+t162_e;
REAL t163_e_1=fabs(t163);
REAL t163_e=eps*t163_e_1;
REAL  t163_h, t163_l_n;
if(t11_h<=0){
if(t162_h<=0){
REAL t11_h_n=-t11_h;
REAL t163_l_n_1=t11_h_n*t162_h;
REAL t163_h_1=t11_l_n*t162_l_n;
t163_h=t163_h_1+t163_e;
t163_l_n=t163_l_n_1+t163_e;
}
else if(t162_l_n>=0){
REAL t163_l_n_1=t11_l_n*t162_h;
REAL t163_h_1=t11_l_n*t162_l_n;
t163_h=t163_h_1+t163_e;
t163_l_n=t163_l_n_1+t163_e;
}
else{
REAL t163_l_n_1=t11_l_n*t162_h;
REAL t11_h_n=-t11_h;
REAL t163_h_1=t11_h_n*t162_l_n;
t163_h=t163_h_1+t163_e;
t163_l_n=t163_l_n_1+t163_e;
}
}
else if(t11_l_n>=0){
if(t162_h<=0){
REAL t163_l_n_1=t11_h*t162_l_n;
REAL t163_h_1=t11_l_n*t162_l_n;
t163_h=t163_h_1+t163_e;
t163_l_n=t163_l_n_1+t163_e;
}
else if(t162_l_n>=0){
REAL t163_l_n_1_1=t11_l_n*t162_h;
REAL t163_l_n_1_2=t162_l_n*t11_h;
REAL t163_l_n_1=max(t163_l_n_1_1,t163_l_n_1_2);
REAL t163_h_1_1=t11_l_n*t162_l_n;
REAL t163_h_1_2=t11_h*t162_h;
REAL t163_h_1=max(t163_h_1_1,t163_h_1_2);
t163_h=t163_h_1+t163_e;
t163_l_n=t163_l_n_1+t163_e;
}
else{
REAL t163_l_n_1=t11_l_n*t162_h;
REAL t163_h_1=t11_h*t162_h;
t163_h=t163_h_1+t163_e;
t163_l_n=t163_l_n_1+t163_e;
}
}
else{
if(t162_h<=0){
REAL t163_l_n_1=t162_l_n*t11_h;
REAL t162_h_n=-t162_h;
REAL t163_h_1=t162_h_n*t11_l_n;
t163_h=t163_h_1+t163_e;
t163_l_n=t163_l_n_1+t163_e;
}
else if(t162_l_n>=0){
REAL t163_l_n_1=t162_l_n*t11_h;
REAL t163_h_1=t11_h*t162_h;
t163_h=t163_h_1+t163_e;
t163_l_n=t163_l_n_1+t163_e;
}
else{
REAL t11_l=-t11_l_n;
REAL t163_l_n_1=t11_l*t162_l_n;
REAL t163_h_1=t11_h*t162_h;
t163_h=t163_h_1+t163_e;
t163_l_n=t163_l_n_1+t163_e;
}
}
REAL t164_h=t163_l_n;
REAL t164_l_n=t163_h;
REAL t165_e_1=fabs(t165);
REAL t165_e=eps*t165_e_1;
REAL  t165_h, t165_l_n;
if(t105_h<=0){
if(t115_h<=0){
REAL t105_h_n=-t105_h;
REAL t165_l_n_1=t105_h_n*t115_h;
REAL t165_h_1=t105_l_n*t115_l_n;
t165_h=t165_h_1+t165_e;
t165_l_n=t165_l_n_1+t165_e;
}
else if(t115_l_n>=0){
REAL t165_l_n_1=t105_l_n*t115_h;
REAL t165_h_1=t105_l_n*t115_l_n;
t165_h=t165_h_1+t165_e;
t165_l_n=t165_l_n_1+t165_e;
}
else{
REAL t165_l_n_1=t105_l_n*t115_h;
REAL t105_h_n=-t105_h;
REAL t165_h_1=t105_h_n*t115_l_n;
t165_h=t165_h_1+t165_e;
t165_l_n=t165_l_n_1+t165_e;
}
}
else if(t105_l_n>=0){
if(t115_h<=0){
REAL t165_l_n_1=t105_h*t115_l_n;
REAL t165_h_1=t105_l_n*t115_l_n;
t165_h=t165_h_1+t165_e;
t165_l_n=t165_l_n_1+t165_e;
}
else if(t115_l_n>=0){
REAL t165_l_n_1_1=t105_l_n*t115_h;
REAL t165_l_n_1_2=t115_l_n*t105_h;
REAL t165_l_n_1=max(t165_l_n_1_1,t165_l_n_1_2);
REAL t165_h_1_1=t105_l_n*t115_l_n;
REAL t165_h_1_2=t105_h*t115_h;
REAL t165_h_1=max(t165_h_1_1,t165_h_1_2);
t165_h=t165_h_1+t165_e;
t165_l_n=t165_l_n_1+t165_e;
}
else{
REAL t165_l_n_1=t105_l_n*t115_h;
REAL t165_h_1=t105_h*t115_h;
t165_h=t165_h_1+t165_e;
t165_l_n=t165_l_n_1+t165_e;
}
}
else{
if(t115_h<=0){
REAL t165_l_n_1=t115_l_n*t105_h;
REAL t115_h_n=-t115_h;
REAL t165_h_1=t115_h_n*t105_l_n;
t165_h=t165_h_1+t165_e;
t165_l_n=t165_l_n_1+t165_e;
}
else if(t115_l_n>=0){
REAL t165_l_n_1=t115_l_n*t105_h;
REAL t165_h_1=t105_h*t115_h;
t165_h=t165_h_1+t165_e;
t165_l_n=t165_l_n_1+t165_e;
}
else{
REAL t105_l=-t105_l_n;
REAL t165_l_n_1=t105_l*t115_l_n;
REAL t165_h_1=t105_h*t115_h;
t165_h=t165_h_1+t165_e;
t165_l_n=t165_l_n_1+t165_e;
}
}
REAL t166_e_1=fabs(t166);
REAL t166_e=eps*t166_e_1;
REAL  t166_h, t166_l_n;
if(t113_h<=0){
if(t117_h<=0){
REAL t113_h_n=-t113_h;
REAL t166_l_n_1=t113_h_n*t117_h;
REAL t166_h_1=t113_l_n*t117_l_n;
t166_h=t166_h_1+t166_e;
t166_l_n=t166_l_n_1+t166_e;
}
else if(t117_l_n>=0){
REAL t166_l_n_1=t113_l_n*t117_h;
REAL t166_h_1=t113_l_n*t117_l_n;
t166_h=t166_h_1+t166_e;
t166_l_n=t166_l_n_1+t166_e;
}
else{
REAL t166_l_n_1=t113_l_n*t117_h;
REAL t113_h_n=-t113_h;
REAL t166_h_1=t113_h_n*t117_l_n;
t166_h=t166_h_1+t166_e;
t166_l_n=t166_l_n_1+t166_e;
}
}
else if(t113_l_n>=0){
if(t117_h<=0){
REAL t166_l_n_1=t113_h*t117_l_n;
REAL t166_h_1=t113_l_n*t117_l_n;
t166_h=t166_h_1+t166_e;
t166_l_n=t166_l_n_1+t166_e;
}
else if(t117_l_n>=0){
REAL t166_l_n_1_1=t113_l_n*t117_h;
REAL t166_l_n_1_2=t117_l_n*t113_h;
REAL t166_l_n_1=max(t166_l_n_1_1,t166_l_n_1_2);
REAL t166_h_1_1=t113_l_n*t117_l_n;
REAL t166_h_1_2=t113_h*t117_h;
REAL t166_h_1=max(t166_h_1_1,t166_h_1_2);
t166_h=t166_h_1+t166_e;
t166_l_n=t166_l_n_1+t166_e;
}
else{
REAL t166_l_n_1=t113_l_n*t117_h;
REAL t166_h_1=t113_h*t117_h;
t166_h=t166_h_1+t166_e;
t166_l_n=t166_l_n_1+t166_e;
}
}
else{
if(t117_h<=0){
REAL t166_l_n_1=t117_l_n*t113_h;
REAL t117_h_n=-t117_h;
REAL t166_h_1=t117_h_n*t113_l_n;
t166_h=t166_h_1+t166_e;
t166_l_n=t166_l_n_1+t166_e;
}
else if(t117_l_n>=0){
REAL t166_l_n_1=t117_l_n*t113_h;
REAL t166_h_1=t113_h*t117_h;
t166_h=t166_h_1+t166_e;
t166_l_n=t166_l_n_1+t166_e;
}
else{
REAL t113_l=-t113_l_n;
REAL t166_l_n_1=t113_l*t117_l_n;
REAL t166_h_1=t113_h*t117_h;
t166_h=t166_h_1+t166_e;
t166_l_n=t166_l_n_1+t166_e;
}
}
REAL t168_e_1=fabs(t168);
REAL t168_e=eps*t168_e_1;
REAL t168_h_1=t165_h+t166_l_n;
REAL t168_h=t168_h_1+t168_e;
REAL t168_l_n_1=t166_h+t165_l_n;
REAL t168_l_n=t168_l_n_1+t168_e;
REAL t169_e_1=fabs(t169);
REAL t169_e=eps*t169_e_1;
REAL  t169_h, t169_l_n;
if(t168_h<=0){
if(t19_h<=0){
REAL t168_h_n=-t168_h;
REAL t169_l_n_1=t168_h_n*t19_h;
REAL t169_h_1=t168_l_n*t19_l_n;
t169_h=t169_h_1+t169_e;
t169_l_n=t169_l_n_1+t169_e;
}
else if(t19_l_n>=0){
REAL t169_l_n_1=t168_l_n*t19_h;
REAL t169_h_1=t168_l_n*t19_l_n;
t169_h=t169_h_1+t169_e;
t169_l_n=t169_l_n_1+t169_e;
}
else{
REAL t169_l_n_1=t168_l_n*t19_h;
REAL t168_h_n=-t168_h;
REAL t169_h_1=t168_h_n*t19_l_n;
t169_h=t169_h_1+t169_e;
t169_l_n=t169_l_n_1+t169_e;
}
}
else if(t168_l_n>=0){
if(t19_h<=0){
REAL t169_l_n_1=t168_h*t19_l_n;
REAL t169_h_1=t168_l_n*t19_l_n;
t169_h=t169_h_1+t169_e;
t169_l_n=t169_l_n_1+t169_e;
}
else if(t19_l_n>=0){
REAL t169_l_n_1_1=t168_l_n*t19_h;
REAL t169_l_n_1_2=t19_l_n*t168_h;
REAL t169_l_n_1=max(t169_l_n_1_1,t169_l_n_1_2);
REAL t169_h_1_1=t168_l_n*t19_l_n;
REAL t169_h_1_2=t168_h*t19_h;
REAL t169_h_1=max(t169_h_1_1,t169_h_1_2);
t169_h=t169_h_1+t169_e;
t169_l_n=t169_l_n_1+t169_e;
}
else{
REAL t169_l_n_1=t168_l_n*t19_h;
REAL t169_h_1=t168_h*t19_h;
t169_h=t169_h_1+t169_e;
t169_l_n=t169_l_n_1+t169_e;
}
}
else{
if(t19_h<=0){
REAL t169_l_n_1=t19_l_n*t168_h;
REAL t19_h_n=-t19_h;
REAL t169_h_1=t19_h_n*t168_l_n;
t169_h=t169_h_1+t169_e;
t169_l_n=t169_l_n_1+t169_e;
}
else if(t19_l_n>=0){
REAL t169_l_n_1=t19_l_n*t168_h;
REAL t169_h_1=t168_h*t19_h;
t169_h=t169_h_1+t169_e;
t169_l_n=t169_l_n_1+t169_e;
}
else{
REAL t168_l=-t168_l_n;
REAL t169_l_n_1=t168_l*t19_l_n;
REAL t169_h_1=t168_h*t19_h;
t169_h=t169_h_1+t169_e;
t169_l_n=t169_l_n_1+t169_e;
}
}
REAL t171_e_1=fabs(t171);
REAL t171_e=eps*t171_e_1;
REAL  t171_h, t171_l_n;
if(t105_h<=0){
if(t98_h<=0){
REAL t105_h_n=-t105_h;
REAL t171_l_n_1=t105_h_n*t98_h;
REAL t171_h_1=t105_l_n*t98_l_n;
t171_h=t171_h_1+t171_e;
t171_l_n=t171_l_n_1+t171_e;
}
else if(t98_l_n>=0){
REAL t171_l_n_1=t105_l_n*t98_h;
REAL t171_h_1=t105_l_n*t98_l_n;
t171_h=t171_h_1+t171_e;
t171_l_n=t171_l_n_1+t171_e;
}
else{
REAL t171_l_n_1=t105_l_n*t98_h;
REAL t105_h_n=-t105_h;
REAL t171_h_1=t105_h_n*t98_l_n;
t171_h=t171_h_1+t171_e;
t171_l_n=t171_l_n_1+t171_e;
}
}
else if(t105_l_n>=0){
if(t98_h<=0){
REAL t171_l_n_1=t105_h*t98_l_n;
REAL t171_h_1=t105_l_n*t98_l_n;
t171_h=t171_h_1+t171_e;
t171_l_n=t171_l_n_1+t171_e;
}
else if(t98_l_n>=0){
REAL t171_l_n_1_1=t105_l_n*t98_h;
REAL t171_l_n_1_2=t98_l_n*t105_h;
REAL t171_l_n_1=max(t171_l_n_1_1,t171_l_n_1_2);
REAL t171_h_1_1=t105_l_n*t98_l_n;
REAL t171_h_1_2=t105_h*t98_h;
REAL t171_h_1=max(t171_h_1_1,t171_h_1_2);
t171_h=t171_h_1+t171_e;
t171_l_n=t171_l_n_1+t171_e;
}
else{
REAL t171_l_n_1=t105_l_n*t98_h;
REAL t171_h_1=t105_h*t98_h;
t171_h=t171_h_1+t171_e;
t171_l_n=t171_l_n_1+t171_e;
}
}
else{
if(t98_h<=0){
REAL t171_l_n_1=t98_l_n*t105_h;
REAL t98_h_n=-t98_h;
REAL t171_h_1=t98_h_n*t105_l_n;
t171_h=t171_h_1+t171_e;
t171_l_n=t171_l_n_1+t171_e;
}
else if(t98_l_n>=0){
REAL t171_l_n_1=t98_l_n*t105_h;
REAL t171_h_1=t105_h*t98_h;
t171_h=t171_h_1+t171_e;
t171_l_n=t171_l_n_1+t171_e;
}
else{
REAL t105_l=-t105_l_n;
REAL t171_l_n_1=t105_l*t98_l_n;
REAL t171_h_1=t105_h*t98_h;
t171_h=t171_h_1+t171_e;
t171_l_n=t171_l_n_1+t171_e;
}
}
REAL t172_e_1=fabs(t172);
REAL t172_e=eps*t172_e_1;
REAL  t172_h, t172_l_n;
if(t109_h<=0){
if(t113_h<=0){
REAL t109_h_n=-t109_h;
REAL t172_l_n_1=t109_h_n*t113_h;
REAL t172_h_1=t109_l_n*t113_l_n;
t172_h=t172_h_1+t172_e;
t172_l_n=t172_l_n_1+t172_e;
}
else if(t113_l_n>=0){
REAL t172_l_n_1=t109_l_n*t113_h;
REAL t172_h_1=t109_l_n*t113_l_n;
t172_h=t172_h_1+t172_e;
t172_l_n=t172_l_n_1+t172_e;
}
else{
REAL t172_l_n_1=t109_l_n*t113_h;
REAL t109_h_n=-t109_h;
REAL t172_h_1=t109_h_n*t113_l_n;
t172_h=t172_h_1+t172_e;
t172_l_n=t172_l_n_1+t172_e;
}
}
else if(t109_l_n>=0){
if(t113_h<=0){
REAL t172_l_n_1=t109_h*t113_l_n;
REAL t172_h_1=t109_l_n*t113_l_n;
t172_h=t172_h_1+t172_e;
t172_l_n=t172_l_n_1+t172_e;
}
else if(t113_l_n>=0){
REAL t172_l_n_1_1=t109_l_n*t113_h;
REAL t172_l_n_1_2=t113_l_n*t109_h;
REAL t172_l_n_1=max(t172_l_n_1_1,t172_l_n_1_2);
REAL t172_h_1_1=t109_l_n*t113_l_n;
REAL t172_h_1_2=t109_h*t113_h;
REAL t172_h_1=max(t172_h_1_1,t172_h_1_2);
t172_h=t172_h_1+t172_e;
t172_l_n=t172_l_n_1+t172_e;
}
else{
REAL t172_l_n_1=t109_l_n*t113_h;
REAL t172_h_1=t109_h*t113_h;
t172_h=t172_h_1+t172_e;
t172_l_n=t172_l_n_1+t172_e;
}
}
else{
if(t113_h<=0){
REAL t172_l_n_1=t113_l_n*t109_h;
REAL t113_h_n=-t113_h;
REAL t172_h_1=t113_h_n*t109_l_n;
t172_h=t172_h_1+t172_e;
t172_l_n=t172_l_n_1+t172_e;
}
else if(t113_l_n>=0){
REAL t172_l_n_1=t113_l_n*t109_h;
REAL t172_h_1=t109_h*t113_h;
t172_h=t172_h_1+t172_e;
t172_l_n=t172_l_n_1+t172_e;
}
else{
REAL t109_l=-t109_l_n;
REAL t172_l_n_1=t109_l*t113_l_n;
REAL t172_h_1=t109_h*t113_h;
t172_h=t172_h_1+t172_e;
t172_l_n=t172_l_n_1+t172_e;
}
}
REAL t174_e_1=fabs(t174);
REAL t174_e=eps*t174_e_1;
REAL t174_h_1=t171_h+t172_l_n;
REAL t174_h=t174_h_1+t174_e;
REAL t174_l_n_1=t172_h+t171_l_n;
REAL t174_l_n=t174_l_n_1+t174_e;
REAL t175_e_1=fabs(t175);
REAL t175_e=eps*t175_e_1;
REAL t175_h=t175+t175_e;
REAL t175_l_n=t175_e-t175;
REAL t176_e_1=fabs(t176);
REAL t176_e=eps*t176_e_1;
REAL t176_h=t176+t176_e;
REAL t176_l_n=t176_e-t176;
REAL t177_e_1=fabs(t177);
REAL t177_e=eps*t177_e_1;
REAL  t177_h, t177_l_n;
if(t175_h<=0){
if(t176_h<=0){
REAL t175_h_n=-t175_h;
REAL t177_l_n_1=t175_h_n*t176_h;
REAL t177_h_1=t175_l_n*t176_l_n;
t177_h=t177_h_1+t177_e;
t177_l_n=t177_l_n_1+t177_e;
}
else if(t176_l_n>=0){
REAL t177_l_n_1=t175_l_n*t176_h;
REAL t177_h_1=t175_l_n*t176_l_n;
t177_h=t177_h_1+t177_e;
t177_l_n=t177_l_n_1+t177_e;
}
else{
REAL t177_l_n_1=t175_l_n*t176_h;
REAL t175_h_n=-t175_h;
REAL t177_h_1=t175_h_n*t176_l_n;
t177_h=t177_h_1+t177_e;
t177_l_n=t177_l_n_1+t177_e;
}
}
else if(t175_l_n>=0){
if(t176_h<=0){
REAL t177_l_n_1=t175_h*t176_l_n;
REAL t177_h_1=t175_l_n*t176_l_n;
t177_h=t177_h_1+t177_e;
t177_l_n=t177_l_n_1+t177_e;
}
else if(t176_l_n>=0){
REAL t177_l_n_1_1=t175_l_n*t176_h;
REAL t177_l_n_1_2=t176_l_n*t175_h;
REAL t177_l_n_1=max(t177_l_n_1_1,t177_l_n_1_2);
REAL t177_h_1_1=t175_l_n*t176_l_n;
REAL t177_h_1_2=t175_h*t176_h;
REAL t177_h_1=max(t177_h_1_1,t177_h_1_2);
t177_h=t177_h_1+t177_e;
t177_l_n=t177_l_n_1+t177_e;
}
else{
REAL t177_l_n_1=t175_l_n*t176_h;
REAL t177_h_1=t175_h*t176_h;
t177_h=t177_h_1+t177_e;
t177_l_n=t177_l_n_1+t177_e;
}
}
else{
if(t176_h<=0){
REAL t177_l_n_1=t176_l_n*t175_h;
REAL t176_h_n=-t176_h;
REAL t177_h_1=t176_h_n*t175_l_n;
t177_h=t177_h_1+t177_e;
t177_l_n=t177_l_n_1+t177_e;
}
else if(t176_l_n>=0){
REAL t177_l_n_1=t176_l_n*t175_h;
REAL t177_h_1=t175_h*t176_h;
t177_h=t177_h_1+t177_e;
t177_l_n=t177_l_n_1+t177_e;
}
else{
REAL t175_l=-t175_l_n;
REAL t177_l_n_1=t175_l*t176_l_n;
REAL t177_h_1=t175_h*t176_h;
t177_h=t177_h_1+t177_e;
t177_l_n=t177_l_n_1+t177_e;
}
}
REAL t178_e_1=fabs(t178);
REAL t178_e=eps*t178_e_1;
REAL t178_h=t178+t178_e;
REAL t178_l_n=t178_e-t178;
REAL t179_e_1=fabs(t179);
REAL t179_e=eps*t179_e_1;
REAL t179_h=t179+t179_e;
REAL t179_l_n=t179_e-t179;
REAL t180_e_1=fabs(t180);
REAL t180_e=eps*t180_e_1;
REAL  t180_h, t180_l_n;
if(t178_h<=0){
if(t179_h<=0){
REAL t178_h_n=-t178_h;
REAL t180_l_n_1=t178_h_n*t179_h;
REAL t180_h_1=t178_l_n*t179_l_n;
t180_h=t180_h_1+t180_e;
t180_l_n=t180_l_n_1+t180_e;
}
else if(t179_l_n>=0){
REAL t180_l_n_1=t178_l_n*t179_h;
REAL t180_h_1=t178_l_n*t179_l_n;
t180_h=t180_h_1+t180_e;
t180_l_n=t180_l_n_1+t180_e;
}
else{
REAL t180_l_n_1=t178_l_n*t179_h;
REAL t178_h_n=-t178_h;
REAL t180_h_1=t178_h_n*t179_l_n;
t180_h=t180_h_1+t180_e;
t180_l_n=t180_l_n_1+t180_e;
}
}
else if(t178_l_n>=0){
if(t179_h<=0){
REAL t180_l_n_1=t178_h*t179_l_n;
REAL t180_h_1=t178_l_n*t179_l_n;
t180_h=t180_h_1+t180_e;
t180_l_n=t180_l_n_1+t180_e;
}
else if(t179_l_n>=0){
REAL t180_l_n_1_1=t178_l_n*t179_h;
REAL t180_l_n_1_2=t179_l_n*t178_h;
REAL t180_l_n_1=max(t180_l_n_1_1,t180_l_n_1_2);
REAL t180_h_1_1=t178_l_n*t179_l_n;
REAL t180_h_1_2=t178_h*t179_h;
REAL t180_h_1=max(t180_h_1_1,t180_h_1_2);
t180_h=t180_h_1+t180_e;
t180_l_n=t180_l_n_1+t180_e;
}
else{
REAL t180_l_n_1=t178_l_n*t179_h;
REAL t180_h_1=t178_h*t179_h;
t180_h=t180_h_1+t180_e;
t180_l_n=t180_l_n_1+t180_e;
}
}
else{
if(t179_h<=0){
REAL t180_l_n_1=t179_l_n*t178_h;
REAL t179_h_n=-t179_h;
REAL t180_h_1=t179_h_n*t178_l_n;
t180_h=t180_h_1+t180_e;
t180_l_n=t180_l_n_1+t180_e;
}
else if(t179_l_n>=0){
REAL t180_l_n_1=t179_l_n*t178_h;
REAL t180_h_1=t178_h*t179_h;
t180_h=t180_h_1+t180_e;
t180_l_n=t180_l_n_1+t180_e;
}
else{
REAL t178_l=-t178_l_n;
REAL t180_l_n_1=t178_l*t179_l_n;
REAL t180_h_1=t178_h*t179_h;
t180_h=t180_h_1+t180_e;
t180_l_n=t180_l_n_1+t180_e;
}
}
REAL t182_e_1=fabs(t182);
REAL t182_e=eps*t182_e_1;
REAL t182_h_1=t177_h+t180_l_n;
REAL t182_h=t182_h_1+t182_e;
REAL t182_l_n_1=t180_h+t177_l_n;
REAL t182_l_n=t182_l_n_1+t182_e;
REAL t183_e_1=fabs(t183);
REAL t183_e=eps*t183_e_1;
REAL  t183_h, t183_l_n;
if(t174_h<=0){
if(t182_h<=0){
REAL t174_h_n=-t174_h;
REAL t183_l_n_1=t174_h_n*t182_h;
REAL t183_h_1=t174_l_n*t182_l_n;
t183_h=t183_h_1+t183_e;
t183_l_n=t183_l_n_1+t183_e;
}
else if(t182_l_n>=0){
REAL t183_l_n_1=t174_l_n*t182_h;
REAL t183_h_1=t174_l_n*t182_l_n;
t183_h=t183_h_1+t183_e;
t183_l_n=t183_l_n_1+t183_e;
}
else{
REAL t183_l_n_1=t174_l_n*t182_h;
REAL t174_h_n=-t174_h;
REAL t183_h_1=t174_h_n*t182_l_n;
t183_h=t183_h_1+t183_e;
t183_l_n=t183_l_n_1+t183_e;
}
}
else if(t174_l_n>=0){
if(t182_h<=0){
REAL t183_l_n_1=t174_h*t182_l_n;
REAL t183_h_1=t174_l_n*t182_l_n;
t183_h=t183_h_1+t183_e;
t183_l_n=t183_l_n_1+t183_e;
}
else if(t182_l_n>=0){
REAL t183_l_n_1_1=t174_l_n*t182_h;
REAL t183_l_n_1_2=t182_l_n*t174_h;
REAL t183_l_n_1=max(t183_l_n_1_1,t183_l_n_1_2);
REAL t183_h_1_1=t174_l_n*t182_l_n;
REAL t183_h_1_2=t174_h*t182_h;
REAL t183_h_1=max(t183_h_1_1,t183_h_1_2);
t183_h=t183_h_1+t183_e;
t183_l_n=t183_l_n_1+t183_e;
}
else{
REAL t183_l_n_1=t174_l_n*t182_h;
REAL t183_h_1=t174_h*t182_h;
t183_h=t183_h_1+t183_e;
t183_l_n=t183_l_n_1+t183_e;
}
}
else{
if(t182_h<=0){
REAL t183_l_n_1=t182_l_n*t174_h;
REAL t182_h_n=-t182_h;
REAL t183_h_1=t182_h_n*t174_l_n;
t183_h=t183_h_1+t183_e;
t183_l_n=t183_l_n_1+t183_e;
}
else if(t182_l_n>=0){
REAL t183_l_n_1=t182_l_n*t174_h;
REAL t183_h_1=t174_h*t182_h;
t183_h=t183_h_1+t183_e;
t183_l_n=t183_l_n_1+t183_e;
}
else{
REAL t174_l=-t174_l_n;
REAL t183_l_n_1=t174_l*t182_l_n;
REAL t183_h_1=t174_h*t182_h;
t183_h=t183_h_1+t183_e;
t183_l_n=t183_l_n_1+t183_e;
}
}
REAL t185_e_1=fabs(t185);
REAL t185_e=eps*t185_e_1;
REAL t185_h_1=t164_h+t169_l_n;
REAL t185_h=t185_h_1+t185_e;
REAL t185_l_n_1=t169_h+t164_l_n;
REAL t185_l_n=t185_l_n_1+t185_e;
REAL t186_e_1=fabs(t186);
REAL t186_e=eps*t186_e_1;
REAL t186_h_1=t185_h+t183_l_n;
REAL t186_h=t186_h_1+t186_e;
REAL t186_l_n_1=t183_h+t185_l_n;
REAL t186_l_n=t186_l_n_1+t186_e;
REAL t187_e_1=fabs(t187);
REAL t187_e=eps*t187_e_1;
if(t186_l_n*t186_h>=0) inferr=1;
REAL t186_l_n_inv=1/t186_l_n;
REAL t186_h_inv=1/t186_h;
REAL  t187_h, t187_l_n;
if(t158_h<=0){
if(t186_l_n_inv<=0){
REAL t158_h_n=-t158_h;
REAL t187_l_n_1=t158_h_n*t186_l_n_inv;
REAL t187_h_1=t158_l_n*t186_h_inv;
t187_h=t187_h_1+t187_e;
t187_l_n=t187_l_n_1+t187_e;
}
else if(t186_h_inv>=0){
REAL t187_l_n_1=t158_l_n*t186_l_n_inv;
REAL t187_h_1=t158_l_n*t186_h_inv;
t187_h=t187_h_1+t187_e;
t187_l_n=t187_l_n_1+t187_e;
}
else{
REAL t187_l_n_1=t158_l_n*t186_l_n_inv;
REAL t158_h_n=-t158_h;
REAL t187_h_1=t158_h_n*t186_h_inv;
t187_h=t187_h_1+t187_e;
t187_l_n=t187_l_n_1+t187_e;
}
}
else if(t158_l_n>=0){
if(t186_l_n_inv<=0){
REAL t187_l_n_1=t158_h*t186_h_inv;
REAL t187_h_1=t158_l_n*t186_h_inv;
t187_h=t187_h_1+t187_e;
t187_l_n=t187_l_n_1+t187_e;
}
else if(t186_h_inv>=0){
REAL t187_l_n_1_1=t158_l_n*t186_l_n_inv;
REAL t187_l_n_1_2=t186_h_inv*t158_h;
REAL t187_l_n_1=max(t187_l_n_1_1,t187_l_n_1_2);
REAL t187_h_1_1=t158_l_n*t186_h_inv;
REAL t187_h_1_2=t158_h*t186_l_n_inv;
REAL t187_h_1=max(t187_h_1_1,t187_h_1_2);
t187_h=t187_h_1+t187_e;
t187_l_n=t187_l_n_1+t187_e;
}
else{
REAL t187_l_n_1=t158_l_n*t186_l_n_inv;
REAL t187_h_1=t158_h*t186_l_n_inv;
t187_h=t187_h_1+t187_e;
t187_l_n=t187_l_n_1+t187_e;
}
}
else{
if(t186_l_n_inv<=0){
REAL t187_l_n_1=t186_h_inv*t158_h;
REAL t186_l_n_inv_n=-t186_l_n_inv;
REAL t187_h_1=t186_l_n_inv_n*t158_l_n;
t187_h=t187_h_1+t187_e;
t187_l_n=t187_l_n_1+t187_e;
}
else if(t186_h_inv>=0){
REAL t187_l_n_1=t186_h_inv*t158_h;
REAL t187_h_1=t158_h*t186_l_n_inv;
t187_h=t187_h_1+t187_e;
t187_l_n=t187_l_n_1+t187_e;
}
else{
REAL t158_l=-t158_l_n;
REAL t187_l_n_1=t158_l*t186_h_inv;
REAL t187_h_1=t158_h*t186_l_n_inv;
t187_h=t187_h_1+t187_e;
t187_l_n=t187_l_n_1+t187_e;
}
}
fesetround(roundmode);
*result=t187;
if (inferr==0) {
*hbound=t187_h;
*lbound=-t187_l_n;
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
int range = 30;
REAL varlist;
REAL res, err,herr,lerr;
REAL FPStotaltime=0; 
REAL IAtotaltime=0;
REAL FPStime=0; 
REAL IAtime=0;
clock_t before, after;
REAL ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3;
int i,j;
for(j=0; j<sample; j++){
FPStime=0;
IAtime=0;


ax1=randfl(range);
ax2=randfl(range);
ax3=randfl(range);
ay1=randfl(range);
ay2=randfl(range);
ay3=randfl(range);
az1=randfl(range);
az2=randfl(range);
az3=randfl(range);
bx1=randfl(range);
bx2=randfl(range);
bx3=randfl(range);
by1=randfl(range);
by2=randfl(range);
by3=randfl(range);
bz1=randfl(range);
bz2=randfl(range);
bz3=randfl(range);
cx1=randfl(range);
cx2=randfl(range);
cx3=randfl(range);
cy1=randfl(range);
cy2=randfl(range);
cy3=randfl(range);
cz1=randfl(range);
cz2=randfl(range);
cz3=randfl(range);

before = clock();
FPSyn_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&err);
FPSyn_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&err);
FPSyn_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&err);
FPSyn_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&err);
FPSyn_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&err);
FPSyn_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&err);
FPSyn_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&err);
FPSyn_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&err);
FPSyn_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&err);
FPSyn_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&err);
FPSyn_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&err);
FPSyn_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&err);
FPSyn_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&err);
FPSyn_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&err);
FPSyn_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&err);
FPSyn_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&err);
FPSyn_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&err);
FPSyn_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&err);
FPSyn_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&err);
FPSyn_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&err);
FPSyn_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&err);
FPSyn_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&err);
FPSyn_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&err);
FPSyn_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&err);
FPSyn_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&err);
FPSyn_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&err);
FPSyn_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&err);
FPSyn_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&err);
FPSyn_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&err);
FPSyn_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&err);
FPSyn_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&err);
FPSyn_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&err);
FPSyn_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&err);
FPSyn_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&err);
FPSyn_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&err);
FPSyn_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&err);
FPSyn_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&err);
FPSyn_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&err);
FPSyn_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&err);
FPSyn_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&err);
FPSyn_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&err);
FPSyn_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&err);
FPSyn_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&err);
FPSyn_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&err);
FPSyn_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&err);
FPSyn_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&err);
FPSyn_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&err);
FPSyn_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&err);
FPSyn_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&err);
FPSyn_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&err);
FPSyn_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&err);
FPSyn_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&err);
FPSyn_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&err);
FPSyn_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&err);
FPSyn_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&err);
FPSyn_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&err);
FPSyn_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&err);
FPSyn_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&err);
FPSyn_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&err);
FPSyn_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&err);
FPSyn_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&err);
FPSyn_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&err);
FPSyn_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&err);
FPSyn_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&err);
FPSyn_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&err);
FPSyn_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&err);
FPSyn_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&err);
FPSyn_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&err);
FPSyn_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&err);
FPSyn_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&err);
FPSyn_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&err);
FPSyn_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&err);
FPSyn_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&err);
FPSyn_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&err);
FPSyn_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&err);
FPSyn_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&err);
FPSyn_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&err);
FPSyn_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&err);
FPSyn_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&err);
FPSyn_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&err);
FPSyn_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&err);
FPSyn_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&err);
FPSyn_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&err);
FPSyn_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&err);
FPSyn_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&err);
FPSyn_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&err);
FPSyn_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&err);
FPSyn_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&err);
FPSyn_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&err);
FPSyn_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&err);
FPSyn_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&err);
FPSyn_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&err);
FPSyn_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&err);
FPSyn_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&err);
FPSyn_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&err);
FPSyn_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&err);
FPSyn_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&err);
FPSyn_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&err);
FPSyn_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&err);
FPSyn_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&err);
FPSyn_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&err);
FPSyn_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&err);
FPSyn_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&err);
FPSyn_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&err);
FPSyn_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&err);
FPSyn_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&err);
FPSyn_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&err);
FPSyn_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&err);
FPSyn_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&err);
FPSyn_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&err);
FPSyn_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&err);
FPSyn_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&err);
FPSyn_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&err);
FPSyn_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&err);
FPSyn_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&err);
FPSyn_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&err);
FPSyn_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&err);
FPSyn_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&err);
FPSyn_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&err);
FPSyn_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&err);
FPSyn_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&err);
FPSyn_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&err);
FPSyn_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&err);
FPSyn_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&err);
FPSyn_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&err);
FPSyn_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&err);
FPSyn_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&err);
FPSyn_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&err);
FPSyn_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&err);
FPSyn_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&err);
FPSyn_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&err);
FPSyn_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&err);
FPSyn_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&err);
FPSyn_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&err);
FPSyn_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&err);
FPSyn_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&err);
FPSyn_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&err);
FPSyn_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&err);
FPSyn_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&err);
FPSyn_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&err);
FPSyn_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&err);
FPSyn_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&err);
FPSyn_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&err);
FPSyn_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&err);
FPSyn_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&err);
FPSyn_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&err);
FPSyn_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&err);
FPSyn_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&err);
FPSyn_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&err);
FPSyn_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&err);
FPSyn_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&err);
FPSyn_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&err);
FPSyn_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&err);
FPSyn_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&err);
FPSyn_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&err);
FPSyn_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&err);
FPSyn_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&err);
FPSyn_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&err);
FPSyn_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&err);
FPSyn_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&err);
FPSyn_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&err);
FPSyn_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&err);
FPSyn_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&err);
FPSyn_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&err);
FPSyn_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&err);
FPSyn_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&err);
FPSyn_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&err);
FPSyn_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&err);
FPSyn_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&err);
FPSyn_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&err);
FPSyn_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&err);
FPSyn_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&err);
FPSyn_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&err);
FPSyn_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&err);
FPSyn_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&err);
FPSyn_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&err);
FPSyn_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&err);
FPSyn_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&err);
FPSyn_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&err);
FPSyn_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&err);
FPSyn_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&err);
FPSyn_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&err);
FPSyn_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&err);
FPSyn_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&err);
FPSyn_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&err);
FPSyn_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&err);
FPSyn_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&err);
FPSyn_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&err);
FPSyn_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&err);
FPSyn_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&err);
FPSyn_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&err);
FPSyn_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&err);
FPSyn_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&err);
FPSyn_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&err);
FPSyn_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&err);
FPSyn_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&err);
FPSyn_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&err);
FPSyn_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&err);
FPSyn_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&err);
FPSyn_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&err);
FPSyn_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&err);
FPSyn_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&err);
FPSyn_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&err);
FPSyn_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&err);
FPSyn_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&err);
FPSyn_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&err);
FPSyn_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&err);
FPSyn_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&err);
FPSyn_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&err);
FPSyn_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&err);
FPSyn_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&err);
FPSyn_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&err);
FPSyn_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&err);
FPSyn_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&err);
FPSyn_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&err);
FPSyn_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&err);
FPSyn_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&err);
FPSyn_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&err);
FPSyn_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&err);
FPSyn_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&err);
FPSyn_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&err);
FPSyn_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&err);
FPSyn_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&err);
FPSyn_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&err);
FPSyn_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&err);
FPSyn_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&err);
FPSyn_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&err);
FPSyn_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&err);
FPSyn_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&err);
FPSyn_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&err);
FPSyn_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&err);
FPSyn_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&err);
FPSyn_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&err);
FPSyn_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&err);
FPSyn_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&err);
FPSyn_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&err);
FPSyn_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&err);
FPSyn_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&err);
FPSyn_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&err);
FPSyn_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&err);
FPSyn_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&err);
FPSyn_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&err);
FPSyn_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&err);
FPSyn_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&err);
FPSyn_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&err);
FPSyn_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&err);
FPSyn_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&err);
FPSyn_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&err);
FPSyn_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&err);
FPSyn_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&err);
FPSyn_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&err);
FPSyn_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&err);
FPSyn_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&err);
FPSyn_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&err);
FPSyn_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&err);
FPSyn_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&err);
FPSyn_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&err);
FPSyn_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&err);
FPSyn_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&err);
FPSyn_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&err);
FPSyn_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&err);
FPSyn_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&err);
FPSyn_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&err);
FPSyn_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&err);
FPSyn_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&err);
FPSyn_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&err);
FPSyn_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&err);
FPSyn_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&err);
FPSyn_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&err);
FPSyn_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&err);
FPSyn_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&err);
FPSyn_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&err);
FPSyn_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&err);
FPSyn_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&err);
FPSyn_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&err);
FPSyn_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&err);
FPSyn_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&err);
FPSyn_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&err);
FPSyn_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&err);
FPSyn_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&err);
FPSyn_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&err);
FPSyn_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&err);
FPSyn_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&err);
FPSyn_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&err);
FPSyn_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&err);
FPSyn_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&err);
FPSyn_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&err);
FPSyn_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&err);
FPSyn_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&err);
FPSyn_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&err);
FPSyn_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&err);
FPSyn_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&err);
FPSyn_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&err);
FPSyn_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&err);
FPSyn_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&err);
FPSyn_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&err);
FPSyn_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&err);
FPSyn_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&err);
FPSyn_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&err);
FPSyn_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&err);
FPSyn_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&err);
FPSyn_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&err);
FPSyn_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&err);
FPSyn_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&err);
FPSyn_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&err);
FPSyn_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&err);
FPSyn_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&err);
FPSyn_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&err);
FPSyn_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&err);
FPSyn_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&err);
FPSyn_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&err);
FPSyn_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&err);
FPSyn_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&err);
FPSyn_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&err);
FPSyn_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&err);
FPSyn_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&err);
FPSyn_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&err);
FPSyn_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&err);
FPSyn_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&err);
FPSyn_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&err);
FPSyn_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&err);
FPSyn_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&err);
FPSyn_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&err);
FPSyn_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&err);
FPSyn_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&err);
FPSyn_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&err);
FPSyn_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&err);
FPSyn_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&err);
FPSyn_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&err);
FPSyn_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&err);
FPSyn_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&err);
FPSyn_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&err);
FPSyn_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&err);
FPSyn_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&err);
FPSyn_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&err);
FPSyn_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&err);
FPSyn_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&err);
FPSyn_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&err);
FPSyn_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&err);
FPSyn_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&err);
FPSyn_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&err);
FPSyn_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&err);
FPSyn_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&err);
FPSyn_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&err);
FPSyn_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&err);
FPSyn_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&err);
FPSyn_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&err);
FPSyn_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&err);
FPSyn_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&err);
FPSyn_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&err);
FPSyn_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&err);
FPSyn_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&err);
FPSyn_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&err);
FPSyn_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&err);
FPSyn_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&err);
FPSyn_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&err);
FPSyn_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&err);
FPSyn_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&err);
FPSyn_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&err);
FPSyn_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&err);
FPSyn_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&err);
FPSyn_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&err);
FPSyn_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&err);
FPSyn_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&err);
FPSyn_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&err);
FPSyn_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&err);
FPSyn_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&err);
FPSyn_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&err);
FPSyn_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&err);
FPSyn_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&err);
FPSyn_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&err);
FPSyn_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&err);
FPSyn_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&err);
FPSyn_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&err);
FPSyn_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&err);
FPSyn_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&err);
FPSyn_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&err);
FPSyn_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&err);
FPSyn_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&err);
FPSyn_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&err);
FPSyn_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&err);
FPSyn_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&err);
FPSyn_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&err);
FPSyn_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&err);
FPSyn_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&err);
FPSyn_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&err);
FPSyn_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&err);
FPSyn_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&err);
FPSyn_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&err);
FPSyn_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&err);
FPSyn_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&err);
FPSyn_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&err);
FPSyn_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&err);
FPSyn_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&err);
FPSyn_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&err);
FPSyn_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&err);
FPSyn_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&err);
FPSyn_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&err);
FPSyn_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&err);
FPSyn_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&err);
FPSyn_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&err);
FPSyn_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&err);
FPSyn_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&err);
FPSyn_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&err);
FPSyn_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&err);
FPSyn_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&err);
FPSyn_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&err);
FPSyn_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&err);
FPSyn_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&err);
FPSyn_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&err);
FPSyn_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&err);
FPSyn_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&err);
FPSyn_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&err);
FPSyn_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&err);
FPSyn_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&err);
FPSyn_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&err);
FPSyn_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&err);
FPSyn_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&err);
FPSyn_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&err);
FPSyn_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&err);
FPSyn_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&err);
FPSyn_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&err);
FPSyn_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&err);
FPSyn_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&err);
FPSyn_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&err);
FPSyn_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&err);
FPSyn_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&err);
FPSyn_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&err);
FPSyn_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&err);
FPSyn_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&err);
FPSyn_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&err);
FPSyn_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&err);
FPSyn_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&err);
FPSyn_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&err);
FPSyn_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&err);
FPSyn_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&err);
FPSyn_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&err);
FPSyn_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&err);
FPSyn_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&err);
FPSyn_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&err);
FPSyn_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&err);
FPSyn_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&err);
FPSyn_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&err);
FPSyn_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&err);
FPSyn_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&err);
FPSyn_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&err);
FPSyn_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&err);
FPSyn_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&err);
FPSyn_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&err);
FPSyn_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&err);
FPSyn_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&err);
FPSyn_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&err);
FPSyn_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&err);
FPSyn_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&err);
FPSyn_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&err);
FPSyn_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&err);
FPSyn_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&err);
FPSyn_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&err);
FPSyn_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&err);
FPSyn_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&err);
FPSyn_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&err);
FPSyn_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&err);
FPSyn_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&err);
FPSyn_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&err);
FPSyn_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&err);
FPSyn_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&err);
FPSyn_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&err);
FPSyn_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&err);
FPSyn_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&err);
FPSyn_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&err);
FPSyn_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&err);
FPSyn_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&err);
FPSyn_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&err);
FPSyn_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&err);
FPSyn_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&err);
FPSyn_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&err);
FPSyn_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&err);
FPSyn_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&err);
FPSyn_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&err);
FPSyn_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&err);
FPSyn_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&err);
FPSyn_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&err);
FPSyn_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&err);
FPSyn_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&err);
FPSyn_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&err);
FPSyn_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&err);
FPSyn_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&err);
FPSyn_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&err);
FPSyn_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&err);
FPSyn_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&err);
FPSyn_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&err);
FPSyn_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&err);
FPSyn_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&err);
FPSyn_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&err);
FPSyn_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&err);
FPSyn_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&err);
FPSyn_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&err);
FPSyn_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&err);
FPSyn_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&err);
FPSyn_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&err);
FPSyn_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&err);
FPSyn_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&err);
FPSyn_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&err);
FPSyn_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&err);
FPSyn_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&err);
FPSyn_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&err);
FPSyn_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&err);
FPSyn_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&err);
FPSyn_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&err);
FPSyn_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&err);
FPSyn_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&err);
FPSyn_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&err);
FPSyn_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&err);
FPSyn_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&err);
FPSyn_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&err);
FPSyn_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&err);
FPSyn_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&err);
FPSyn_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&err);
FPSyn_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&err);
FPSyn_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&err);
FPSyn_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&err);
FPSyn_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&err);
FPSyn_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&err);
FPSyn_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&err);
FPSyn_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&err);
FPSyn_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&err);
FPSyn_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&err);
FPSyn_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&err);
FPSyn_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&err);
FPSyn_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&err);
FPSyn_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&err);
FPSyn_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&err);
FPSyn_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&err);
FPSyn_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&err);
FPSyn_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&err);
FPSyn_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&err);
FPSyn_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&err);
FPSyn_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&err);
FPSyn_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&err);
FPSyn_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&err);
FPSyn_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&err);
FPSyn_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&err);
FPSyn_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&err);
FPSyn_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&err);
FPSyn_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&err);
FPSyn_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&err);
FPSyn_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&err);
FPSyn_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&err);
FPSyn_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&err);
FPSyn_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&err);
FPSyn_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&err);
FPSyn_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&err);
FPSyn_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&err);
FPSyn_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&err);
FPSyn_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&err);
FPSyn_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&err);
FPSyn_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&err);
FPSyn_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&err);
FPSyn_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&err);
FPSyn_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&err);
FPSyn_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&err);
FPSyn_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&err);
FPSyn_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&err);
FPSyn_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&err);
FPSyn_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&err);
FPSyn_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&err);
FPSyn_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&err);
FPSyn_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&err);
FPSyn_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&err);
FPSyn_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&err);
FPSyn_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&err);
FPSyn_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&err);
FPSyn_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&err);
FPSyn_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&err);
FPSyn_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&err);
FPSyn_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&err);
FPSyn_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&err);
FPSyn_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&err);
FPSyn_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&err);
FPSyn_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&err);
FPSyn_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&err);
FPSyn_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&err);
FPSyn_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&err);
FPSyn_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&err);
FPSyn_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&err);
FPSyn_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&err);
FPSyn_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&err);
FPSyn_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&err);
FPSyn_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&err);
FPSyn_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&err);
FPSyn_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&err);
FPSyn_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&err);
FPSyn_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&err);
FPSyn_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&err);
FPSyn_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&err);
FPSyn_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&err);
FPSyn_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&err);
FPSyn_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&err);
FPSyn_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&err);
FPSyn_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&err);
FPSyn_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&err);
FPSyn_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&err);
FPSyn_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&err);
FPSyn_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&err);
FPSyn_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&err);
FPSyn_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&err);
FPSyn_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&err);
FPSyn_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&err);
FPSyn_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&err);
FPSyn_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&err);
FPSyn_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&err);
FPSyn_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&err);
FPSyn_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&err);
FPSyn_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&err);
FPSyn_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&err);
FPSyn_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&err);
FPSyn_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&err);
FPSyn_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&err);
FPSyn_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&err);
FPSyn_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&err);
FPSyn_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&err);
FPSyn_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&err);
FPSyn_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&err);
FPSyn_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&err);
FPSyn_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&err);
FPSyn_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&err);
FPSyn_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&err);
FPSyn_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&err);
FPSyn_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&err);
FPSyn_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&err);
FPSyn_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&err);
FPSyn_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&err);
FPSyn_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&err);
FPSyn_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&err);
FPSyn_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&err);
FPSyn_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&err);
FPSyn_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&err);
FPSyn_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&err);
FPSyn_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&err);
FPSyn_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&err);
FPSyn_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&err);
FPSyn_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&err);
FPSyn_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&err);
FPSyn_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&err);
FPSyn_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&err);
FPSyn_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&err);
FPSyn_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&err);
FPSyn_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&err);
FPSyn_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&err);
FPSyn_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&err);
FPSyn_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&err);
FPSyn_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&err);
FPSyn_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&err);
FPSyn_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&err);
FPSyn_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&err);
FPSyn_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&err);
FPSyn_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&err);
FPSyn_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&err);
FPSyn_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&err);
FPSyn_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&err);
FPSyn_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&err);
FPSyn_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&err);
FPSyn_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&err);
FPSyn_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&err);
FPSyn_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&err);
FPSyn_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&err);
FPSyn_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&err);
FPSyn_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&err);
FPSyn_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&err);
FPSyn_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&err);
FPSyn_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&err);
FPSyn_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&err);
FPSyn_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&err);
FPSyn_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&err);
FPSyn_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&err);
FPSyn_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&err);
FPSyn_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&err);
FPSyn_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&err);
FPSyn_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&err);
FPSyn_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&err);
FPSyn_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&err);
FPSyn_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&err);
FPSyn_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&err);
FPSyn_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&err);
FPSyn_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&err);
FPSyn_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&err);
FPSyn_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&err);
FPSyn_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&err);
FPSyn_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&err);
FPSyn_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&err);
FPSyn_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&err);
FPSyn_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&err);
FPSyn_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&err);
FPSyn_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&err);
FPSyn_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&err);
FPSyn_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&err);
FPSyn_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&err);
FPSyn_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&err);
FPSyn_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&err);
FPSyn_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&err);
FPSyn_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&err);
FPSyn_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&err);
FPSyn_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&err);
FPSyn_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&err);
FPSyn_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&err);
FPSyn_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&err);
FPSyn_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&err);
FPSyn_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&err);
FPSyn_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&err);
FPSyn_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&err);
FPSyn_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&err);
FPSyn_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&err);
FPSyn_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&err);
FPSyn_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&err);
FPSyn_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&err);
FPSyn_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&err);
FPSyn_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&err);
FPSyn_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&err);
FPSyn_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&err);
FPSyn_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&err);
FPSyn_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&err);
FPSyn_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&err);
FPSyn_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&err);
FPSyn_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&err);
FPSyn_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&err);
FPSyn_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&err);
FPSyn_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&err);
FPSyn_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&err);
FPSyn_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&err);
FPSyn_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&err);
FPSyn_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&err);
FPSyn_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&err);
FPSyn_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&err);
FPSyn_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&err);
FPSyn_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&err);
FPSyn_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&err);
FPSyn_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&err);
FPSyn_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&err);
FPSyn_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&err);
FPSyn_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&err);
FPSyn_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&err);
FPSyn_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&err);
FPSyn_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&err);
FPSyn_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&err);
FPSyn_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&err);
FPSyn_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&err);
FPSyn_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&err);
FPSyn_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&err);
FPSyn_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&err);
FPSyn_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&err);
FPSyn_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&err);
FPSyn_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&err);
FPSyn_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&err);
FPSyn_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&err);
FPSyn_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&err);
FPSyn_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&err);
FPSyn_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&err);
FPSyn_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&err);
FPSyn_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&err);
FPSyn_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&err);
FPSyn_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&err);
FPSyn_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&err);
FPSyn_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&err);
FPSyn_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&err);
FPSyn_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&err);
FPSyn_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&err);
FPSyn_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&err);
FPSyn_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&err);
FPSyn_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&err);
FPSyn_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&err);
FPSyn_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&err);
FPSyn_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&err);
FPSyn_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&err);
FPSyn_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&err);
FPSyn_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&err);
FPSyn_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&err);
FPSyn_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&err);
FPSyn_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&err);
FPSyn_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&err);
FPSyn_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&err);
FPSyn_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&err);
FPSyn_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&err);
FPSyn_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&err);
FPSyn_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&err);
FPSyn_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&err);
FPSyn_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&err);
FPSyn_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&err);
FPSyn_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&err);
FPSyn_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&err);
FPSyn_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&err);
FPSyn_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&err);
FPSyn_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&err);
FPSyn_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&err);
FPSyn_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&err);
FPSyn_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&err);
FPSyn_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&err);
FPSyn_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&err);
FPSyn_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&err);
FPSyn_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&err);
FPSyn_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&err);
FPSyn_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&err);
FPSyn_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&err);
FPSyn_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&err);
FPSyn_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&err);
FPSyn_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&err);
FPSyn_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&err);
FPSyn_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&err);
FPSyn_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&err);
FPSyn_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&err);
FPSyn_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&err);
FPSyn_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&err);
FPSyn_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&err);
FPSyn_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&err);
FPSyn_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&err);
FPSyn_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&err);
FPSyn_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&err);
FPSyn_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&err);
FPSyn_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&err);
FPSyn_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&err);
FPSyn_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&err);
FPSyn_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&err);
FPSyn_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&err);
FPSyn_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&err);
FPSyn_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&err);
FPSyn_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&err);
FPSyn_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&err);
FPSyn_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&err);
FPSyn_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&err);
FPSyn_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&err);
FPSyn_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&err);
FPSyn_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&err);
FPSyn_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&err);
FPSyn_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&err);
FPSyn_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&err);
FPSyn_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&err);
FPSyn_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&err);
FPSyn_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&err);
FPSyn_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&err);
FPSyn_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&err);
FPSyn_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&err);
FPSyn_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&err);
FPSyn_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&err);
FPSyn_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&err);
FPSyn_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&err);
FPSyn_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&err);
FPSyn_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&err);
FPSyn_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&err);
FPSyn_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&err);
FPSyn_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&err);
FPSyn_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&err);
FPSyn_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&err);
FPSyn_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&err);
FPSyn_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&err);
FPSyn_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&err);
FPSyn_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&err);
FPSyn_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&err);
FPSyn_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&err);
FPSyn_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&err);
FPSyn_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&err);
FPSyn_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&err);
FPSyn_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&err);
FPSyn_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&err);
FPSyn_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&err);
FPSyn_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&err);
FPSyn_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&err);
FPSyn_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&err);
FPSyn_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&err);
FPSyn_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&err);
FPSyn_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&err);
FPSyn_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&err);
FPSyn_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&err);
FPSyn_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&err);
FPSyn_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&err);
FPSyn_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&err);
FPSyn_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&err);
FPSyn_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&err);
FPSyn_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&err);
FPSyn_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&err);
FPSyn_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&err);
FPSyn_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&err);
FPSyn_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&err);
FPSyn_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&err);
FPSyn_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&err);
FPSyn_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&err);
FPSyn_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&err);
FPSyn_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&err);
FPSyn_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&err);
FPSyn_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&err);
FPSyn_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&err);
FPSyn_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&err);
FPSyn_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&err);
FPSyn_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&err);
FPSyn_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&err);
FPSyn_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&err);
FPSyn_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&err);
FPSyn_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&err);
FPSyn_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&err);
FPSyn_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&err);
FPSyn_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&err);
FPSyn_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&err);
FPSyn_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&err);
FPSyn_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&err);
FPSyn_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&err);
FPSyn_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&err);
FPSyn_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&err);
FPSyn_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&err);
FPSyn_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&err);
FPSyn_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&err);
FPSyn_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&err);
FPSyn_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&err);
FPSyn_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&err);
FPSyn_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&err);
FPSyn_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&err);
FPSyn_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&err);
FPSyn_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&err);
FPSyn_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&err);
FPSyn_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&err);
FPSyn_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&err);
FPSyn_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&err);
FPSyn_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&err);
FPSyn_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&err);
FPSyn_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&err);
FPSyn_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&err);
FPSyn_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&err);
FPSyn_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&err);
FPSyn_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&err);
FPSyn_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&err);
FPSyn_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&err);
FPSyn_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&err);
FPSyn_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&err);
FPSyn_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&err);
FPSyn_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&err);
FPSyn_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&err);
FPSyn_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&err);
FPSyn_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&err);
FPSyn_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&err);
FPSyn_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&err);
FPSyn_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&err);
FPSyn_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&err);
FPSyn_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&err);
FPSyn_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&err);
FPSyn_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&err);
FPSyn_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&err);
FPSyn_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&err);
FPSyn_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&err);
FPSyn_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&err);
FPSyn_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&err);
FPSyn_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&err);
FPSyn_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&err);
FPSyn_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&err);
FPSyn_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&err);
FPSyn_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&err);
FPSyn_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&err);
FPSyn_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&err);
FPSyn_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&err);
FPSyn_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&err);
FPSyn_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&err);
FPSyn_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&err);
FPSyn_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&err);
FPSyn_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&err);
FPSyn_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&err);
FPSyn_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&err);
FPSyn_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&err);
FPSyn_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&err);
FPSyn_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&err);
FPSyn_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&err);
FPSyn_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&err);
FPSyn_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&err);
FPSyn_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&err);
FPSyn_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&err);
FPSyn_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&err);
FPSyn_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&err);
FPSyn_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&err);
FPSyn_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&err);
FPSyn_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&err);
FPSyn_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&err);
FPSyn_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&err);
FPSyn_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&err);
FPSyn_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&err);
FPSyn_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&err);
FPSyn_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&err);
FPSyn_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&err);
FPSyn_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&err);
FPSyn_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&err);
FPSyn_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&err);
FPSyn_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&err);
FPSyn_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&err);
FPSyn_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&err);
FPSyn_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&err);
FPSyn_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&err);
FPSyn_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&err);
FPSyn_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&err);
FPSyn_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&err);
FPSyn_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&err);
FPSyn_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&err);
FPSyn_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&err);
FPSyn_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&err);
FPSyn_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&err);
FPSyn_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&err);
FPSyn_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&err);
after = clock();
FPStime = ((double)(after - before))/CLOCKS_PER_SEC;
FPStotaltime = FPStotaltime + FPStime;


before = clock();
IA_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&herr,&lerr);
IA_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&herr,&lerr);
IA_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&herr,&lerr);
IA_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&herr,&lerr);
IA_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&herr,&lerr);
IA_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&herr,&lerr);
IA_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&herr,&lerr);
IA_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&herr,&lerr);
IA_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&herr,&lerr);
IA_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&herr,&lerr);
IA_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&herr,&lerr);
IA_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&herr,&lerr);
IA_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&herr,&lerr);
IA_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&herr,&lerr);
IA_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&herr,&lerr);
IA_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&herr,&lerr);
IA_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&herr,&lerr);
IA_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&herr,&lerr);
IA_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&herr,&lerr);
IA_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&herr,&lerr);
IA_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&herr,&lerr);
IA_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&herr,&lerr);
IA_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&herr,&lerr);
IA_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&herr,&lerr);
IA_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&herr,&lerr);
IA_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&herr,&lerr);
IA_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&herr,&lerr);
IA_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&herr,&lerr);
IA_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&herr,&lerr);
IA_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&herr,&lerr);
IA_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&herr,&lerr);
IA_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&herr,&lerr);
IA_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&herr,&lerr);
IA_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&herr,&lerr);
IA_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&herr,&lerr);
IA_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&herr,&lerr);
IA_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&herr,&lerr);
IA_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&herr,&lerr);
IA_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&herr,&lerr);
IA_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&herr,&lerr);
IA_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&herr,&lerr);
IA_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&herr,&lerr);
IA_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&herr,&lerr);
IA_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&herr,&lerr);
IA_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&herr,&lerr);
IA_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&herr,&lerr);
IA_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&herr,&lerr);
IA_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&herr,&lerr);
IA_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&herr,&lerr);
IA_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&herr,&lerr);
IA_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&herr,&lerr);
IA_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&herr,&lerr);
IA_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&herr,&lerr);
IA_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&herr,&lerr);
IA_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&herr,&lerr);
IA_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&herr,&lerr);
IA_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&herr,&lerr);
IA_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&herr,&lerr);
IA_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&herr,&lerr);
IA_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&herr,&lerr);
IA_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&herr,&lerr);
IA_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&herr,&lerr);
IA_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&herr,&lerr);
IA_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&herr,&lerr);
IA_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&herr,&lerr);
IA_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&herr,&lerr);
IA_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&herr,&lerr);
IA_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&herr,&lerr);
IA_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&herr,&lerr);
IA_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&herr,&lerr);
IA_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&herr,&lerr);
IA_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&herr,&lerr);
IA_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&herr,&lerr);
IA_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&herr,&lerr);
IA_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&herr,&lerr);
IA_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&herr,&lerr);
IA_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&herr,&lerr);
IA_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&herr,&lerr);
IA_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&herr,&lerr);
IA_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&herr,&lerr);
IA_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&herr,&lerr);
IA_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&herr,&lerr);
IA_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&herr,&lerr);
IA_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&herr,&lerr);
IA_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&herr,&lerr);
IA_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&herr,&lerr);
IA_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&herr,&lerr);
IA_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&herr,&lerr);
IA_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&herr,&lerr);
IA_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&herr,&lerr);
IA_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&herr,&lerr);
IA_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&herr,&lerr);
IA_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&herr,&lerr);
IA_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&herr,&lerr);
IA_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&herr,&lerr);
IA_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&herr,&lerr);
IA_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&herr,&lerr);
IA_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&herr,&lerr);
IA_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&herr,&lerr);
IA_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&herr,&lerr);
IA_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&herr,&lerr);
IA_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&herr,&lerr);
IA_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&herr,&lerr);
IA_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&herr,&lerr);
IA_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&herr,&lerr);
IA_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&herr,&lerr);
IA_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&herr,&lerr);
IA_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&herr,&lerr);
IA_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&herr,&lerr);
IA_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&herr,&lerr);
IA_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&herr,&lerr);
IA_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&herr,&lerr);
IA_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&herr,&lerr);
IA_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&herr,&lerr);
IA_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&herr,&lerr);
IA_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&herr,&lerr);
IA_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&herr,&lerr);
IA_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&herr,&lerr);
IA_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&herr,&lerr);
IA_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&herr,&lerr);
IA_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&herr,&lerr);
IA_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&herr,&lerr);
IA_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&herr,&lerr);
IA_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&herr,&lerr);
IA_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&herr,&lerr);
IA_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&herr,&lerr);
IA_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&herr,&lerr);
IA_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&herr,&lerr);
IA_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&herr,&lerr);
IA_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&herr,&lerr);
IA_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&herr,&lerr);
IA_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&herr,&lerr);
IA_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&herr,&lerr);
IA_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&herr,&lerr);
IA_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&herr,&lerr);
IA_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&herr,&lerr);
IA_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&herr,&lerr);
IA_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&herr,&lerr);
IA_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&herr,&lerr);
IA_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&herr,&lerr);
IA_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&herr,&lerr);
IA_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&herr,&lerr);
IA_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&herr,&lerr);
IA_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&herr,&lerr);
IA_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&herr,&lerr);
IA_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&herr,&lerr);
IA_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&herr,&lerr);
IA_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&herr,&lerr);
IA_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&herr,&lerr);
IA_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&herr,&lerr);
IA_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&herr,&lerr);
IA_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&herr,&lerr);
IA_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&herr,&lerr);
IA_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&herr,&lerr);
IA_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&herr,&lerr);
IA_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&herr,&lerr);
IA_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&herr,&lerr);
IA_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&herr,&lerr);
IA_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&herr,&lerr);
IA_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&herr,&lerr);
IA_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&herr,&lerr);
IA_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&herr,&lerr);
IA_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&herr,&lerr);
IA_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&herr,&lerr);
IA_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&herr,&lerr);
IA_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&herr,&lerr);
IA_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&herr,&lerr);
IA_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&herr,&lerr);
IA_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&herr,&lerr);
IA_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&herr,&lerr);
IA_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&herr,&lerr);
IA_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&herr,&lerr);
IA_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&herr,&lerr);
IA_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&herr,&lerr);
IA_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&herr,&lerr);
IA_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&herr,&lerr);
IA_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&herr,&lerr);
IA_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&herr,&lerr);
IA_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&herr,&lerr);
IA_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&herr,&lerr);
IA_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&herr,&lerr);
IA_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&herr,&lerr);
IA_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&herr,&lerr);
IA_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&herr,&lerr);
IA_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&herr,&lerr);
IA_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&herr,&lerr);
IA_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&herr,&lerr);
IA_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&herr,&lerr);
IA_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&herr,&lerr);
IA_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&herr,&lerr);
IA_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&herr,&lerr);
IA_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&herr,&lerr);
IA_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&herr,&lerr);
IA_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&herr,&lerr);
IA_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&herr,&lerr);
IA_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&herr,&lerr);
IA_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&herr,&lerr);
IA_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&herr,&lerr);
IA_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&herr,&lerr);
IA_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&herr,&lerr);
IA_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&herr,&lerr);
IA_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&herr,&lerr);
IA_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&herr,&lerr);
IA_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&herr,&lerr);
IA_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&herr,&lerr);
IA_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&herr,&lerr);
IA_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&herr,&lerr);
IA_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&herr,&lerr);
IA_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&herr,&lerr);
IA_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&herr,&lerr);
IA_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&herr,&lerr);
IA_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&herr,&lerr);
IA_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&herr,&lerr);
IA_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&herr,&lerr);
IA_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&herr,&lerr);
IA_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&herr,&lerr);
IA_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&herr,&lerr);
IA_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&herr,&lerr);
IA_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&herr,&lerr);
IA_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&herr,&lerr);
IA_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&herr,&lerr);
IA_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&herr,&lerr);
IA_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&herr,&lerr);
IA_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&herr,&lerr);
IA_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&herr,&lerr);
IA_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&herr,&lerr);
IA_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&herr,&lerr);
IA_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&herr,&lerr);
IA_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&herr,&lerr);
IA_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&herr,&lerr);
IA_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&herr,&lerr);
IA_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&herr,&lerr);
IA_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&herr,&lerr);
IA_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&herr,&lerr);
IA_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&herr,&lerr);
IA_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&herr,&lerr);
IA_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&herr,&lerr);
IA_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&herr,&lerr);
IA_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&herr,&lerr);
IA_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&herr,&lerr);
IA_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&herr,&lerr);
IA_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&herr,&lerr);
IA_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&herr,&lerr);
IA_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&herr,&lerr);
IA_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&herr,&lerr);
IA_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&herr,&lerr);
IA_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&herr,&lerr);
IA_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&herr,&lerr);
IA_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&herr,&lerr);
IA_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&herr,&lerr);
IA_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&herr,&lerr);
IA_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&herr,&lerr);
IA_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&herr,&lerr);
IA_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&herr,&lerr);
IA_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&herr,&lerr);
IA_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&herr,&lerr);
IA_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&herr,&lerr);
IA_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&herr,&lerr);
IA_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&herr,&lerr);
IA_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&herr,&lerr);
IA_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&herr,&lerr);
IA_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&herr,&lerr);
IA_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&herr,&lerr);
IA_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&herr,&lerr);
IA_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&herr,&lerr);
IA_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&herr,&lerr);
IA_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&herr,&lerr);
IA_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&herr,&lerr);
IA_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&herr,&lerr);
IA_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&herr,&lerr);
IA_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&herr,&lerr);
IA_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&herr,&lerr);
IA_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&herr,&lerr);
IA_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&herr,&lerr);
IA_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&herr,&lerr);
IA_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&herr,&lerr);
IA_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&herr,&lerr);
IA_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&herr,&lerr);
IA_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&herr,&lerr);
IA_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&herr,&lerr);
IA_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&herr,&lerr);
IA_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&herr,&lerr);
IA_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&herr,&lerr);
IA_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&herr,&lerr);
IA_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&herr,&lerr);
IA_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&herr,&lerr);
IA_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&herr,&lerr);
IA_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&herr,&lerr);
IA_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&herr,&lerr);
IA_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&herr,&lerr);
IA_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&herr,&lerr);
IA_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&herr,&lerr);
IA_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&herr,&lerr);
IA_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&herr,&lerr);
IA_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&herr,&lerr);
IA_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&herr,&lerr);
IA_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&herr,&lerr);
IA_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&herr,&lerr);
IA_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&herr,&lerr);
IA_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&herr,&lerr);
IA_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&herr,&lerr);
IA_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&herr,&lerr);
IA_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&herr,&lerr);
IA_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&herr,&lerr);
IA_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&herr,&lerr);
IA_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&herr,&lerr);
IA_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&herr,&lerr);
IA_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&herr,&lerr);
IA_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&herr,&lerr);
IA_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&herr,&lerr);
IA_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&herr,&lerr);
IA_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&herr,&lerr);
IA_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&herr,&lerr);
IA_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&herr,&lerr);
IA_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&herr,&lerr);
IA_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&herr,&lerr);
IA_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&herr,&lerr);
IA_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&herr,&lerr);
IA_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&herr,&lerr);
IA_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&herr,&lerr);
IA_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&herr,&lerr);
IA_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&herr,&lerr);
IA_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&herr,&lerr);
IA_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&herr,&lerr);
IA_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&herr,&lerr);
IA_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&herr,&lerr);
IA_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&herr,&lerr);
IA_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&herr,&lerr);
IA_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&herr,&lerr);
IA_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&herr,&lerr);
IA_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&herr,&lerr);
IA_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&herr,&lerr);
IA_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&herr,&lerr);
IA_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&herr,&lerr);
IA_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&herr,&lerr);
IA_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&herr,&lerr);
IA_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&herr,&lerr);
IA_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&herr,&lerr);
IA_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&herr,&lerr);
IA_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&herr,&lerr);
IA_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&herr,&lerr);
IA_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&herr,&lerr);
IA_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&herr,&lerr);
IA_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&herr,&lerr);
IA_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&herr,&lerr);
IA_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&herr,&lerr);
IA_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&herr,&lerr);
IA_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&herr,&lerr);
IA_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&herr,&lerr);
IA_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&herr,&lerr);
IA_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&herr,&lerr);
IA_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&herr,&lerr);
IA_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&herr,&lerr);
IA_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&herr,&lerr);
IA_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&herr,&lerr);
IA_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&herr,&lerr);
IA_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&herr,&lerr);
IA_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&herr,&lerr);
IA_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&herr,&lerr);
IA_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&herr,&lerr);
IA_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&herr,&lerr);
IA_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&herr,&lerr);
IA_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&herr,&lerr);
IA_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&herr,&lerr);
IA_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&herr,&lerr);
IA_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&herr,&lerr);
IA_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&herr,&lerr);
IA_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&herr,&lerr);
IA_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&herr,&lerr);
IA_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&herr,&lerr);
IA_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&herr,&lerr);
IA_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&herr,&lerr);
IA_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&herr,&lerr);
IA_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&herr,&lerr);
IA_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&herr,&lerr);
IA_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&herr,&lerr);
IA_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&herr,&lerr);
IA_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&herr,&lerr);
IA_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&herr,&lerr);
IA_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&herr,&lerr);
IA_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&herr,&lerr);
IA_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&herr,&lerr);
IA_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&herr,&lerr);
IA_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&herr,&lerr);
IA_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&herr,&lerr);
IA_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&herr,&lerr);
IA_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&herr,&lerr);
IA_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&herr,&lerr);
IA_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&herr,&lerr);
IA_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&herr,&lerr);
IA_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&herr,&lerr);
IA_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&herr,&lerr);
IA_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&herr,&lerr);
IA_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&herr,&lerr);
IA_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&herr,&lerr);
IA_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&herr,&lerr);
IA_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&herr,&lerr);
IA_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&herr,&lerr);
IA_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&herr,&lerr);
IA_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&herr,&lerr);
IA_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&herr,&lerr);
IA_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&herr,&lerr);
IA_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&herr,&lerr);
IA_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&herr,&lerr);
IA_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&herr,&lerr);
IA_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&herr,&lerr);
IA_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&herr,&lerr);
IA_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&herr,&lerr);
IA_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&herr,&lerr);
IA_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&herr,&lerr);
IA_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&herr,&lerr);
IA_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&herr,&lerr);
IA_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&herr,&lerr);
IA_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&herr,&lerr);
IA_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&herr,&lerr);
IA_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&herr,&lerr);
IA_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&herr,&lerr);
IA_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&herr,&lerr);
IA_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&herr,&lerr);
IA_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&herr,&lerr);
IA_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&herr,&lerr);
IA_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&herr,&lerr);
IA_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&herr,&lerr);
IA_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&herr,&lerr);
IA_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&herr,&lerr);
IA_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&herr,&lerr);
IA_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&herr,&lerr);
IA_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&herr,&lerr);
IA_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&herr,&lerr);
IA_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&herr,&lerr);
IA_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&herr,&lerr);
IA_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&herr,&lerr);
IA_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&herr,&lerr);
IA_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&herr,&lerr);
IA_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&herr,&lerr);
IA_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&herr,&lerr);
IA_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&herr,&lerr);
IA_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&herr,&lerr);
IA_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&herr,&lerr);
IA_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&herr,&lerr);
IA_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&herr,&lerr);
IA_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&herr,&lerr);
IA_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&herr,&lerr);
IA_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&herr,&lerr);
IA_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&herr,&lerr);
IA_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&herr,&lerr);
IA_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&herr,&lerr);
IA_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&herr,&lerr);
IA_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&herr,&lerr);
IA_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&herr,&lerr);
IA_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&herr,&lerr);
IA_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&herr,&lerr);
IA_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&herr,&lerr);
IA_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&herr,&lerr);
IA_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&herr,&lerr);
IA_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&herr,&lerr);
IA_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&herr,&lerr);
IA_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&herr,&lerr);
IA_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&herr,&lerr);
IA_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&herr,&lerr);
IA_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&herr,&lerr);
IA_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&herr,&lerr);
IA_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&herr,&lerr);
IA_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&herr,&lerr);
IA_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&herr,&lerr);
IA_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&herr,&lerr);
IA_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&herr,&lerr);
IA_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&herr,&lerr);
IA_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&herr,&lerr);
IA_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&herr,&lerr);
IA_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&herr,&lerr);
IA_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&herr,&lerr);
IA_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&herr,&lerr);
IA_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&herr,&lerr);
IA_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&herr,&lerr);
IA_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&herr,&lerr);
IA_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&herr,&lerr);
IA_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&herr,&lerr);
IA_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&herr,&lerr);
IA_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&herr,&lerr);
IA_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&herr,&lerr);
IA_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&herr,&lerr);
IA_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&herr,&lerr);
IA_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&herr,&lerr);
IA_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&herr,&lerr);
IA_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&herr,&lerr);
IA_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&herr,&lerr);
IA_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&herr,&lerr);
IA_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&herr,&lerr);
IA_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&herr,&lerr);
IA_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&herr,&lerr);
IA_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&herr,&lerr);
IA_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&herr,&lerr);
IA_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&herr,&lerr);
IA_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&herr,&lerr);
IA_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&herr,&lerr);
IA_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&herr,&lerr);
IA_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&herr,&lerr);
IA_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&herr,&lerr);
IA_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&herr,&lerr);
IA_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&herr,&lerr);
IA_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&herr,&lerr);
IA_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&herr,&lerr);
IA_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&herr,&lerr);
IA_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&herr,&lerr);
IA_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&herr,&lerr);
IA_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&herr,&lerr);
IA_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&herr,&lerr);
IA_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&herr,&lerr);
IA_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&herr,&lerr);
IA_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&herr,&lerr);
IA_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&herr,&lerr);
IA_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&herr,&lerr);
IA_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&herr,&lerr);
IA_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&herr,&lerr);
IA_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&herr,&lerr);
IA_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&herr,&lerr);
IA_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&herr,&lerr);
IA_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&herr,&lerr);
IA_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&herr,&lerr);
IA_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&herr,&lerr);
IA_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&herr,&lerr);
IA_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&herr,&lerr);
IA_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&herr,&lerr);
IA_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&herr,&lerr);
IA_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&herr,&lerr);
IA_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&herr,&lerr);
IA_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&herr,&lerr);
IA_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&herr,&lerr);
IA_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&herr,&lerr);
IA_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&herr,&lerr);
IA_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&herr,&lerr);
IA_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&herr,&lerr);
IA_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&herr,&lerr);
IA_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&herr,&lerr);
IA_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&herr,&lerr);
IA_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&herr,&lerr);
IA_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&herr,&lerr);
IA_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&herr,&lerr);
IA_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&herr,&lerr);
IA_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&herr,&lerr);
IA_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&herr,&lerr);
IA_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&herr,&lerr);
IA_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&herr,&lerr);
IA_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&herr,&lerr);
IA_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&herr,&lerr);
IA_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&herr,&lerr);
IA_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&herr,&lerr);
IA_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&herr,&lerr);
IA_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&herr,&lerr);
IA_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&herr,&lerr);
IA_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&herr,&lerr);
IA_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&herr,&lerr);
IA_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&herr,&lerr);
IA_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&herr,&lerr);
IA_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&herr,&lerr);
IA_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&herr,&lerr);
IA_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&herr,&lerr);
IA_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&herr,&lerr);
IA_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&herr,&lerr);
IA_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&herr,&lerr);
IA_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&herr,&lerr);
IA_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&herr,&lerr);
IA_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&herr,&lerr);
IA_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&herr,&lerr);
IA_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&herr,&lerr);
IA_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&herr,&lerr);
IA_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&herr,&lerr);
IA_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&herr,&lerr);
IA_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&herr,&lerr);
IA_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&herr,&lerr);
IA_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&herr,&lerr);
IA_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&herr,&lerr);
IA_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&herr,&lerr);
IA_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&herr,&lerr);
IA_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&herr,&lerr);
IA_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&herr,&lerr);
IA_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&herr,&lerr);
IA_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&herr,&lerr);
IA_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&herr,&lerr);
IA_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&herr,&lerr);
IA_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&herr,&lerr);
IA_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&herr,&lerr);
IA_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&herr,&lerr);
IA_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&herr,&lerr);
IA_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&herr,&lerr);
IA_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&herr,&lerr);
IA_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&herr,&lerr);
IA_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&herr,&lerr);
IA_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&herr,&lerr);
IA_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&herr,&lerr);
IA_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&herr,&lerr);
IA_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&herr,&lerr);
IA_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&herr,&lerr);
IA_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&herr,&lerr);
IA_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&herr,&lerr);
IA_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&herr,&lerr);
IA_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&herr,&lerr);
IA_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&herr,&lerr);
IA_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&herr,&lerr);
IA_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&herr,&lerr);
IA_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&herr,&lerr);
IA_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&herr,&lerr);
IA_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&herr,&lerr);
IA_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&herr,&lerr);
IA_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&herr,&lerr);
IA_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&herr,&lerr);
IA_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&herr,&lerr);
IA_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&herr,&lerr);
IA_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&herr,&lerr);
IA_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&herr,&lerr);
IA_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&herr,&lerr);
IA_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&herr,&lerr);
IA_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&herr,&lerr);
IA_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&herr,&lerr);
IA_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&herr,&lerr);
IA_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&herr,&lerr);
IA_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&herr,&lerr);
IA_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&herr,&lerr);
IA_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&herr,&lerr);
IA_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&herr,&lerr);
IA_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&herr,&lerr);
IA_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&herr,&lerr);
IA_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&herr,&lerr);
IA_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&herr,&lerr);
IA_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&herr,&lerr);
IA_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&herr,&lerr);
IA_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&herr,&lerr);
IA_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&herr,&lerr);
IA_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&herr,&lerr);
IA_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&herr,&lerr);
IA_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&herr,&lerr);
IA_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&herr,&lerr);
IA_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&herr,&lerr);
IA_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&herr,&lerr);
IA_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&herr,&lerr);
IA_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&herr,&lerr);
IA_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&herr,&lerr);
IA_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&herr,&lerr);
IA_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&herr,&lerr);
IA_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&herr,&lerr);
IA_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&herr,&lerr);
IA_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&herr,&lerr);
IA_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&herr,&lerr);
IA_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&herr,&lerr);
IA_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&herr,&lerr);
IA_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&herr,&lerr);
IA_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&herr,&lerr);
IA_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&herr,&lerr);
IA_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&herr,&lerr);
IA_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&herr,&lerr);
IA_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&herr,&lerr);
IA_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&herr,&lerr);
IA_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&herr,&lerr);
IA_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&herr,&lerr);
IA_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&herr,&lerr);
IA_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&herr,&lerr);
IA_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&herr,&lerr);
IA_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&herr,&lerr);
IA_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&herr,&lerr);
IA_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&herr,&lerr);
IA_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&herr,&lerr);
IA_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&herr,&lerr);
IA_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&herr,&lerr);
IA_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&herr,&lerr);
IA_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&herr,&lerr);
IA_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&herr,&lerr);
IA_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&herr,&lerr);
IA_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&herr,&lerr);
IA_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&herr,&lerr);
IA_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&herr,&lerr);
IA_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&herr,&lerr);
IA_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&herr,&lerr);
IA_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&herr,&lerr);
IA_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&herr,&lerr);
IA_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&herr,&lerr);
IA_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&herr,&lerr);
IA_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&herr,&lerr);
IA_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&herr,&lerr);
IA_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&herr,&lerr);
IA_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&herr,&lerr);
IA_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&herr,&lerr);
IA_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&herr,&lerr);
IA_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&herr,&lerr);
IA_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&herr,&lerr);
IA_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&herr,&lerr);
IA_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&herr,&lerr);
IA_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&herr,&lerr);
IA_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&herr,&lerr);
IA_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&herr,&lerr);
IA_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&herr,&lerr);
IA_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&herr,&lerr);
IA_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&herr,&lerr);
IA_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&herr,&lerr);
IA_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&herr,&lerr);
IA_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&herr,&lerr);
IA_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&herr,&lerr);
IA_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&herr,&lerr);
IA_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&herr,&lerr);
IA_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&herr,&lerr);
IA_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&herr,&lerr);
IA_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&herr,&lerr);
IA_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&herr,&lerr);
IA_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&herr,&lerr);
IA_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&herr,&lerr);
IA_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&herr,&lerr);
IA_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&herr,&lerr);
IA_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&herr,&lerr);
IA_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&herr,&lerr);
IA_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&herr,&lerr);
IA_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&herr,&lerr);
IA_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&herr,&lerr);
IA_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&herr,&lerr);
IA_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&herr,&lerr);
IA_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&herr,&lerr);
IA_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&herr,&lerr);
IA_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&herr,&lerr);
IA_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&herr,&lerr);
IA_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&herr,&lerr);
IA_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&herr,&lerr);
IA_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&herr,&lerr);
IA_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&herr,&lerr);
IA_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&herr,&lerr);
IA_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&herr,&lerr);
IA_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&herr,&lerr);
IA_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&herr,&lerr);
IA_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&herr,&lerr);
IA_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&herr,&lerr);
IA_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&herr,&lerr);
IA_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&herr,&lerr);
IA_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&herr,&lerr);
IA_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&herr,&lerr);
IA_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&herr,&lerr);
IA_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&herr,&lerr);
IA_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&herr,&lerr);
IA_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&herr,&lerr);
IA_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&herr,&lerr);
IA_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&herr,&lerr);
IA_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&herr,&lerr);
IA_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&herr,&lerr);
IA_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&herr,&lerr);
IA_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&herr,&lerr);
IA_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&herr,&lerr);
IA_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&herr,&lerr);
IA_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&herr,&lerr);
IA_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&herr,&lerr);
IA_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&herr,&lerr);
IA_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&herr,&lerr);
IA_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&herr,&lerr);
IA_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&herr,&lerr);
IA_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&herr,&lerr);
IA_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&herr,&lerr);
IA_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&herr,&lerr);
IA_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&herr,&lerr);
IA_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&herr,&lerr);
IA_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&herr,&lerr);
IA_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&herr,&lerr);
IA_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&herr,&lerr);
IA_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&herr,&lerr);
IA_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&herr,&lerr);
IA_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&herr,&lerr);
IA_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&herr,&lerr);
IA_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&herr,&lerr);
IA_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&herr,&lerr);
IA_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&herr,&lerr);
IA_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&herr,&lerr);
IA_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&herr,&lerr);
IA_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&herr,&lerr);
IA_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&herr,&lerr);
IA_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&herr,&lerr);
IA_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&herr,&lerr);
IA_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&herr,&lerr);
IA_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&herr,&lerr);
IA_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&herr,&lerr);
IA_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&herr,&lerr);
IA_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&herr,&lerr);
IA_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&herr,&lerr);
IA_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&herr,&lerr);
IA_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&herr,&lerr);
IA_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&herr,&lerr);
IA_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&herr,&lerr);
IA_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&herr,&lerr);
IA_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&herr,&lerr);
IA_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&herr,&lerr);
IA_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&herr,&lerr);
IA_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&herr,&lerr);
IA_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&herr,&lerr);
IA_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&herr,&lerr);
IA_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&herr,&lerr);
IA_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&herr,&lerr);
IA_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&herr,&lerr);
IA_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&herr,&lerr);
IA_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&herr,&lerr);
IA_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&herr,&lerr);
IA_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&herr,&lerr);
IA_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&herr,&lerr);
IA_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&herr,&lerr);
IA_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&herr,&lerr);
IA_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&herr,&lerr);
IA_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&herr,&lerr);
IA_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&herr,&lerr);
IA_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&herr,&lerr);
IA_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&herr,&lerr);
IA_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&herr,&lerr);
IA_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&herr,&lerr);
IA_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&herr,&lerr);
IA_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&herr,&lerr);
IA_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&herr,&lerr);
IA_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&herr,&lerr);
IA_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&herr,&lerr);
IA_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&herr,&lerr);
IA_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&herr,&lerr);
IA_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&herr,&lerr);
IA_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&herr,&lerr);
IA_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&herr,&lerr);
IA_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&herr,&lerr);
IA_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&herr,&lerr);
IA_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&herr,&lerr);
IA_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&herr,&lerr);
IA_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&herr,&lerr);
IA_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&herr,&lerr);
IA_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&herr,&lerr);
IA_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&herr,&lerr);
IA_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&herr,&lerr);
IA_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&herr,&lerr);
IA_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&herr,&lerr);
IA_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&herr,&lerr);
IA_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&herr,&lerr);
IA_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&herr,&lerr);
IA_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&herr,&lerr);
IA_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&herr,&lerr);
IA_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&herr,&lerr);
IA_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&herr,&lerr);
IA_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&herr,&lerr);
IA_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&herr,&lerr);
IA_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&herr,&lerr);
IA_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&herr,&lerr);
IA_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&herr,&lerr);
IA_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&herr,&lerr);
IA_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&herr,&lerr);
IA_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&herr,&lerr);
IA_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&herr,&lerr);
IA_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&herr,&lerr);
IA_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&herr,&lerr);
IA_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&herr,&lerr);
IA_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&herr,&lerr);
IA_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&herr,&lerr);
IA_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&herr,&lerr);
IA_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&herr,&lerr);
IA_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&herr,&lerr);
IA_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&herr,&lerr);
IA_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&herr,&lerr);
IA_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&herr,&lerr);
IA_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&herr,&lerr);
IA_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&herr,&lerr);
IA_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&herr,&lerr);
IA_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&herr,&lerr);
IA_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&herr,&lerr);
IA_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&herr,&lerr);
IA_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&herr,&lerr);
IA_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&herr,&lerr);
IA_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&herr,&lerr);
IA_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&herr,&lerr);
IA_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&herr,&lerr);
IA_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&herr,&lerr);
IA_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&herr,&lerr);
IA_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&herr,&lerr);
IA_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&herr,&lerr);
IA_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&herr,&lerr);
IA_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&herr,&lerr);
IA_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&herr,&lerr);
IA_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&herr,&lerr);
IA_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&herr,&lerr);
IA_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&herr,&lerr);
IA_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&herr,&lerr);
IA_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&herr,&lerr);
IA_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&herr,&lerr);
IA_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&herr,&lerr);
IA_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&herr,&lerr);
IA_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&herr,&lerr);
IA_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&herr,&lerr);
IA_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&herr,&lerr);
IA_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&herr,&lerr);
IA_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&herr,&lerr);
IA_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&herr,&lerr);
IA_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&herr,&lerr);
IA_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&herr,&lerr);
IA_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&herr,&lerr);
IA_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&herr,&lerr);
IA_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&herr,&lerr);
IA_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&herr,&lerr);
IA_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&herr,&lerr);
IA_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&herr,&lerr);
IA_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&herr,&lerr);
IA_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&herr,&lerr);
IA_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&herr,&lerr);
IA_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&herr,&lerr);
IA_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&herr,&lerr);
IA_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&herr,&lerr);
IA_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&herr,&lerr);
IA_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&herr,&lerr);
IA_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&herr,&lerr);
IA_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&herr,&lerr);
IA_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&herr,&lerr);
IA_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&herr,&lerr);
IA_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&herr,&lerr);
IA_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&herr,&lerr);
IA_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&herr,&lerr);
IA_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&herr,&lerr);
IA_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&herr,&lerr);
IA_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&herr,&lerr);
IA_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&herr,&lerr);
IA_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&herr,&lerr);
IA_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&herr,&lerr);
IA_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&herr,&lerr);
IA_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&herr,&lerr);
IA_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&herr,&lerr);
IA_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&herr,&lerr);
IA_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&herr,&lerr);
IA_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&herr,&lerr);
IA_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&herr,&lerr);
IA_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&herr,&lerr);
IA_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&herr,&lerr);
IA_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&herr,&lerr);
IA_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&herr,&lerr);
IA_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&herr,&lerr);
IA_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&herr,&lerr);
IA_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&herr,&lerr);
IA_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&herr,&lerr);
IA_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&herr,&lerr);
IA_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&herr,&lerr);
IA_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&herr,&lerr);
IA_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&herr,&lerr);
IA_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&herr,&lerr);
IA_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&herr,&lerr);
IA_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&herr,&lerr);
IA_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&herr,&lerr);
IA_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&herr,&lerr);
IA_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&herr,&lerr);
IA_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&herr,&lerr);
IA_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&herr,&lerr);
IA_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&herr,&lerr);
IA_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&herr,&lerr);
IA_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&herr,&lerr);
IA_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&herr,&lerr);
IA_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&herr,&lerr);
IA_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&herr,&lerr);
IA_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&herr,&lerr);
IA_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&herr,&lerr);
IA_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&herr,&lerr);
IA_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&herr,&lerr);
IA_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&herr,&lerr);
IA_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&herr,&lerr);
IA_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&herr,&lerr);
IA_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&herr,&lerr);
IA_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&herr,&lerr);
IA_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&herr,&lerr);
IA_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&herr,&lerr);
IA_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&herr,&lerr);
IA_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&herr,&lerr);
IA_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&herr,&lerr);
IA_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&herr,&lerr);
IA_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&herr,&lerr);
IA_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&herr,&lerr);
IA_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&herr,&lerr);
IA_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&herr,&lerr);
IA_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&herr,&lerr);
IA_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&herr,&lerr);
IA_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&herr,&lerr);
IA_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&herr,&lerr);
IA_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&herr,&lerr);
IA_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&herr,&lerr);
IA_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&herr,&lerr);
IA_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&herr,&lerr);
IA_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&herr,&lerr);
IA_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&herr,&lerr);
IA_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&herr,&lerr);
IA_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&herr,&lerr);
IA_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&herr,&lerr);
IA_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&herr,&lerr);
IA_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&herr,&lerr);
IA_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&herr,&lerr);
IA_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&herr,&lerr);
IA_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&herr,&lerr);
IA_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&herr,&lerr);
IA_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&herr,&lerr);
IA_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&herr,&lerr);
IA_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&herr,&lerr);
IA_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&herr,&lerr);
IA_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&herr,&lerr);
IA_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&herr,&lerr);
IA_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&herr,&lerr);
IA_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&herr,&lerr);
IA_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&herr,&lerr);
IA_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&herr,&lerr);
IA_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&herr,&lerr);
IA_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&herr,&lerr);
IA_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&herr,&lerr);
IA_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&herr,&lerr);
IA_Intersection3D(ax1,ax2,ax3,ay1,ay2,ay3,az1,az2,az3,bx1,bx2,bx3,by1,by2,by3,bz1,bz2,bz3,cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,&res,&herr,&lerr);
after = clock();
IAtime = ((double)(after - before))/CLOCKS_PER_SEC;
IAtotaltime = IAtotaltime + IAtime;

printf("Total time FPS: %e \n", FPStime);
printf("Total time IA: %e \n\n", IAtime);

fprintf(save,"%e\n",FPStime);
fprintf(save,"%e\n",IAtime);
}

printf("Average FPS time : %f\n",FPStotaltime/sample);
printf("Average IA time : %f\n",IAtotaltime/sample);


fprintf(report,"Average FPS time : %f\n",FPStotaltime/sample);
fprintf(report,"Average IA time : %f\n",IAtotaltime/sample);

fclose(save);
fclose(report);
}