%% Q2_1&2. ghabl az run kardane haryek az section haye Q3.1 va Q3.2, ebteda in section run shavad
clear;clc;close all;
syms s;
%defining transfer functions
H1(s)=1/(s^3+40*s^2+10*s+500);
H2(s)=1/(s^4+12.5*s^3+10*s^2+10*s+1);
H3(s)=1/(s^5+125*s^4+100*s^3+100*s^2+20*s+10);
H4(s)=1/(s^6+5*s^5+125*s^4+100*s^3+100*s^2+20*s+10);
%get numerator and denominators
[H1_num(s),H1_den(s)] = numden(H1);
[H2_num(s),H2_den(s)] = numden(H2);
[H3_num(s),H3_den(s)] = numden(H3);
[H4_num(s),H4_den(s)] = numden(H4);
%get numerator and denominators coefficients
H1_den_coeffs=double(coeffs(H1_den(s),'All'));
H2_den_coeffs=double(coeffs(H2_den(s),'All'));
H3_den_coeffs=double(coeffs(H3_den(s),'All'));
H4_den_coeffs=double(coeffs(H4_den(s),'All'));
%creating transfer function
H1_tf=tf(1,H1_den_coeffs);
H2_tf=tf(1,H2_den_coeffs);
H3_tf=tf(1,H3_den_coeffs);
H4_tf=tf(1,H4_den_coeffs);

%% Q2.1
clc;

H1_zeros=zero(H1_tf);
H2_zeros=zero(H2_tf);
H3_zeros=zero(H3_tf);
H4_zeros=zero(H4_tf);

H1_poles=pole(H1_tf)
H2_poles=pole(H2_tf)
H3_poles=pole(H3_tf)
H4_poles=pole(H4_tf)

%plotting zeros and poles
figure('Name','zero-pole plots');
subplot(2,2,1);
zplane(H1_zeros,H1_poles);
title('H1(s)');
subplot(2,2,2);
zplane(H2_zeros,H2_poles);
title('H2(s)');
subplot(2,2,3);
zplane(H3_zeros,H3_poles);
title('H3(s)');
subplot(2,2,4);
zplane(H4_zeros,H4_poles);
title('H4(s)');

%% Q2.2
clc;
syms t;
%laplace of u(t}
X(s)=laplace(heaviside(t));
%laplace of output=laplace of input*transfer function
Y1(s)=X(s)*H1(s);
Y2(s)=X(s)*H2(s);
Y3(s)=X(s)*H3(s);
Y4(s)=X(s)*H4(s);
%inverse laplace to get outputs
y1(t)=ilaplace(Y1)*heaviside(t);
y2(t)=ilaplace(Y2)*heaviside(t);
y3(t)=ilaplace(Y3)*heaviside(t);
y4(t)=ilaplace(Y4)*heaviside(t);
%defining plots ranges
t1=linspace(-1,100,600);
t2=linspace(-1,100);
t3=linspace(-1,100);
t4=linspace(-1,100);
%plotting unit step responses
figure('Name','unit step responses');
subplot(2,2,1);
plot(t1,y1(t1));
title('y1(t)');
xlabel('t');
ylabel('y1(t)');
grid;
subplot(2,2,2);
plot(t2,y2(t2));
title('y2(t)');
xlabel('t');
ylabel('y2(t)');
grid;
subplot(2,2,3);
plot(t3,y3(t3));
title('y3(t)');
xlabel('t');
ylabel('y3(t)');
grid;
subplot(2,2,4);
plot(t4,y4(t4));
title('y4(t)');
xlabel('t');
ylabel('y4(t)');
grid;

%% Q2.4
clear;clc;close all;
syms s t;
%defining transfer function
G1(s)=(s+1)/(s^2+6*s+8);
%inputs
x1(t)=dirac(t);
x2(t)=heaviside(t);
x3(t)=sin(3*t)*heaviside(t);
x4(t)=exp(-0.5*t)*heaviside(t);
x5(t)=exp(-0.2*t)*cos(2*t)*heaviside(t);
x6(t)=t^4*exp(-0.5*t)*heaviside(t);
%laplace transform
X1(s)=laplace(x1(t));
X2(s)=laplace(x2(t));
X3(s)=laplace(x3(t));
X4(s)=laplace(x4(t));
X5(s)=laplace(x5(t));
X6(s)=laplace(x6(t));
%output laplace
Y1(s)=X1(s)*G1(s);
Y2(s)=X2(s)*G1(s);
Y3(s)=X3(s)*G1(s);
Y4(s)=X4(s)*G1(s);
Y5(s)=X5(s)*G1(s);
Y6(s)=X6(s)*G1(s);
%output
y1(t)=ilaplace(Y1)*heaviside(t);
y2(t)=ilaplace(Y2)*heaviside(t);
y3(t)=ilaplace(Y3)*heaviside(t);
y4(t)=ilaplace(Y4)*heaviside(t);
y5(t)=ilaplace(Y5)*heaviside(t);
y6(t)=ilaplace(Y6)*heaviside(t);
%plots ranges
t1=linspace(-1,5);
t2=linspace(-1,10);
t3=linspace(-1,10);
t4=linspace(-1,20);
t5=linspace(-1,35);
t6=linspace(-1,40);
%plot
figure('Name','responses to system G1');
subplot(2,3,1);
plot(t1,y1(t1));
title('y1(t)');
xlabel('t');
ylabel('y1(t)');
grid;
subplot(2,3,2);
plot(t2,y2(t2));
title('y2(t)');
xlabel('t');
ylabel('y2(t)');
grid;
subplot(2,3,3);
plot(t3,y3(t3));
title('y3(t)');
xlabel('t');
ylabel('y3(t)');
grid;
subplot(2,3,4);
plot(t4,y4(t4));
title('y4(t)');
xlabel('t');
ylabel('y4(t)');
grid;
subplot(2,3,5);
plot(t5,y5(t5));
title('y5(t)');
xlabel('t');
ylabel('y5(t)');
grid;
subplot(2,3,6);
plot(t6,y6(t6));
title('y6(t)');
xlabel('t');
ylabel('y6(t)');
grid;

%print y1(t) to y6(t)
y1(t)
y2(t)
y3(t)
y4(t)
y5(t)
y6(t)

%% Q2.5
clear;clc;close all;

syms s t;
%define transfer functions
G21(s)=(2*s+1)/(s^2+4*s+7);
G22(s)=(2*s+1)/(s^2+5*s+7);
G23(s)=(2*s+1)/(s^2+6*s+7);
%input laplace
X(s)=laplace(heaviside(t));
%outputs laplace
Y1(s)=X(s)*G21(s);
Y2(s)=X(s)*G22(s);
Y3(s)=X(s)*G23(s);
%outputs
y1(t)=ilaplace(Y1)*heaviside(t);
y2(t)=ilaplace(Y2)*heaviside(t);
y3(t)=ilaplace(Y3)*heaviside(t);
%plots ranges
t1=linspace(-1,10,1500);
t2=linspace(-1,10,1500);
t3=linspace(-1,10,1500);

figure('Name','unit response to systems G2');
subplot(1,3,1);
plot(t1,y1(t1));
title('y1(t)');
xlabel('t');
ylabel('y1(t)');
grid;
subplot(1,3,2);
plot(t2,y2(t2));
title('y2(t)');
xlabel('t');
ylabel('y2(t)');
grid;
subplot(1,3,3);
plot(t3,y3(t3));
title('y3(t)');
xlabel('t');
ylabel('y3(t)');
grid;

y1(t)
y2(t)
y3(t)

%final values for y1(t),y2(t),y3(t) at t-->Inf
final1=limit(y1,t,Inf);
final2=limit(y2,t,Inf);
final3=limit(y3,t,Inf);
finalVals=double([final1,final2,final3]);
%max values for y1(t),y2(t),y3(t)
maxVal1=max(y1(t1));
maxVal2=max(y2(t2));
maxVal3=max(y3(t3));
maxVals=double([maxVal1,maxVal2,maxVal3]);
%points on which max occurs
maxPoint1=solve(y1(t)-maxVal1);
maxPoint2=solve(y2(t)-maxVal2);
maxPoint3=solve(y3(t)-maxVal3);
maxValPoints=double([maxPoint1,maxPoint2,maxPoint3]);
%points on which y reaches final/2
assume(t>0);
point1=min(solve(y1(t)-(final1)/2));
point2=min(solve(y2(t)-(final2)/2));
point3=min(solve(y3(t)-(final3)/2));
halfPoints=double([point1,point2,point3]);

%% Q2.6
clear;clc;close all;

syms s t;

%defining transfer functions
H1(s)=(s+1)/(s^2+3*s+4);
H2(s)=(s+1)/(s^3+3*s^2+4*s);
H3(s)=(s+1)/(s^4+3*s^3+4*s^2);
%get numerator and denominators
[H1_num(s),H1_den(s)] = numden(H1);
[H2_num(s),H2_den(s)] = numden(H2);
[H3_num(s),H3_den(s)] = numden(H3);
%get numerator and denominators coefficients
H1_num_coeffs=double(coeffs(H1_num(s),'All'));
H2_num_coeffs=double(coeffs(H2_num(s),'All'));
H3_num_coeffs=double(coeffs(H3_num(s),'All'));
H1_den_coeffs=double(coeffs(H1_den(s),'All'));
H2_den_coeffs=double(coeffs(H2_den(s),'All'));
H3_den_coeffs=double(coeffs(H3_den(s),'All'));
%creating transfer function
H1_tf=tf(H1_num_coeffs,H1_den_coeffs);
H2_tf=tf(H2_num_coeffs,H2_den_coeffs);
H3_tf=tf(H3_num_coeffs,H3_den_coeffs);

%defining feedback system as fbsys(s)=1
fbsys = tf([1],[1]);
%adding negative feedback to systems
H1_fb_tf=feedback(H1_tf,fbsys,-1);
H2_fb_tf=feedback(H2_tf,fbsys,-1);
H3_fb_tf=feedback(H3_tf,fbsys,-1);
%collecting feedback system coeffs as vectors
[H1_num_fb,H1_den_fb]= tfdata(H1_fb_tf,'v');
[H2_num_fb,H2_den_fb]= tfdata(H2_fb_tf,'v');
[H3_num_fb,H3_den_fb]= tfdata(H3_fb_tf,'v');
%creating systems
H1_fb(s)=poly2sym(H1_num_fb,s)/poly2sym(H1_den_fb,s);
H2_fb(s)=poly2sym(H2_num_fb,s)/poly2sym(H2_den_fb,s);
H3_fb(s)=poly2sym(H3_num_fb,s)/poly2sym(H3_den_fb,s);

%Y1i: unit step response, Y2i: unit slope response
X1(s)=laplace(heaviside(t));
X2(s)=laplace(t*heaviside(t));
Y11(s)=H1_fb(s)*X1(s);
Y12(s)=H2_fb(s)*X1(s);
Y13(s)=H3_fb(s)*X1(s);
Y21(s)=H1_fb(s)*X2(s);
Y22(s)=H2_fb(s)*X2(s);
Y23(s)=H3_fb(s)*X2(s);
%inverse laplace
y11(t)=ilaplace(Y11(s))*heaviside(t);
y12(t)=ilaplace(Y12(s))*heaviside(t);
y13(t)=ilaplace(Y13(s))*heaviside(t);
y21(t)=ilaplace(Y21(s))*heaviside(t);
y22(t)=ilaplace(Y22(s))*heaviside(t);
y23(t)=ilaplace(Y23(s))*heaviside(t);

t1=linspace(-1,100);
t2=linspace(-1,100);
t3=linspace(-1,100);
t4=linspace(-1,100);
t5=linspace(-1,100);
t6=linspace(-1,100);

figure('Name','unit step & slope response to feedback systems');
subplot(2,3,1);
plot(t1,y11(t1));
title('y11(t)');
xlabel('t');
ylabel('y11(t)');
grid;
subplot(2,3,2);
plot(t2,y12(t2));
title('y12(t)');
xlabel('t');
ylabel('y12(t)');
grid;
subplot(2,3,3);
plot(t3,y13(t3));
title('y13(t)');
xlabel('t');
ylabel('y13(t)');
grid;
subplot(2,3,4);
plot(t4,y21(t4));
title('y21(t)');
xlabel('t');
ylabel('y21(t)');
grid;
subplot(2,3,5);
plot(t5,y22(t5));
title('y22(t)');
xlabel('t');
ylabel('y22(t)');
grid;
subplot(2,3,6);
plot(t6,y23(t6));
title('y23(t)');
xlabel('t');
ylabel('y23(t)');
grid;










