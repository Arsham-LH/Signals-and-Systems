%% Q3.1.1
clear;clc;close all;
syms n z;
oldParam=sympref('HeavisideAtOrigin',1);
n1=[-2:35];
n2=[-1:45];
n3=[-1:10];
n4=[-1:10];

h1(n) = ((0.5).^n).*(heaviside(n)-heaviside(n-30));
h2(n) = ((2).^n).*(heaviside(n)-heaviside(n-40));
h3(n) = ((0.5).^n).*heaviside(n);
h4(n) = ((2).^n).*heaviside(n);

figure('Name','descrete signals');
subplot(2,2,1);
stem(n1,h1(n1),'markerFaceColor','Blue');
title('h1[n]');
xlabel('n');
ylabel('h1[n]');
grid;
subplot(2,2,2);
stem(n2,h2(n2),'markerFaceColor','Blue');
title('h2[n]');
xlabel('n');
ylabel('h2[n]');
grid;
subplot(2,2,3);
stem(n3,h3(n3),'markerFaceColor','Blue');
title('h3[n]');
xlabel('n');
ylabel('h3[n]');
grid;
subplot(2,2,4);
stem(n4,h4(n4),'markerFaceColor','Blue');
title('h4[n]');
xlabel('n');
ylabel('h4[n]');
grid;

H1(z)=ztrans(h1,z);
H2(z)=ztrans(h2,z);
H3(z)=ztrans(h3,z);
H4(z)=ztrans(h4,z);

[H1_num(z),H1_den(z)] = numden(H1);
[H2_num(z),H2_den(z)] = numden(H2);
[H3_num(z),H3_den(z)] = numden(H3);
[H4_num(z),H4_den(z)] = numden(H4);

H1_num_coeffs=double(coeffs(H1_num(z),'All'));
H2_num_coeffs=double(coeffs(H2_num(z),'All'));
H3_num_coeffs=double(coeffs(H3_num(z),'All'));
H4_num_coeffs=double(coeffs(H4_num(z),'All'));

H1_den_coeffs=double(coeffs(H1_den(z),'All'));
H2_den_coeffs=double(coeffs(H2_den(z),'All'));
H3_den_coeffs=double(coeffs(H3_den(z),'All'));
H4_den_coeffs=double(coeffs(H4_den(z),'All'));

H1_tf=tf(H1_num_coeffs,H1_den_coeffs);
H2_tf=tf(H2_num_coeffs,H2_den_coeffs);
H3_tf=tf(H3_num_coeffs,H3_den_coeffs);
H4_tf=tf(H4_num_coeffs,H4_den_coeffs);

H1_zeros=zero(H1_tf);
H2_zeros=zero(H2_tf);
H3_zeros=zero(H3_tf);
H4_zeros=zero(H4_tf);

H1_poles=pole(H1_tf);
H2_poles=pole(H2_tf);
H3_poles=pole(H3_tf);
H4_poles=pole(H4_tf);

figure('Name','zero-pole plots');
subplot(2,2,1);
zplane(H1_zeros,H1_poles);
title('H1(z)');
subplot(2,2,2);
zplane(H2_zeros,H2_poles);
title('H2(z)');
subplot(2,2,3);
zplane(H3_zeros,H3_poles);
title('H3(z)');
subplot(2,2,4);
zplane(H4_zeros,H4_poles);
title('H4(z)');

%% Q3.1.2
clear;clc;close all;
oldParam=sympref('HeavisideAtOrigin',1);

syms n z;
X1(z)=1/(1-0.5*z^(-1));
X2(z)=1/(1-1.1*z^(-1));
X3(z)=1/(1-z^(-1));
X4(z)=1/(1-0.9*z^(-1));
X5(z)=1/(1-5*z^(-1));
X6(z)=1/((1-z^(-1))^2);

[X1_num(z),X1_den(z)] = numden(X1);
[X2_num(z),X2_den(z)] = numden(X2);
[X3_num(z),X3_den(z)] = numden(X3);
[X4_num(z),X4_den(z)] = numden(X4);
[X5_num(z),X5_den(z)] = numden(X5);
[X6_num(z),X6_den(z)] = numden(X6);

X1_num_coeffs=double(coeffs(X1_num(z),'All'));
X2_num_coeffs=double(coeffs(X2_num(z),'All'));
X3_num_coeffs=double(coeffs(X3_num(z),'All'));
X4_num_coeffs=double(coeffs(X4_num(z),'All'));
X5_num_coeffs=double(coeffs(X5_num(z),'All'));
X6_num_coeffs=double(coeffs(X6_num(z),'All'));

X1_den_coeffs=double(coeffs(X1_den(z),'All'));
X2_den_coeffs=double(coeffs(X2_den(z),'All'));
X3_den_coeffs=double(coeffs(X3_den(z),'All'));
X4_den_coeffs=double(coeffs(X4_den(z),'All'));
X5_den_coeffs=double(coeffs(X5_den(z),'All'));
X6_den_coeffs=double(coeffs(X6_den(z),'All'));

X1_tf=tf(X1_num_coeffs,X1_den_coeffs);
X2_tf=tf(X2_num_coeffs,X2_den_coeffs);
X3_tf=tf(X3_num_coeffs,X3_den_coeffs);
X4_tf=tf(X4_num_coeffs,X4_den_coeffs);
X5_tf=tf(X5_num_coeffs,X5_den_coeffs);
X6_tf=tf(X6_num_coeffs,X6_den_coeffs);

X1_zeros=zero(X1_tf);
X2_zeros=zero(X2_tf);
X3_zeros=zero(X3_tf);
X4_zeros=zero(X4_tf);
X5_zeros=zero(X5_tf);
X6_zeros=zero(X6_tf);

X1_poles=pole(X1_tf);
X2_poles=pole(X2_tf);
X3_poles=pole(X3_tf);
X4_poles=pole(X4_tf);
X5_poles=pole(X5_tf);
X6_poles=pole(X6_tf);

figure('Name','zero-pole&signal plots');
subplot(2,6,1);
zplane(X1_zeros,X1_poles);
title('X1(z)');
subplot(2,6,2);
zplane(X2_zeros,X2_poles);
title('X2(z)');
subplot(2,6,3);
zplane(X3_zeros,X3_poles);
title('X3(z)');
subplot(2,6,4);
zplane(X4_zeros,X4_poles);
title('X4(z)');
subplot(2,6,5);
zplane(X5_zeros,X5_poles);
title('X5(z)');
subplot(2,6,6);
zplane(X6_zeros,X6_poles);
title('X6(z)');

x1(n)=iztrans(X1)*heaviside(n);
x2(n)=iztrans(X2)*heaviside(n);
x3(n)=iztrans(X3)*heaviside(n);
x4(n)=iztrans(X4)*heaviside(n);
x5(n)=iztrans(X5)*heaviside(n);
x6(n)=iztrans(X6)*heaviside(n);

r = [-1:10];
subplot(2,6,7);
stem(r,x1(r),'markerFaceColor','Blue');
title('x1[n]');
xlabel('n');
ylabel('x1[n]');
grid;
subplot(2,6,8);
stem(r,x2(r),'markerFaceColor','Blue');
title('x2[n]');
xlabel('n');
ylabel('x2[n]');
grid;
subplot(2,6,9);
stem(r,x3(r),'markerFaceColor','Blue');
title('x3[n]');
xlabel('n');
ylabel('x3[n]');
grid;
subplot(2,6,10);
stem(r,x4(r),'markerFaceColor','Blue');
title('x4[n]');
xlabel('n');
ylabel('x4[n]');
grid;
subplot(2,6,11);
stem(r,x5(r),'markerFaceColor','Blue');
title('x5[n]');
xlabel('n');
ylabel('x5[n]');
grid;
subplot(2,6,12);
stem(r,x6(r),'markerFaceColor','Blue');
title('x6[n]');
xlabel('n');
ylabel('x6[n]');
grid;

%% Q3.1.3
clear;clc;close all;
oldParam=sympref('HeavisideAtOrigin',1);

syms n z;
X7(z)=1/(1-sqrt(2)*z^(-1)+z^(-2));
X8(z)=1/(1-z^(-1)+2*z^(-2));
X9(z)=1/(2-z^(-1)+z^(-2));

[X7_num(z),X7_den(z)] = numden(X7);
[X8_num(z),X8_den(z)] = numden(X8);
[X9_num(z),X9_den(z)] = numden(X9);

X7_num_coeffs=double(coeffs(X7_num(z),'All'));
X8_num_coeffs=double(coeffs(X8_num(z),'All'));
X9_num_coeffs=double(coeffs(X9_num(z),'All'));

X7_den_coeffs=double(coeffs(X7_den(z),'All'));
X8_den_coeffs=double(coeffs(X8_den(z),'All'));
X9_den_coeffs=double(coeffs(X9_den(z),'All'));

X7_tf=tf(X7_num_coeffs,X7_den_coeffs);
X8_tf=tf(X8_num_coeffs,X8_den_coeffs);
X9_tf=tf(X9_num_coeffs,X9_den_coeffs);

X7_zeros=zero(X7_tf);
X8_zeros=zero(X8_tf);
X9_zeros=zero(X9_tf);

X7_poles=pole(X7_tf);
X8_poles=pole(X8_tf);
X9_poles=pole(X9_tf);

figure('Name','zero-pole&signal plots');
subplot(2,3,1);
zplane(X7_zeros,X7_poles);
title('X7(z)');
subplot(2,3,2);
zplane(X8_zeros,X8_poles);
title('X8(z)');
subplot(2,3,3);
zplane(X9_zeros,X9_poles);
title('X9(z)');

x7(n)=iztrans(X7)*heaviside(n);
x8(n)=iztrans(X8)*heaviside(n);
x9(n)=iztrans(X9)*heaviside(n);

r = [-1:10];
subplot(2,3,4);
stem(r,x7(r),'markerFaceColor','Blue');
title('x7[n]');
xlabel('n');
ylabel('x7[n]');
grid;
subplot(2,3,5);
stem(r,x8(r),'markerFaceColor','Blue');
title('x8[n]');
xlabel('n');
ylabel('x8[n]');
grid;
subplot(2,3,6);
stem(r,x9(r),'markerFaceColor','Blue');
title('x9[n]');
xlabel('n');
ylabel('x9[n]');
grid;

