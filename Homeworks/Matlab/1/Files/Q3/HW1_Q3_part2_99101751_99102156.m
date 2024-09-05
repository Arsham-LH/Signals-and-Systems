%% Q3.2
%part 1
clear;clc;close all;
syms n z;
oldParam=sympref('HeavisideAtOrigin',1); %to change the value at origin to 1
x(n) = cospi(n/6)*heaviside(n); %cospi(n/6)=cos(n*pi/6)
x_z(n) = cospi(n/6); %because ztrans is not bilateral
X(z)=ztrans(x_z,z);
disp(expand(X(z)));
[X_num,X_den] = numden(X);
X_num_coeffs=double(coeffs(X_num(z),'All'));
X_den_coeffs=double(coeffs(X_den(z),'All'));
X_tf = tf(X_num_coeffs,X_den_coeffs);
X_zeros = zero(X_tf);
X_poles = pole(X_tf);
figure('Name','x[n]&zero-pole plots for X(z) (or X0(z))');
r0 = [-1:13];
subplot(1,2,1);
stem(r0,x(r0),'markerFaceColor','Blue');
title("x[n]");
xlabel('n');
ylabel('x[n]');
grid;
subplot(1,2,2);
zplane(X_zeros,X_poles);
title('zero-pole plot');
%part 2
X0(z) = X(z);
X1(z)=X0(2*z);
[X1_num,X1_den] = numden(X1);
X1_num_coeffs=double(coeffs(X1_num(z),'All'));
X1_den_coeffs=double(coeffs(X1_den(z),'All'));
X1_tf = tf(X1_num_coeffs,X1_den_coeffs);
X1_zeros = zero(X1_tf);
X1_poles = pole(X1_tf);

x1(n) = iztrans(X1(z))*heaviside(n);
figure('Name','x1[n]&zero-pole plots for X1(z)');
r1 = [-1:13];
subplot(1,2,1);
stem(r1,x1(r1),'markerFaceColor','Blue');
title("x1[n]");
xlabel('n');
ylabel('x1[n]');
grid;
subplot(1,2,2);
zplane(X1_zeros,X1_poles);
title('zero-pole plot');
grid;

%part 3
X2(z)=X0(z^3);
[X2_num,X2_den] = numden(X2);
X2_num_coeffs=double(coeffs(X2_num(z),'All'));
X2_den_coeffs=double(coeffs(X2_den(z),'All'));
X2_tf = tf(X2_num_coeffs,X2_den_coeffs);
X2_zeros = zero(X2_tf);
X2_poles = pole(X2_tf);

x2(n) = iztrans(X2(z))*heaviside(n);
figure('Name','x2[n]&zero-pole plots for X2(z)');
r2 = [-1:37];
subplot(1,2,1);
stem(r2,x2(r2),'markerFaceColor','Blue');
title("x2[n]");
xlabel('n');
ylabel('x2[n]');
grid;
subplot(1,2,2);
zplane(X2_zeros,X2_poles);
title('zero-pole plot');
grid;




