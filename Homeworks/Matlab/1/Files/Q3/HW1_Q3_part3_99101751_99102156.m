%% Q3.3
%part 1
clear;clc;close all;
syms z n;
H1(z) = (1-z^(-1))/(1-z^(-1)+0.5*z^(-2));
H2(z) = z^(-1)/(1-z^(-1)+0.5*z^(-2));

[H1_num,H1_den] = numden(H1);
H1_num_coeffs=double(coeffs(H1_num(z),'All'));
H1_den_coeffs=double(coeffs(H1_den(z),'All'));
H1_tf = tf(H1_num_coeffs,H1_den_coeffs);
H1_zeros = zero(H1_tf);
H1_poles = pole(H1_tf);

[H2_num,H2_den] = numden(H2);
H2_num_coeffs=double(coeffs(H2_num(z),'All'));
H2_den_coeffs=double(coeffs(H2_den(z),'All'));
H2_tf = tf(H2_num_coeffs,H2_den_coeffs);
H2_zeros = zero(H2_tf);
H2_poles = pole(H2_tf);

figure('Name','zero-pole plots');
subplot(1,2,1);
zplane(H1_zeros,H1_poles);
title('H1(z)');
grid;
subplot(1,2,2);
zplane(H2_zeros,H2_poles);
title('H2(z)');
grid;

%part 2
b1 = [1 -1]; %numerator coeffs of H1(z)
a1 = [1 -1 0.5]; %denominator coeffs of H2(z)
[r1,p1,k1]=residuez(b1,a1);

b2 = [0 1];
a2 = [1 -1 0.5];
[r2,p2,k2]=residuez(b2,a2);
%part 3
h1(n) = simplify(iztrans(H1(z)));
h2(n) = simplify(iztrans(H2(z)));
