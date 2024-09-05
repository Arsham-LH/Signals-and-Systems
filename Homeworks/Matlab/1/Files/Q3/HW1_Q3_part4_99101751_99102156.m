%% Q3.4.1
clear;clc;close all;
syms z n;
oldParam=sympref('HeavisideAtOrigin',1); %to change the value at origin to 1

H(z) = (2-z^(-1))/(1-0.7*z^(-1)+0.49*z^(-2));
b=[2 -1]; %numerator coeffs of H(z)
a=[1 -0.7 0.49]; %denominator coeffs of H(z)
[r,p,k]=residuez(b,a);
disp(r);
disp(p);
disp(k);
h(n) = r(1)*(p(1))^n+r(2)*(p(2))^n;
h_simp(n) = 2.06*(0.7)^n*cos(0.24+n*pi/3); %h[n], after simplification, with small approximation
disp(h(n));
figure('Name','h[n]');
n1=[-1:6];
stem(n1,h_simp(n1).*heaviside(n1),'markerFaceColor','Blue');
title('h[n]');
xlabel('n');
ylabel('h[n]');
grid;
%% Q3.4.2
clear;clc;close all;
syms n;
h(n) = (1+3^(1/2)/7*j)*(7/20*(1+3^(1/2)*j))^n+(1-3^(1/2)/7*j)*(7/20*(1-3^(1/2)*j))^n;
figure('Name','h[n]');
n1=[-1:6];
stem(n1,h(n1).*heaviside(n1),'markerFaceColor','Blue');
title('h[n]');
xlabel('n');
ylabel('h[n]');
grid;
%% Q3.4.3
clear;clc;close all;
b = [2 -1];
a=[1 -0.7 0.49];
t=linspace(-1,15,17);
x=(t==0);
h=filter(b,a,x);
figure('Name','h[n] using filter');
stem(t,h,'markerFaceColor','Blue');
title('h[n]');
xlabel('n');
ylabel('h[n]');
grid;









