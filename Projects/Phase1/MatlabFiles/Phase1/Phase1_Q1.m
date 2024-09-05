%% part1
clear;clc;close all;
y=cell2mat(struct2cell(load('y.mat'))); %loading file and converting to vector
y_dft=fft(y); %calculating fft
y_dft=fftshift(y_dft); %shifting to the center
N=size(y,2); %vector size
omega=-(N-1)/2:(N-1)/2; %since N is odd
mag=abs(y_dft/N); %fft magnitude
phas=angle(y_dft); %fft phase

%plots
figure('Name','Magnitude&Phase diagrams for y');
subplot(1,2,1);
plot(omega,mag);
title('fft Magnitude');
subplot(1,2,2);
plot(omega,phas);
title('fft phase');