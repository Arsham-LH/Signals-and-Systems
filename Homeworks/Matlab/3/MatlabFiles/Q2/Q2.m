%% Q2.1
clear;clc;close all;
[m,fs]=audioread('Old_Mc.Donald.wav');
m=m.';
L=length(m);
Ac=0.1;
f0=20000;
u=0.7;
t=linspace(0,1,L)*L/fs;
x=Ac*(1+u*m).*cos(2*pi*f0*t);
figure;
subplot(1,2,1);
plot(t,x);
title('x(t) vs t');
xlabel('t');
ylabel('x(t)');

fft_x = fft(x); %take fourier transform
fft_x_shift = fftshift(fft_x); %shift zero-frequency component to center of spectrum
f = fs/2*linspace(-1,1,L); %the frequency interval for ploting the response, after shifting the fft, from -fs/2 to +fs/2
subplot(1,2,2);
plot(f, abs(fft_x_shift));
title('frequency response for x(t)');
xlabel('Frequency (Hz)');
ylabel('magnitude');

%% Q2.3
clc;close all;
x_noise=awgn(x,70);
figure;
subplot(1,1,1);
plot(t,x_noise);
title('noise added to x(t)');
xlabel('t');
ylabel('noisy x(t)');
audiowrite('noisy song.wav',x_noise,fs);
%% Q2.4
clc;close all;
x_noisy1=awgn(x,40);
x_noisy2=awgn(x,20);
x_noisy3=awgn(x,15);
powerbw(x_noisy1);
wpass1=1.953*pi*10^(-6);
wpass2=1.984*pi*10^(-6);
wpass3=1.902*pi*10^(-6);

y11=abs(x_noisy1);
y12=lowpass(y11,wpass1,fs);
y21=abs(x_noisy2);
y22=lowpass(y21,wpass2,fs);
y31=abs(x_noisy3);
y32=lowpass(y31,wpass3,fs);

audiowrite('Old_Mc.Donald_noisy_SNR_dB1.wav',y12,fs);
audiowrite('Old_Mc.Donald_noisy_SNR_dB2.wav',y22,fs);
audiowrite('Old_Mc.Donald_noisy_SNR_dB3.wav',y32,fs);

figure;
% frequency response for 40dB
fft_x = fft(y12); %take fourier transform
fft_x_shift = fftshift(fft_x); %shift zero-frequency component to center of spectrum
f = fs/2*linspace(-1,1,L); %the frequency interval for ploting the response, after shifting the fft, from -fs/2 to +fs/2
plot(f, abs(fft_x_shift));
title('frequency response for y12(t)');
xlabel('Frequency (Hz)');
ylabel('magnitude');


