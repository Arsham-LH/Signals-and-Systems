%% phase 1_2 
clear;
clc;


f = input ('Enter the analog frequencies: '); %% wm
fs = input ('Enter the Sampling frequency: ');  %% ws


Ts = 1/fs ; %% T of sampeling
Ns = 512  ; %% size of signal 
t = [0: Ts :Ts*(Ns-1)]; %% time


x = sin(f*t*2*pi) + cos(pi*f*t); % signal
plotfft(x,fs,1) % using function

%% functions
function plotfft(InputSignal, Fs, flag)
Ns = length(InputSignal) ; 
fourierTransform = fft(InputSignal,Ns).*Fs ;
fourierTransform = fftshift(fourierTransform) ;
if flag == 1
    df = Fs/Ns;
    f = -Fs/2:df:Fs/2-df;
    plot(f,abs(fourierTransform))
else
    Y = [fourierTransform,fourierTransform] ;
    df = Fs/Ns ;
    f = 0:df:Fs-df ;
    f = f.*(2*pi/Fs) ;
    plot(f,abs(Y((Ns/2+1):(3*Ns)/2)))
end
end
