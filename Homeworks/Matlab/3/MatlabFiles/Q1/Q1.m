%% Q1.1
clear;clc;close all;
oldParam=sympref('HeavisideAtOrigin',1);
syms t;
notes=["F5","F5","F5","C5","D5","D5","C5","A5","A5","G5","G5","F5"];
durations=[1,1,1,1,1,1,2,1,1,1,1,4];
% notes=["F5","F5"];
% durations=[1,1];
y=digital_piano(notes,durations);
% totalTime = sum(durations)*dur+(length(durations)-1)*space; %total music duration
audiowrite('Old_Mc.Donald.wav',y,44100);

%% Q1.2
close all;clc;
L=length(y); %total number of samples
fs = 44100; %sampling frequency, equal to fs in digital_piano function
fft_y = fft(y); %take fourier transform
fft_y_shift = fftshift(fft_y); %shift zero-frequency component to center of spectrum
f = fs/2*linspace(-1,1,L); %the frequency interval for ploting the response, after shifting the fft, from -fs/2 to +fs/2
figure;
plot(f, abs(fft_y_shift));
title('frequency response for song signal');
xlabel('Frequency (Hz)');
ylabel('magnitude');

figure;
powerbw(y,fs); %calculating power 3dB bandwidth
%% Q1.7
clc;close all;
[echo_y,fs]=echo('Old_Mc.Donald.wav',3);
audiowrite('Old_Mc.Donald_echo.wav',echo_y,fs);

[reverb_y,fs]=reverb('Old_Mc.Donald.wav',70);
audiowrite('Old_Mc.Donald_reverb_70.wav',reverb_y,fs);

%% Q1.8
clc;close all;

L_echo=length(echo_y); %total number of samples
L_reverb=length(reverb_y); %total number of samples

fs = 44100; %sampling frequency, equal to fs in digital_piano function

fft_y_echo = fft(echo_y); %take fourier transform
fft_y_echo_shift = fftshift(fft_y_echo); %shift zero-frequency component to center of spectrum
fft_y_reverb = fft(reverb_y); %take fourier transform
fft_y_reverb_shift = fftshift(fft_y_reverb); %shift zero-frequency component to center of spectrum

f_echo = fs/2*linspace(-1,1,L_echo); %the frequency interval for ploting the response, after shifting the fft, from -fs/2 to +fs/2
f_reverb = fs/2*linspace(-1,1,L_reverb); %the frequency interval for ploting the response, after shifting the fft, from -fs/2 to +fs/2

figure;
subplot(1,2,1);
plot(f_echo, abs(fft_y_echo_shift));
title('frequency response for echo signal');
xlabel('Frequency (Hz)');
ylabel('magnitude');
subplot(1,2,2);
plot(f_reverb, abs(fft_y_reverb_shift));
title('frequency response for reveberation signal');
xlabel('Frequency (Hz)');
ylabel('magnitude');




%% functions
syms t;
function y=digital_piano(notes,durations)
    Fs=44100; %sampling frequency (number of samples per second)
    Ts=1/Fs;
    dur=0.20; %duration of a single hit (seconds)
    space=0.300; %space duration between two signals(seconds)
    M=2000; %domain of sinosoid signals
    totalTime = sum(durations)*dur+(length(durations))*space; %total music duration
    t=[0:Ts:totalTime]; %time (horizontal axis)
    for i=1:length(notes)
        switch notes(i) %determining frequency for each note in the vector
            case 'C3'
                f0=130.8;
            case 'D3'
                f0=146.8;
            case 'E3'
                f0=164.8;
            case 'F3'
                f0=174.6;
            case 'G3'
                f0=196.0;
            case 'A3'
                f0=220.0;
            case 'B3'
                f0=246.9;
            case 'C4'
                f0=261.6;
            case 'D4'
                f0=293.7;
            case 'E4'
                f0=329.6;
            case 'F4'
                f0=349.2;
            case 'G4'
                f0=392.0;
            case 'A4'
                f0=440.0;
            case 'B4'
                f0=493.9;
            case 'C5'
                f0=523.3;
            case 'D5'
                f0=587.3;
            case 'E5'
                f0=659.3;
            case 'F5'
                f0=698.46;
            case 'G5'
                f0=784.0;
            case 'A5'
                f0=880.0;
            case 'B5'
                f0=987.8;
            case 'C6'
                f0=1047;
            otherwise
                f0=0;
        end
        w0=2*pi*f0;
        signalDur = durations(i)*dur; %current signal duration
        if i==1
            signalStart=0;
        else
            signalStart=sum(durations(1:(i-1)))*dur+(i-1)*space; % where to start current signal(harmonic)
        end
        signalFinish = signalStart+signalDur; % where to finish current signal(harmonic)
        partFinish=signalFinish+space;% where to finish the space after releasing the key
        %creating new signal, and adding to the previuos:
        if i==1
            y=M*(sin(w0*t).*(heaviside(t-signalStart)-heaviside(t-signalFinish))+exp(-0.0075*w0*(t-signalFinish)).*sin(w0*t).*(heaviside(t-signalFinish)-heaviside(t-partFinish)));
        else
            y=y+M*(sin(w0*t).*(heaviside(t-signalStart)-heaviside(t-signalFinish))+exp(-0.0075*w0*(t-signalFinish)).*sin(w0*t).*(heaviside(t-signalFinish)-heaviside(t-partFinish)));
        end
    end
end

% for Q1.5
function [y_eff,fs]=echo(filename,stages)
    [y,fs]=audioread(filename); %getting audio file signal and samoling frequency
    gain=0.4;
    y=y.';
    delayTime=0.20; %delay time between hearing the sound and it's echo (seconds)
    extendedSamplesSize=round(delayTime*fs); %extra samples after echo
    y_extended=cat(2,y,zeros(1,stages*extendedSamplesSize));
    y_eff=y_extended;
    for i=1:stages
        y_shifted=(gain^i)*cat(2,zeros(1,extendedSamplesSize*i), y,zeros(1,extendedSamplesSize*(stages-i)));
        y_eff=y_eff+y_shifted;
    end
end

% for Q1.6
function [y_eff,fs]=reverb(filename,stages)
    [y,fs]=audioread(filename); %getting audio file signal and samoling frequency
    gain=0.60;
    y=y.';
    delayTime=0.090; %delay time between hearing the sound and it's echo (seconds)
    extendedSamplesSize=round(delayTime*fs); %extra samples after echo
    y_extended=cat(2,y,zeros(1,stages*extendedSamplesSize));
    y_eff=y_extended;
    for i=1:stages
        y_eff_shifted=gain*cat(2,zeros(1,extendedSamplesSize), y_eff(1:end-extendedSamplesSize));
        y_eff=y_extended+y_eff_shifted;
    end
end
    
    

    
    
    
    
    
