%% part1
clear;clc;close all;
sub=cell2mat(struct2cell(load('Subject1.mat')));
fs=256.41;

figure('Name','Fourier Transform for channels');
subplot(2,4,1);
plotfft(sub(2,:),fs,1);
title('Channel1');

subplot(2,4,2);
plotfft(sub(3,:),fs,1);
title('Channel2');

subplot(2,4,3);
plotfft(sub(4,:),fs,1);
title('Channel3');

subplot(2,4,4);
plotfft(sub(5,:),fs,1);
title('Channel4');

subplot(2,4,5);
plotfft(sub(6,:),fs,1);
title('Channel5');

subplot(2,4,6);
plotfft(sub(7,:),fs,1);
title('Channel6');

subplot(2,4,7);
plotfft(sub(8,:),fs,1);
title('Channel7');

subplot(2,4,8);
plotfft(sub(9,:),fs,1);
title('Channel8');



%% part2: calculating cutoff freq for lowpass filter, using Energy
clear;clc;close all;
sub=cell2mat(struct2cell(load('Subject1.mat')));
fs=256.41;


%calculating cutoff freq for channel1
ch1=sub(2,:);
N1=size(ch1,2); %total size of channel1 signal
cut1=15000; %for more than 95 percent of energy, almost equal to f=59 Hz
y1=fft(ch1);
y1_cut=y1(1,1:cut1);
y1_cut(1,1)=0; %removing DC effect
y1(1,1)=0; %removing DC effect

cutEn1=sum(y1_cut.*conj(y1_cut))/N1; %energy of cut signal
totalEn1=(sum(y1.*conj(y1))/N1)/2; %total energy on positive frequency
r1=cutEn1/totalEn1; %cut energy to total energy ratio

%calculating cutoff freq for channel2
ch2=sub(3,:);
N2=size(ch2,2); %total size of channel1 signal
cut2=15000; %for more than 95 percent of energy, almost equal to f=59 Hz
y2=fft(ch2);
y2_cut=y2(1,1:cut2);
y2_cut(1,1)=0; %removing DC effect
y2(1,1)=0; %removing DC effect


cutEn2=sum(y2_cut.*conj(y2_cut))/N2; %energy of cut signal
totalEn2=(sum(y2.*conj(y2))/N2)/2; %total energy on positive frequency
r2=cutEn2/totalEn2; %cut energy to total energy ratio

%calculating cutoff freq for channel3
ch3=sub(4,:);
N3=size(ch3,2); %total size of channel1 signal
cut3=15000; %for more than 95 percent of energy, almost equal to f=59 Hz
y3=fft(ch3);
y3_cut=y3(1,1:cut3);
y3_cut(1,1)=0; %removing DC effect
y3(1,1)=0; %removing DC effect


cutEn3=sum(y3_cut.*conj(y3_cut))/N3; %energy of cut signal
totalEn3=(sum(y3.*conj(y3))/N3)/2; %total energy on positive frequency
r3=cutEn3/totalEn3; %cut energy to total energy ratio

%calculating cutoff freq for channel4
ch4=sub(5,:);
N4=size(ch4,2); %total size of channel1 signal
cut4=15000; %for more than 95 percent of energy, almost equal to f=59 Hz
y4=fft(ch4);
y4_cut=y4(1,1:cut4);
y4_cut(1,1)=0; %removing DC effect
y4(1,1)=0; %removing DC effect


cutEn4=sum(y4_cut.*conj(y4_cut))/N4; %energy of cut signal
totalEn4=(sum(y4.*conj(y4))/N4)/2; %total energy on positive frequency
r4=cutEn4/totalEn4; %cut energy to total energy ratio

%calculating cutoff freq for channel5
ch5=sub(6,:);
N5=size(ch5,2); %total size of channel1 signal
cut5=15000; %for more than 95 percent of energy, almost equal to f=59 Hz
y5=fft(ch5);
y5_cut=y5(1,1:cut5);
y5_cut(1,1)=0; %removing DC effect
y5(1,1)=0; %removing DC effect


cutEn5=sum(y5_cut.*conj(y5_cut))/N5; %energy of cut signal
totalEn5=(sum(y5.*conj(y5))/N5)/2; %total energy on positive frequency
r5=cutEn5/totalEn5; %cut energy to total energy ratio

%calculating cutoff freq for channel6
ch6=sub(7,:);
N6=size(ch6,2); %total size of channel1 signal
cut6=15000; %for more than 95 percent of energy, almost equal to f=59 Hz
y6=fft(ch6);
y6_cut=y6(1,1:cut6);
y6_cut(1,1)=0; %removing DC effect
y6(1,1)=0; %removing DC effect


cutEn6=sum(y6_cut.*conj(y6_cut))/N6; %energy of cut signal
totalEn6=(sum(y6.*conj(y6))/N6)/2; %total energy on positive frequency
r6=cutEn6/totalEn6; %cut energy to total energy ratio

%calculating cutoff freq for channel7
ch7=sub(8,:);
N7=size(ch7,2); %total size of channel1 signal
cut7=15000; %for more than 95 percent of energy, almost equal to f=59 Hz
y7=fft(ch7);
y7_cut=y7(1,1:cut7);
y7_cut(1,1)=0; %removing DC effect
y7(1,1)=0; %removing DC effect


cutEn7=sum(y7_cut.*conj(y7_cut))/N7; %energy of cut signal
totalEn7=(sum(y7.*conj(y7))/N7)/2; %total energy on positive frequency
r7=cutEn7/totalEn7; %cut energy to total energy ratio

%calculating cutoff freq for channel8
ch8=sub(9,:);
N8=size(ch8,2); %total size of channel1 signal
cut8=15000; %for more than 95 percent of energy, almost equal to f=59 Hz
y8=fft(ch8);
y8_cut=y8(1,1:cut8);
y8_cut(1,1)=0; %removing DC effect
y8(1,1)=0; %removing DC effect

cutEn8=sum(y8_cut.*conj(y8_cut))/N8; %energy of cut signal
totalEn8=(sum(y8.*conj(y8))/N8)/2; %total energy on positive frequency
r8=cutEn8/totalEn8; %cut energy to total energy ratio

%% part3: band pass filter on channels
clear;clc;close all;
filt=load('objFilt2.mat'); %band pass filter
h=filt.Hd.Numerator; %impulse response in time domain
sub=cell2mat(struct2cell(load('Subject1.mat'))); %main matrix
N=size(sub,2); %total size of channels signal
fs=256.41; %sampling frequency

figure('Name','Spectrums after filtering');

subplot(2,4,1);
ch1=sub(2,:); %channel1
ch1_filt=conv(ch1,h); %channel1 after filtering by BandPass filter
plotfft(ch1_filt,fs,1); %spectrum after filtering
title('channel1');

subplot(2,4,2);
ch2=sub(3,:); %channel2
ch2_filt=conv(ch2,h); %channel2 after filtering by BandPass filter
plotfft(ch2_filt,fs,1); %spectrum after filtering
title('channel2');


subplot(2,4,3);
ch3=sub(4,:); %channel3
ch3_filt=conv(ch3,h); %channel3 after filtering by BandPass filter
plotfft(ch3_filt,fs,1); %spectrum after filtering
title('channel3');

subplot(2,4,4);
ch4=sub(5,:); %channel4
ch4_filt=conv(ch4,h); %channel4 after filtering by BandPass filter
plotfft(ch4_filt,fs,1); %spectrum after filtering
title('channel4');

subplot(2,4,5);
ch5=sub(6,:); %channel5
ch5_filt=conv(ch5,h); %channel5 after filtering by BandPass filter
plotfft(ch5_filt,fs,1); %spectrum after filtering
title('channel5');

subplot(2,4,6);
ch6=sub(7,:); %channel6
ch6_filt=conv(ch6,h); %channel6 after filtering by BandPass filter
plotfft(ch6_filt,fs,1); %spectrum after filtering
title('channel6');

subplot(2,4,7);
ch7=sub(8,:); %channel7
ch7_filt=conv(ch7,h); %channel7 after filtering by BandPass filter
plotfft(ch7_filt,fs,1); %spectrum after filtering
title('channel7');

subplot(2,4,8);
ch8=sub(9,:); %channel8
ch8_filt=conv(ch8,h); %channel8 after filtering by BandPass filter
plotfft(ch8_filt,fs,1); %spectrum after filtering
title('channel8');





%% part4: decreasing sampling rate
clear;clc;close all;
sub=cell2mat(struct2cell(load('Subject1.mat')));
channels=sub(2:9,:);
channels_resamp=channels(:,1:2:end);

%% part5:epoching function
clear;clc;close all;
sub=cell2mat(struct2cell(load('Subject1.mat')));
z=epoching(sub(2:end-2,:),0.2,0.8,sub(end-1,:));

%% part6: frequency bands energy
clear;clc;close all;
sub=cell2mat(struct2cell(load('Subject1.mat')));
Fs=256.41;
z=epoching(sub(2:end-2,:),0.2,0.8,sub(end-1,:));
y=freqband(z,Fs); %filtered epoching matrix, with same size
y2=freqband2(z,Fs,0.008); %filtered epoching matrix, with valid size


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


function z=epoching(input,backSamps,forSamps,stimuli)
channels=size(input,1); %number of channels
excits=nnz(stimuli)/4; %number of excitations (each excit equals to 4 non-zero elements)
Ts=0.0039; %sampling period
backIndex=round(backSamps/Ts); %backward indices
forIndex=round(forSamps/Ts); %forward indices

trialTimeLength=forSamps+backSamps; %time length of each trial
trialLength=round(trialTimeLength/Ts)+1; %lenghth of each trial in matrix z


z=zeros(channels,trialLength,excits);

excit=find(stimuli,excits*4);
excit=excit(1:4:end);

for i=1:excits
    z(:,:,i)=input(:,excit(i)-backIndex:excit(i)+forIndex);
end

end

function y=freqband(z,Fs)
% y=zeros(size(z,1),size(z,2),size(z,3));
filt=load('BPfilter.mat'); %band pass filter
h=filt.Num; %impulse response in time domain
excits=size(z,3); %number of excitations
for i=1:excits
   temp=z(:,:,i);
   x=zeros(size(z,1),size(z,2));
   x(:,:)=temp;
   for j=1:8
       y(j,:,i)=conv((x(j,:)).',h,'same');
   end
end
end

function y=freqband2(z,Fs,th)
% y=zeros(size(z,1),size(z,2),size(z,3));
filt=load('BPfilter.mat'); %band pass filter
h=filt.Num; %impulse response in time domain
zeroInd=find(abs(h)<th); %indices of h in which the value is zero (values less than threshold are assumed to be 0)
%removing small values from h
h(1,zeroInd)=0;
[row,col,h]=find(h);

excits=size(z,3); %number of excitations
for i=1:excits
   temp=z(:,:,i);
   x=zeros(size(z,1),size(z,2));
   x(:,:)=temp;
   for j=1:8
       y(j,:,i)=conv((x(j,:)).',h,'valid');
   end
end
end



