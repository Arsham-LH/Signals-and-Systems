%% Q1.1
clear;clc;close all;
x=audioRead('Q1_1.wav')
%% Q1.2
clear;clc;close all;
[allSamps,fs]=audioread('Q1.wav');
M=10;
y=hop(allSamps,M);
% audiowrite('Q1_2.wav',y,fs);
%% Q1.3
clear;clc;close all;
[allSamps,fs]=audioread('C6.wav');
M=16;
y=hop(allSamps,M);
audiowrite('Q1_3.wav',y,fs);
%% Q1.4
clear;clc;close all;
[allSamps,fs]=audioread('sound2.wav');
M1=3;
y1=upSamp(allSamps,M1);
M2=5;
y2=hop(y1,M2);
audiowrite('Q1_4.wav',y2,fs);
%% functions
function x=audioRead(address) %for Q1.1
start=1401;
finish=1410;
[x,fs]=audioread(address,[start,finish]);
% audiowrite('ten samples.wav',x,fs); %for writing new file
end
function y=hop(samples,M)
y=samples(M:M:end);
end
function y=upSamp(samples,M)
y=upsample(samples,M);
end