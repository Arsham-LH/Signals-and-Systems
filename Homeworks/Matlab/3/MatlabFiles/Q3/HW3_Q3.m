%% script1 : delete rest data
% first loading data
clear;
clc;
data = importdata('data.mat') ; %read data
withoutRestData = data(1251:86250,:)  ; % delete rest data
%% script2 : normalize 
normalizedData = zscore(withoutRestData) ;
% tozih daghigh dar gozaresh 
restData = normalizedData(10001:15000,:) + normalizedData(25001:30000,:) + normalizedData(40001:45000,:) + normalizedData(55001:60000,:) + normalizedData(70001:75000,:) ; % 5 time 
taskData = normalizedData(1:10000,:) + normalizedData(15001:25000,:) + normalizedData(30001:40000,:) + normalizedData(45001:55000,:) + normalizedData(60001:70000,:) + normalizedData(75001:85000,:) ; % 6 time
meanRestData = restData./5 ;
meanTaskData = taskData./6 ;
%% script3 : power of signal 
% length of signals 
Nrest = length(meanRestData) ;
Ntask = length(meanTaskData) ;
% calc sigma x[i]^2 
sigmaRest = meanRestData.*meanRestData ;
sigmaTask = meanTaskData.*meanTaskData ;
% calc power
Ptask = (1./Ntask).*sum(sigmaTask,1) ;
Prest = (1./Nrest).*sum(sigmaRest,1) ;
%% topography_rest
channel_title = { 'FP1' ,'FP2' ,'F7' ,'F3' ,'FZ' , ...
'F4' ,'F8' ,'T7' ,'C3' ,'CZ' , ... 
'C4' ,'T8' ,'P7' ,'P3' , 'PZ' ,'P4' ,'P8' ,'O1' ,'O2' };
plot_topography(cellstr(channel_title),Prest,true,'10-20',true,true,1000)
%% topography_task
plot_topography(cellstr(channel_title),Ptask,false,'10-20',true,true,1000)
%% part two 
%% part1
clear;
clc;
eyeData = importdata('eye.mat') ;
%% part2
Size = size(eyeData,2) ; % size of data
fivehundredStep = Size./500 ; % size of steps
XBlink = zeros(1,1); %define vector
for i = 1:fivehundredStep
    [M,I] = max(eyeData((i-1)*500+1:i*500)) ;% define Max
    % define Blink
    if M > 2.5 
        sizeOfXBlink = size(XBlink,2);
        if XBlink(1) == 0
            sizeOfXBlink=0 ;
        end
        XBlink(sizeOfXBlink+1) = I+(i-1)*500;
        YBlink(sizeOfXBlink+1) = M;
    end
end
numberOfBlink = size(XBlink,2) ;
%% part 3
% graphs
% top graph
subplot(2,1,1)
plot(eyeData) 
subplot(2,1,2)
plot(eyeData) 
hold on
stem(XBlink,YBlink,'filled') % higlight Blinks with stem