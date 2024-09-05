%% part1: filtering&down-sampling 64channel data, and creating correlation matrix R
clear;clc;close all;
data=cell2mat(struct2cell(load('64channeldata.mat'))); %loading data

%selecting channel1 from data, to show spectrum
temp=data(1,:,:); 
ch1=zeros(size(data,2),size(data,3));
ch1(:,:)=temp;

trial=ch1(:,5).'; %5th trial in channel1

fs=600; %sampling frequency

%plotting spectrum for channel1, trial5
figure('Name','Fourier Transform for channels');
subplot(1,2,1);
plotfft(trial,fs,1);
title('Channel1, trial5');

flt=BPfilt_Q4;
h=flt.Numerator;
y1=conv(trial,h,'same');
subplot(1,2,2);
plotfft(y1,fs,1);
title('Channel1, trial5, after filtering');


excits=size(data,3); %number of excitations
channels=size(data,1); %number of channels
sampsPerTr=size(data,2); %number of samples per trial

data_filt=zeros(channels,sampsPerTr,excits); %filtered data

%applying filter to all epochs
for i=1:excits
    for j=1:channels
        data_filt(j,:,i)=conv(data(j,:,i),h,'same');        
    end
end

%decreasing sampling rate by 4
data_filt_resamp=data_filt(:,1:4:end,:);

% correlation matrix
R=zeros(channels,channels);
%calculating R, under main diagonal
for i=1:channels
    for j=1:i-1
        R(i,j)=corrCal(data_filt_resamp(i,:,1),data_filt_resamp(j,:,1));
        R(j,i)=R(i,j);
    end
    %R equals to 1 on main diag
    R(i,i)=1;
end



%% part2:clustering 64channel data
%we don't clear because R is needed from last section
clc;close all;
y=correlationCluster(R,15);
channel_title = ["AFZ", "FP1", "FP2", "AF3", "AF4", "F7", "F3","FZ","F4","F8",...
 "FC5", "FC1", "FC2", "FC6", "T7", "C3", "CZ", "C4", "T8", "CP5",...
 "CP1", "CP2", "CP6", "P7", "P3", "PZ", "P4","P8","PO3",...
 "PO4", "O1", "O2", "FT10", "AF7", "AF8", "F5","F1","F2","F6",...
 "FT7", "FC3", "FCZ", "FC4", "FT8", "C5", "C1","C2","C6","TP7",...
 "CP3", "CPZ", "CP4", "TP8", "P5", "P1", "P2","P6","PO7","POZ",...
 "PO8", "OZ", "TP9","TP10"];

value = zeros(63,1);
[r, col] = size(y);
for i = 1 : r
    if(isempty(find(y(i,:)==0,1)))
        x=col;
    else
        x=find(y(i,:)==0,1)-1;
    end
    for j = y(i,1:x)
        value(j) = r-i+20;
    end
end
ch_list = cellstr(channel_title);
% run plot_topography for task time based on chanels power
plot_topography(ch_list, value, false,'10-20',true,true,1000)


%% part3:clustering 8 channel data
clear;clc;close all;
filt=load('objFilt2.mat'); %band pass filter
h=filt.Hd.Numerator; %impulse response in time domain
sub=cell2mat(struct2cell(load('Subject1.mat'))); %main matrix
N=size(sub,2); %total size of channels signal
fs=256.41; %sampling frequency

ch1=sub(2,:);
ch2=sub(3,:);
ch3=sub(4,:);
ch4=sub(5,:);
ch5=sub(6,:);
ch6=sub(7,:);
ch7=sub(8,:);
ch8=sub(9,:);

ch1_filt=conv(ch1,h,'same'); %channel1 after filtering by BandPass filter
ch2_filt=conv(ch2,h,'same'); %channel2 after filtering by BandPass filter
ch3_filt=conv(ch3,h,'same'); %channel3 after filtering by BandPass filter
ch4_filt=conv(ch4,h,'same'); %channel4 after filtering by BandPass filter
ch5_filt=conv(ch5,h,'same'); %channel5 after filtering by BandPass filter
ch6_filt=conv(ch6,h,'same'); %channel6 after filtering by BandPass filter
ch7_filt=conv(ch7,h,'same'); %channel7 after filtering by BandPass filter
ch8_filt=conv(ch8,h,'same'); %channel8 after filtering by BandPass filter

channel_filt=[ch1_filt;ch2_filt;ch3_filt;ch4_filt;ch5_filt;ch6_filt;ch7_filt;ch8_filt];
% channel_filt_resamp=channel_filt(:,1:2:end);

z=epoching(channel_filt,0.2,0.8,sub(10,:)); %creating epoch matrix
z_resamp=z(:,1:2:end,:); %down-sampling

channels=size(z_resamp,1);

% correlation matrix
R=zeros(channels,channels);
%calculating R, under main diagonal
for i=1:channels
    for j=1:i-1
        R(i,j)=corrCal(z_resamp(i,:,1),z_resamp(j,:,1));
        R(j,i)=R(i,j);
    end
    %R equals to 1 on main diag
    R(i,i)=1;
end

y=correlationCluster(R,4);

channel_title = ["FZ","CZ","PO7","PO8","P3","OZ","PZ","P4"];

value = zeros(8,1);
[r, col] = size(y);
for i = 1 : r
    if(isempty(find(y(i,:)==0,1)))
        x=col;
    else
        x=find(y(i,:)==0,1)-1;
    end
    for j = y(i,1:x)
        value(j) = r-i+20;
    end
end
ch_list = cellstr(channel_title);
% run plot_topography for task time based on chanels power
plot_topography(cellstr(channel_title), value, false,'10-20',true,true,1000)


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


function r=corrCal(x1,x2) %calculating correlation for 
num=sum(x1.*x2);
den=sqrt(sum(x1.*x1)*sum(x2.*x2));
r=num/den;
end

function y=correlationCluster(corrMat,finalClusts)
channels=size(corrMat,1); %number of channels
cluster=zeros(channels,channels+1); %column1: is there a cluster yet?
cluster(:,1)=ones(channels,1); %at first, we have 63 clusters
cluster(:,2)=(1:channels).'; %at first, each cluster contains exactly one channel

clusters=nnz(cluster(:,1));
while clusters>finalClusts %while number of remaining clusters exceeds the desired number
    clustersCorr=ones(clusters,clusters)*(-1); %initializing clusters correlation matrix
    for i=1:clusters
%         disp('i=');
%         disp(i);
        for j=1:i-1
%             disp('j=');
%             disp(j);
            clustersCorr(i,j)=clusCorr(cluster(i,:),cluster(j,:),corrMat);
        end
    end
    [row,col]=find(clustersCorr==max(max(clustersCorr))); %row&col, for maximum correlation(minimum distance) between clusters
    rowChannels=nnz(cluster(row,2:end));
    colChannels=nnz(cluster(col,2:end));
    combinedChannels=rowChannels+colChannels; %number of channels in the combined cluster
    cluster(row,2:combinedChannels+1)=cat(2,cluster(row,2:rowChannels+1),cluster(col,2:colChannels+1)); %combining two clusters
    cluster(col,:)=[]; %removing empty cluster
    clusters=nnz(cluster(:,1));
end
y=cluster(:,2:end);
end




function D=clusCorr(clus1,clus2,R)
channels1=nnz(clus1)-1; %number of channels belonging to cluster1 (minus 1 because column1 isn't a channel)
channels2=nnz(clus2)-1; %number of channels belonging to cluster2 (minus 1 because column1 isn't a channel)
% corr=-1; %correlation between each channel of clus1, and each channel of clus2
corr=1; %correlation between each channel of clus1, and each channel of clus2
for i=clus1(2:channels1+1)
   for j=clus2(2:channels2+1)
%        disp(j);
       corr=cat(2,corr,R(i,j));
   end    
end
D=min(corr); %maximum correlation equals to minimum distance between channels
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




function Hd = BPfilt_Q4
%BPFILT_Q4 Returns a discrete-time filter object.

% MATLAB Code
% Generated by MATLAB(R) 9.10 and DSP System Toolbox 9.12.
% Generated on: 13-Jul-2022 00:31:43

% Equiripple Bandpass filter designed using the FIRPM function.

% All frequency values are in Hz.
Fs = 600;  % Sampling Frequency

Fstop1 = 0.5;             % First Stopband Frequency
Fpass1 = 1.5;             % First Passband Frequency
Fpass2 = 63;              % Second Passband Frequency
Fstop2 = 64;              % Second Stopband Frequency
Dstop1 = 0.001;           % First Stopband Attenuation
Dpass  = 0.057501127785;  % Passband Ripple
Dstop2 = 0.0001;          % Second Stopband Attenuation
dens   = 20;              % Density Factor

% Calculate the order from the parameters using FIRPMORD.
[N, Fo, Ao, W] = firpmord([Fstop1 Fpass1 Fpass2 Fstop2]/(Fs/2), [0 1 ...
                          0], [Dstop1 Dpass Dstop2]);

% Calculate the coefficients using the FIRPM function.
b  = firpm(N, Fo, Ao, W, {dens});
Hd = dfilt.dffir(b);

% [EOF]
end





% ----------------------------------------------------------------------- %
% Function 'plot_topography' plots a topographical map of the head over the
% desired points given by 'ch_list' and their assigned 'values'.          %
%                                                                         %
%   Input parameters:                                                     %
%       - ch_list:      Channel list in cell-array. Use the string 'all'  %
%                       for displaying all channels available. Note that  %
%                       'Z' indicator should be in upper case.            %
%                       Example: ch_list = {'FP1','FZ','CZ','PZ','OZS'};   %
%       - values:       Numeric vector that contains the values assigned  %
%                       to each channel.                                  %
%       - make_contour: (Optional, default: false) Boolean that controls if
%                       the contour lines should be plotted.              %
%       - system:       (Optional) Measurement system as a string:        %
%           * '10-20':  (Default) 10-20 System, 81 electrodes are available
%           * '10-10':  10-10 system, 47 electrodes are available.        %
%           * 'yokogawa': MEG system by Yokogawa (in testing).            %
%           * table:    Specify a table containing the custom locations,  %
%                       following the structure that is detailed below.   %
%           * path:     Specify a path of a .mat file containing custom   %
%                       locations. The file must have a MATLAB table named%
%                       'locations', with 3 mandatory columns:            %
%                           - labels: contains the name of the electrodes.%
%                           - theta: angle of polar coordinate (in degrees)
%                           - radius: radius of polar coordinate (0-1).   %
%                       Example of 'locations' table:                     %
%                           labels      theta       radius                %
%                       --------------------------------------            %
%                           'Fpz'       0           0.511                 %
%                           'Cz'        90          0                     %
%                           'Oz'        180         0.511                 %
%       - plot_channels:(Optional, default: false) Boolean that controls if
%                       the electrodes should be plotted.                 %
%       - plot_clabels: (Optional, default: false) Boolean that controls if
%                       the text labels of each electrode should be plotted.
%       - INTERP_POINTS:(Optional, default: 1000) No. of interpolation    %
%                       points. The lower N, the lower resolution and     %
%                       faster computation.                               %
%                                                                         %
%   Output variables:                                                     %
%       - h:            Figure handle.                                    %
%                                                                         %
%   Notes:                                                                %
%       - This code was intended to have a small number of input parameters
%       to favor the readability. Therefore, feel free to modify aspects  %
%       such as electrode markers, line widths, colormaps and so on. The  %
%       following lines indicate the key points to modify these aspects:  %
%           * Line 112: Global parameters of number of interpolation      %
%           points or head radius. INTERP_POINTS is set to 1000, if the   %
%           number is increased, the output will be smoother but it will  %
%           take more time. HEAD_RADIUS is fixed to the optimal value in  %
%           order to be anatomically correct if 10-20 system is used.     %
%           * Line 139: Type of interpolation, now 'v4' is used in order to
%           interpolate the surface over a rectangle from -1 to 1.        %
%           * Line 145: Interpolation pcolor.                             %
%           * Line 155: Head plot.                                        %
%           * Line 166: Nose plot.                                        %
%           * Line 183: Ear plots.                                        %
%           * Line 187: Electrode plots.                                  %
% ----------------------------------------------------------------------- %
%   Versions:                                                             %
%       - 1.0:          (18/01/2019) Original script.                     %
%       - 1.1:          (23/01/2019) Yokogawa system added.               %
%       - 1.2:          (04/02/2019) Plot channels added.                 %
%       - 1.3:          (11/09/2019) K-nearest neighbors interpolation.   %
%       - 1.4:          (21/10/2019) Now locations can be directly speci- %
%                       fied as an input table and channel points can be  %
%                       hidden. No. of interp. points may be also passed as
%                       a parameter.                                      %
% ----------------------------------------------------------------------- %
%   Script information:                                                   %
%       - Version:      1.4.                                              %
%       - Author:       V. Martínez-Cagigal                               %
%       - Date:         21/10/2019                                        %
% ----------------------------------------------------------------------- %
%   Example of use:                                                       %
%       plot_topography('all', rand(1,81));                               %
% ----------------------------------------------------------------------- %
function h = plot_topography(ch_list, values, make_contour, system, ...
    plot_channels, plot_clabels, INTERP_POINTS)

    % Error detection
    if nargin < 2, error('[plot_topography] Not enough parameters.');
    else
        if ~iscell(ch_list) && ~ischar(ch_list)
            error('[plot_topography] ch_list must be "all" or a cell array.');
        end
        if ~isnumeric(values)
            error('[plot_topography] values must be a numeric vector.');
        end
    end
    if nargin < 3, make_contour = false;
    else
        if make_contour~=1 && make_contour~=0
            error('[plot_topography] make_contour must be a boolean (true or false).');
        end
    end
    if nargin < 4, system = '10-20';
    else
        if ~ischar(system) && ~istable(system)
            error('[plot_topography] system must be a string or a table.');
        end
    end
    if nargin < 5, plot_channels = true;
    else
        if plot_channels~=1 && plot_channels~=0
            error('[plot_topography] plot_channels must be a boolean (true or false).');
        end
    end
    if nargin < 5, plot_clabels = false;
    else
        if plot_clabels~=1 && plot_clabels~=0
            error('[plot_topography] plot_clabels must be a boolean (true or false).');
        end
    end
    if nargin < 6, INTERP_POINTS = 1000;
    else
        if ~isnumeric(INTERP_POINTS)
            error('[plot_topography] N must be an integer.');
        else
            if mod(INTERP_POINTS,1) ~= 0
                error('[plot_topography] N must be an integer.');
            end
        end
    end
    
    % Loading electrode locations
    if ischar(system)
        switch system
            case '10-20'
                % 10-20 system
                load('Standard_10-20_81ch.mat', 'locations');
            case '10-10'
                % 10-10 system
                load('Standard_10-10_47ch.mat', 'locations');
            case 'yokogawa'
                % Yokogawa MEG system
                load('MEG_Yokogawa-440ag.mat', 'locations');
            otherwise
                % Custom path
                load(system, 'locations');
        end
    else
        % Custom table
        locations = system;
    end
    
    % Finding the desired electrodes
    if ~iscell(ch_list)
        if strcmp(ch_list,'all')
            idx = 1:length(locations.labels);
            if length(values) ~= length(idx)
                error('[plot_topography] There must be a value for each of the %i channels.', length(idx));
            end
        else, error('[plot_topography] ch_list must be "all" or a cell array.');
        end
    else
        if length(values) ~= length(ch_list)
            error('[plot_topography] values must have the same length as ch_list.');
        end
        idx = NaN(length(ch_list),1);
        for ch = 1:length(ch_list)
            if isempty(find(strcmp(locations.labels,ch_list{ch})))
                warning('[plot_topography] Cannot find the %s electrode.',ch_list{ch});
                ch_list{ch} = [];
                values(ch)  = [];
                idx(ch)     = [];
            else
                idx(ch) = find(strcmp(locations.labels,ch_list{ch}));
            end
        end
    end
    values = values(:);
    
    % Global parameters
    %   Note: Head radius should be set as 0.6388 if the 10-20 system is used.
    %   This number was calculated taking into account that the distance from Fpz
    %   to Oz is d=2*0.511. Thus, if the circle head must cross the nasion and
    %   the inion, it should be set at 5d/8 = 0.6388.
    %   Note2: When the number of interpolation points rises, the plots become
    %   smoother and more accurate, however, computational time also rises.
    HEAD_RADIUS     = 5*2*0.511/8;  % 1/2  of the nasion-inion distance
    HEAD_EXTRA      = 1*2*0.511/8;  % 1/10 of the nasion-inion distance
    k = 4;                          % Number of nearest neighbors for interpolation
    
    % Interpolating input data
        % Creating the rectangle grid (-1,1)
        [ch_x, ch_y] = pol2cart((pi/180).*((-1).*locations.theta(idx)+90), ...
                                locations.radius(idx));     % X, Y channel coords
        % Points out of the head to reach more natural interpolation
        r_ext_points = 1.2;
        [add_x, add_y] = pol2cart(0:pi/4:7*pi/4,r_ext_points*ones(1,8));
        linear_grid = linspace(-r_ext_points,r_ext_points,INTERP_POINTS);         % Linear grid (-1,1)
        [interp_x, interp_y] = meshgrid(linear_grid, linear_grid);
        
        % Interpolate and create the mask
        outer_rho = max(locations.radius(idx));
        if outer_rho > HEAD_RADIUS, mask_radius = outer_rho + HEAD_EXTRA;
        else,                       mask_radius = HEAD_RADIUS;
        end
        mask = (sqrt(interp_x.^2 + interp_y.^2) <= mask_radius); 
        add_values = compute_nearest_values([add_x(:), add_y(:)], [ch_x(:), ch_y(:)], values(:), k);
        interp_z = griddata([ch_x(:); add_x(:)], [ch_y(:); add_y(:)], [values; add_values(:)], interp_x, interp_y, 'natural');
        interp_z(mask == 0) = NaN;

        % Plotting the final interpolation
        pcolor(interp_x, interp_y, interp_z);
        shading interp;
        hold on;
        
        % Contour
        if make_contour
            [~, hfigc] = contour(interp_x, interp_y, interp_z); 
            set(hfigc, 'LineWidth',0.75, 'Color', [0.2 0.2 0.2]); 
            hold on;
        end

    % Plotting the head limits as a circle         
    head_rho    = HEAD_RADIUS;                      % Head radius
    if strcmp(system,'yokogawa'), head_rho = 0.45; end
    head_theta  = linspace(0,2*pi,INTERP_POINTS);   % From 0 to 360Âº
    head_x      = head_rho.*cos(head_theta);        % Cartesian X of the head
    head_y      = head_rho.*sin(head_theta);        % Cartesian Y of the head
    plot(head_x, head_y, 'Color', 'k', 'LineWidth',4);
    hold on;

    % Plotting the nose
    nt = 0.15;      % Half-nose width (in percentage of pi/2)
    nr = 0.22;      % Nose length (in radius units)
    nose_rho   = [head_rho, head_rho+head_rho*nr, head_rho];
    nose_theta = [(pi/2)+(nt*pi/2), pi/2, (pi/2)-(nt*pi/2)];
    nose_x     = nose_rho.*cos(nose_theta);
    nose_y     = nose_rho.*sin(nose_theta);
    plot(nose_x, nose_y, 'Color', 'k', 'LineWidth',4);
    hold on;

    % Plotting the ears as ellipses
    ellipse_a = 0.08;                               % Horizontal exentricity
    ellipse_b = 0.16;                               % Vertical exentricity
    ear_angle = 0.9*pi/8;                           % Mask angle
    offset    = 0.05*HEAD_RADIUS;                   % Ear offset
    ear_rho   = @(ear_theta) 1./(sqrt(((cos(ear_theta).^2)./(ellipse_a^2)) ...
        +((sin(ear_theta).^2)./(ellipse_b^2))));    % Ellipse formula in polar coords
    ear_theta_right = linspace(-pi/2-ear_angle,pi/2+ear_angle,INTERP_POINTS);
    ear_theta_left  = linspace(pi/2-ear_angle,3*pi/2+ear_angle,INTERP_POINTS);
    ear_x_right = ear_rho(ear_theta_right).*cos(ear_theta_right);          
    ear_y_right = ear_rho(ear_theta_right).*sin(ear_theta_right); 
    ear_x_left  = ear_rho(ear_theta_left).*cos(ear_theta_left);         
    ear_y_left  = ear_rho(ear_theta_left).*sin(ear_theta_left); 
    plot(ear_x_right+head_rho+offset, ear_y_right, 'Color', 'k', 'LineWidth',4); hold on;
    plot(ear_x_left-head_rho-offset, ear_y_left, 'Color', 'k', 'LineWidth',4); hold on;

    % Plotting the electrodes
    % [ch_x, ch_y] = pol2cart((pi/180).*(locations.theta(idx)+90), locations.radius(idx));
    if plot_channels, he = scatter(ch_x, ch_y, 60,'k', 'LineWidth',1.5); end
    if plot_clabels, text(ch_x, ch_y, ch_list); end
    if strcmp(system,'yokogawa'), delete(he); plot(ch_x, ch_y, '.k'); end
    
    % Last considerations
    max_height = max([max(nose_y), mask_radius]);
    min_height = -mask_radius;
    max_width  = max([max(ear_x_right+head_rho+offset), mask_radius]);
    min_width  = -max_width;
    L = max([min_height, max_height, min_width, max_width]);
    xlim([-L, L]);
    ylim([-L, L]);  
    colorbar;   % Feel free to modify caxis after calling the function
    axis square;
    axis off;
    hold off;
    h = gcf;
end

% This function compute the mean values of the k-nearest neighbors
%   - coor_add:     XY coordinates of the virtual electrodes
%   - coor_neigh:   XY coordinates of the real electrodes
%   - val_neigh:    Values of the real electrodes
%   - k:            Number of neighbors to consider
function add_val = compute_nearest_values(coor_add, coor_neigh, val_neigh, k)
    
    add_val = NaN(size(coor_add,1),1);
    L = length(add_val);
    
    for i = 1:L
        % Distances between the added electrode and the original ones
        target = repmat(coor_add(i,:),size(coor_neigh,1),1);
        d = sqrt(sum((target-coor_neigh).^2,2));
        
        % K-nearest neighbors
        [~, idx] = sort(d,'ascend');
        idx = idx(2:1+k);
        
        % Final value as the mean value of the k-nearest neighbors
        add_val(i) = mean(val_neigh(idx));
    end
    
end


