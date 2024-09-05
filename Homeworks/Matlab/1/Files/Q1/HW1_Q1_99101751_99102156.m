%% Q1.1
clear;clc;close all

dimension = 500; %dimension of the page
C = zeros(dimension,dimension,3);
initial = randi(dimension-1,1,2); %walker first step
C(initial(1),initial(2),:)=255;
prev = initial; %walker previous step
current = initial; %walker current step
total = 100000; %total steps
dir = 0; %direction of walker
for i=1:total
    prev = current;
    dir = randi(8,1);
    if current(1,1)==1 %if walker reaches top
        dir=4;
    end
    if current(1,1)==dimension-1 %if walker reaches down
        dir=3;
    end
    if current(1,2)==1 %if walker reaches left
        dir=1;
    end
    if current(1,2)==dimension-1 %if walker reaches right
        dir=2;
    end
    
    
    if dir==1 %right
        current(1,2)=prev(1,2)+1;
        C(current(1,1),current(1,2),:)=255;
    elseif dir==2 %left
        current(1,2)=prev(1,2)-1;
        C(current(1,1),current(1,2),:)=255;
    elseif dir==3 %up
        current(1,1) = prev(1,1)-1;
        C(current(1,1),current(1,2),:)=255;
    elseif dir==4 %down
        current(1,1)=prev(1,1)+1;
        C(current(1,1),current(1,2),:)=255;
    elseif dir==5 %up right
        current(1,1)=prev(1,1)-1;
        current(1,2)=prev(1,2)+1;
        C(current(1,1),current(1,2),:)=255;
    elseif dir==6 %up left
        current(1,1)=prev(1,1)-1;
        current(1,2)=prev(1,2)-1;
        C(current(1,1),current(1,2),:)=255;
    elseif dir==7 %down right
        current(1,1)=prev(1,1)+1;
        current(1,2)=prev(1,2)+1;
        C(current(1,1),current(1,2),:)=255;
    elseif dir==8 %down left
        current(1,1)=prev(1,1)+1;
        current(1,2)=prev(1,2)-1;
        C(current(1,1),current(1,2),:)=255;
    end
end
imagesc(C);

%% Q1.2

clear;clc;close all;

dimension = 250;
P = zeros(dimension,dimension,3);
center = dimension/2;
P(center,center,:)=255;
imagesc(P);
total_walkers = 2000;
initial2 = [0,0];
prev2=[0,0];
current2=[0,0];

for i=1:total_walkers
    initial2 = (randi(dimension-6,1,2))+[3,3];
    prev2=initial2;
    current2=initial2;
    exit = false;
    while exit~=true
        prev2 = current2;
        dir = randi(8,1);
        if current2(1,1)<=3 %if walker reaches top
            dir=4;
        elseif current2(1,1)>=(dimension-3) %if walker reaches down
            dir=3;
        end
        if current2(1,2)<=3 %if walker reaches left
            dir=1;
        elseif current2(1,2)>=(dimension-3) %if walker reaches right
            dir=2;
        end
        if dir==1 %right
            current2(1,2)=prev2(1,2)+1;
        elseif dir==2 %left
            current2(1,2)=prev2(1,2)-1;
        elseif dir==3 %up
            current2(1,1) = prev2(1,1)-1;
        elseif dir==4 %down
            current2(1,1)=prev2(1,1)+1;
        elseif dir==5 %up right
            current2(1,1)=prev2(1,1)-1;
            current2(1,2)=prev2(1,2)+1;
        elseif dir==6 %up left
            current2(1,1)=prev2(1,1)-1;
            current2(1,2)=prev2(1,2)-1;
        elseif dir==7 %down right
            current2(1,1)=prev2(1,1)+1;
            current2(1,2)=prev2(1,2)+1;
        elseif dir==8 %down left
            current2(1,1)=prev2(1,1)+1;
            current2(1,2)=prev2(1,2)-1;
        end
        if ((P(current2(1,1)-1,current2(1,2)-1,1)==255) | (P(current2(1,1)-1,current2(1,2),1)==255) | (P(current2(1,1)-1,current2(1,2)+1,1)==255) | (P(current2(1,1),current2(1,2)-1,1)==255) | (P(current2(1,1)-1,current2(1,2)+1,1)==255) | (P(current2(1,1)+1,current2(1,2)-1,1)==255) | (P(current2(1,1)+1,current2(1,2),1)==255) | (P(current2(1,1)+1,current2(1,2)+1,1)==255))
            exit = true;
            P(current2(1,1),current2(1,2),:)=255;
             pause(0.0001);
             imagesc(P);
        end
    end
end


    
%% Q1.3
%khat haye 144,217 (imwrite) baraye ijade file .gif hastan ke chon tu baazi masir ha
%baraye neveshtn permission nadarim, comment shodan. baraye ijade file gif
%mitavan anha ra uncomment kard.
clear;clc;close all;

dimension = 250;
P = zeros(dimension,dimension,3);
center = dimension/2;
P(center,center,:)=255;

h=figure;
axis tight manual;
filename='crystal.gif';

imagesc(P);
frame = getframe(h);
im=frame2im(frame);
[imind,cm]=rgb2ind(im,256);
% imwrite(imind,cm,filename,'gif','Loopcount',Inf,'DelayTime',0.001);

total_walkers = 2000;
initial2 = [0,0];
prev2=[0,0];
current2=[0,0];

sum = (dimension/2)*(dimension/2+1)-(50*51); %sum of weights,minus numbers 1 to 50 
w1=[zeros(1,50),(51:dimension/2)/sum]; %wieghts vector1, with probability 0 for margins
w2=fliplr(w1); %wieghts vector2, with probability 0 for margins
W=[w1,w2]; %final wieghts vector

for i=1:total_walkers
    initial2 = randsample(dimension,2,true,W).';
    prev2=initial2;
    current2=initial2;
    exit = false;
    while exit~=true
        prev2 = current2;
        dirW = zeros(1,8);
        if current2(1,1)<=dimension/2 & current2(1,2)<=dimension/2 %walker in quarter up left
            dirW = [2/13,1/13,1/13,2/13,1/13,1/13,4/13,1/13];
        elseif current2(1,1)<=dimension/2 & current2(1,2)>dimension/2 %walker in quarter up right
            dirW = [1/13,2/13,1/13,2/13,1/13,1/13,1/13,4/13];
        elseif current2(1,1)>dimension/2 & current2(1,2)<=dimension/2 %walker in quarter down left
            dirW = [2/13,1/13,2/13,1/13,4/13,1/13,1/13,1/13];
        elseif current2(1,1)>dimension/2 & current2(1,2)>dimension/2 %walker in quarter down right
            dirW = [1/13,2/13,2/13,1/13,1/13,4/13,1/13,1/13];
        end
        
        dir = randsample(8,1,true,dirW);
            
        if current2(1,1)<=3 %if walker reaches top
            dir=4;
        elseif current2(1,1)>=(dimension-3) %if walker reaches down
            dir=3;
        end
        if current2(1,2)<=3 %if walker reaches left
            dir=1;
        elseif current2(1,2)>=(dimension-3) %if walker reaches right
            dir=2;
        end
        if dir==1 %right
            current2(1,2)=prev2(1,2)+1;
        elseif dir==2 %left
            current2(1,2)=prev2(1,2)-1;
        elseif dir==3 %up
            current2(1,1) = prev2(1,1)-1;
        elseif dir==4 %down
            current2(1,1)=prev2(1,1)+1;
        elseif dir==5 %up right
            current2(1,1)=prev2(1,1)-1;
            current2(1,2)=prev2(1,2)+1;
        elseif dir==6 %up left
            current2(1,1)=prev2(1,1)-1;
            current2(1,2)=prev2(1,2)-1;
        elseif dir==7 %down right
            current2(1,1)=prev2(1,1)+1;
            current2(1,2)=prev2(1,2)+1;
        elseif dir==8 %down left
            current2(1,1)=prev2(1,1)+1;
            current2(1,2)=prev2(1,2)-1;
        end
        if ((P(current2(1,1)-1,current2(1,2)-1,1)==255) | (P(current2(1,1)-1,current2(1,2),1)==255) | (P(current2(1,1)-1,current2(1,2)+1,1)==255) | (P(current2(1,1),current2(1,2)-1,1)==255) | (P(current2(1,1)-1,current2(1,2)+1,1)==255) | (P(current2(1,1)+1,current2(1,2)-1,1)==255) | (P(current2(1,1)+1,current2(1,2),1)==255) | (P(current2(1,1)+1,current2(1,2)+1,1)==255))
            exit = true;
            P(current2(1,1),current2(1,2),:)=255;
             pause(0.0001);
 
             imagesc(P);
             drawnow;
             frame = getframe(h);
             im=frame2im(frame);
             [imind,cm]=rgb2ind(im,256);
%              imwrite(imind,cm,filename,'gif','WriteMode','append','DelayTime',0.001);
        end
    end
end




%% Q1.4

clear;clc;close all;

dimension = 250;
P = zeros(dimension,dimension,3);
center = dimension/2;
P(dimension-3:dimension,center,:)=255;
imagesc(P);
total_walkers = 5000;
initial2 = [0,0];
prev2=[0,0];
current2=[0,0];

%horizontal distribution
sum_h = (dimension/2)*(dimension/2+1)-2*(3+2+1); %sum of weights,minus 1,2,3 
wh1=[zeros(1,50),(51:dimension/2)/sum_h]; %wieghts vector1, with probability 0 for margins
wh2=fliplr(wh1); %wieghts vector2, with probability 0 for margins
W_H=[wh1,wh2]; %final wieghts vector, for horizontal position

%vertical distribution
sum_v = (dimension)*(dimension+1)/2-(1+2+3+248+249+250);
W_V=[zeros(1,3),(4:(dimension-3))/sum_v,zeros(1,3)];

for i=1:total_walkers
    initial2(1,1) = randsample(dimension,1,true,W_V);
    initial2(1,2) = randsample(dimension,1,true,W_H);
    prev2=initial2;
    current2=initial2;
    exit = false;
    while exit~=true
        prev2 = current2;
        dir = randi(8,1);
        if current2(1,1)<=3 %if walker reaches top
            dir=4;
        elseif current2(1,1)>=(dimension-3) %if walker reaches down
            dir=3;
        end
        if current2(1,2)<=3 %if walker reaches left
            dir=1;
        elseif current2(1,2)>=(dimension-3) %if walker reaches right
            dir=2;
        end
        if dir==1 %right
            current2(1,2)=prev2(1,2)+1;
        elseif dir==2 %left
            current2(1,2)=prev2(1,2)-1;
        elseif dir==3 %up
            current2(1,1) = prev2(1,1)-1;
        elseif dir==4 %down
            current2(1,1)=prev2(1,1)+1;
        elseif dir==5 %up right
            current2(1,1)=prev2(1,1)-1;
            current2(1,2)=prev2(1,2)+1;
        elseif dir==6 %up left
            current2(1,1)=prev2(1,1)-1;
            current2(1,2)=prev2(1,2)-1;
        elseif dir==7 %down right
            current2(1,1)=prev2(1,1)+1;
            current2(1,2)=prev2(1,2)+1;
        elseif dir==8 %down left
            current2(1,1)=prev2(1,1)+1;
            current2(1,2)=prev2(1,2)-1;
        end
        if ((P(current2(1,1)-1,current2(1,2)-1,1)==255) | (P(current2(1,1)-1,current2(1,2),1)==255) | (P(current2(1,1)-1,current2(1,2)+1,1)==255) | (P(current2(1,1),current2(1,2)-1,1)==255) | (P(current2(1,1)-1,current2(1,2)+1,1)==255) | (P(current2(1,1)+1,current2(1,2)-1,1)==255) | (P(current2(1,1)+1,current2(1,2),1)==255) | (P(current2(1,1)+1,current2(1,2)+1,1)==255))
            exit = true;
            P(current2(1,1),current2(1,2),:)=255;
             pause(0.0001);
             imagesc(P);
        end
    end
end







 
 