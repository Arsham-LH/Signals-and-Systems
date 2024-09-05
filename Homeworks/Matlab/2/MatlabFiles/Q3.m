%%script1
clear;
clc;
H = tf(1,[1,2,-3]) ;
poles = pole(H) ;
%% script2
clear;
clc;
H = tf(1,[1,2,-3]) ;
poles0 = pole(H) ;
H1 = feedback(H,-1) ;
poles1 = pole(H1) ;
%% script3
clear ;
clc;
poles = zeros(2,11) ;
H = tf(1,[1,2,-3]) ;
for k = -10:2:10
    if k~= 0 
        poles(:,(k+10)./2+1) = pole(feedback(H*k,1)) ;
    end
end
stabels = real(poles) < 0;
stabels = sum(stabels,1) ;
stabels = stabels==2 ;
stabels(6) = 1 ;
k = -10:2:10 ;
stem(k,stabels)
xlabel('k')
ylabel('stabel')
%% script4
clear ;
clc;
H = tf(1,[1,2,-3]) ;
rlocus(H,-100:0.01:100)