clear;
clc;

map = rand(400,400) ;
map = map > 0.5 ;
delay = 0.1 ;
steps = 300 ;
fileName = "MyCGOL.gif" ;
designGif(map,delay,steps,fileName)
%% functions 
function designGif(map,delay,steps,FileName)
step = 1 ;
while step <= steps
    map = OneLevelCGOL(map) ;
    map = (map == 0);
    if step == 1
        imwrite(map,FileName,'gif','LoopCount',inf,'DelayTime',delay);
    else
        imwrite(map,FileName,'gif','WriteMode','append','DelayTime',delay);
    end
    step = step + 1 ;
    map = (map == 0);
end
end



function M = OneLevelCGOL(map)

filter = [1 1 1;1 0.5 1;1 1 1];

neighbor = conv2(map,filter,'same');

dead_cells = (neighbor < 2.5 | neighbor > 3.5);
alived_cells = (neighbor >= 2.5 & neighbor <= 3.5);

map(dead_cells) = 0;
map(alived_cells) = 1;
    
M = map ;
end