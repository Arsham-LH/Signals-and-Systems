%% script2
clear;
clc;
Ione = imread('pic1.png');% loading image
Itwo = imread('pic2.png');
FTofIone = fft2(Ione) ; % fourier transform
FTofItwo = fft2(Itwo) ; 
IoneNew = ifft2(abs(FTofIone).*exp(angle(FTofItwo).*1i)); % new pics
ItwoNew = ifft2(abs(FTofItwo).*exp(angle(FTofIone).*1i));
figure;% showing pics 
IoneNew = uint8(IoneNew);
ItwoNew = uint8(ItwoNew);
imshow(IoneNew)
figure;
imshow(ItwoNew)
%% script 4.a
clear;
clc;
in = imread('flower.png'); % loading picture
sz = size(in) ; 
k = 4 ; % number of culster
x = randi(sz(2),k,1); 
y = randi(sz(1),k,1);
%output = cluster_1(in,x,y,k); %hadaf soal
output = cluster(in,x,y,k);
image(uint8(output))
%% script 4.b
clear;
clc;
input = imread('coins.jpg');
imshow(input)
input = rgb2gray(input) ;
figure
imshow(otsu(input))
%% script 5
clear;
clc;
load('fmri.mat') % loading image 
image_1 = image(:,:,1) ;
image_2 = image(:,:,2) ;
cor = 70 ;
cor_2 = 0 ;
output = 0;
Ang = 0 ;
XShift = 0 ;
YShift = 0 ;
% find best cor
for ang = 15 : 0.2 : 22
    for xShift = -4:-0.2:-8
        for yShift = -11:-0.2:-15
            shiftImage = imtranslate(image_2,[xShift,yShift]) ; % shift
            rotateImage =  imrotate(shiftImage,ang) ;
            [row , col] = size(rotateImage) ;
            image_1(:,75:col) = 0 ; % tozih dar gozaresh
            image_1(75:row,:) = 0 ;
            cor_2 = corr2(rotateImage, image_1) ;
            if cor_2 > cor % new best form 
                output = zeros(row,col);
                output = rotateImage ;
                cor = cor_2 ; 
                XShift = xShift ;
                YShift = yShift ;
                Ang = ang ;
            end
        end
    end
end
A = imtranslate(image_2,[XShift,YShift]) ; 
output = imrotate(A,Ang) ;
figure; % show
title('image 1')
imshow(image_1)
figure;
title('image 2')
imshow(output)
%% functions
function output_image = cluster(input_image,x,y,k)
check = 0 ;
first = 1 ;
center = ones(1,3) ; % center of cluster 
sz = size(input_image) ; % size of image
clus = ones(sz(1),sz(2),k); % new clustering 
D = ones(sz(1),sz(2),k) ;
centers_1 = ones(k,3) ;
centers_2 = ones(k,3) ;
while check ~= 1
    for i = 1:k
        
        
        % Identify the center of the i th cluster
        if first == 1
           center = double(input_image(y(i),x(i),:)) ;
           centers_1(i,:) = center ;
        else
            center(1) = mean(nonzeros(clus(:,:,i).*double(input_image(:,:,1))),'all');
            center(2) = mean(nonzeros(clus(:,:,i).*double(input_image(:,:,2))),'all');
            center(3) = mean(nonzeros(clus(:,:,i).*double(input_image(:,:,3))),'all');
            centers_1(i,:) = center ;
        end
        %calc A
        page = double(ones(sz(1),sz(2)));
        page = page.*center(1);
        A = double(input_image(:,:,1))- page ;
        A = A.*A;
        
        
        % calc B
        page = double(ones(sz(1),sz(2)));
        page = page.*center(2);
        B = double(input_image(:,:,2))- page ;
        B = B.*B;
        
        % calc C
        page = double(ones(sz(1),sz(2)));
        page = page.*center(3);
        C = double(input_image(:,:,3))- page ;
        C = C.*C;
        
        % pow 2 Size
        D(:,:,i) = A+B+C;
    end
    Min = min(D,[],3) ;
    for j = 1:k
        clus(:,:,j) = (D(:,:,j) == Min); % clustering 
    end
    if(first == 0 && isequal(centers_1,centers_2))
        check = 1;
    end
    first = 0 ;
    centers_2 = centers_1 ;
end
output_image = zeros(sz(1),sz(2),3) ;
for i = 1 : k
    output_image(:,:,1) = output_image(:,:,1) + double(clus(:,:,i)).*centers_1(i,1) ;
    output_image(:,:,2) = output_image(:,:,2) + double(clus(:,:,i)).*centers_1(i,2) ;
    output_image(:,:,3) = output_image(:,:,3) + double(clus(:,:,i)).*centers_1(i,3) ;
end
end
function output_image = cluster_1(input_image,x,y,k)
check = 0 ;
first = 1 ;
center = ones(1,3) ; % center of cluster 
sz = size(input_image) ; % size of image
clus = ones(sz(1),sz(2),k); % new clustering 
D = ones(sz(1),sz(2),k) ;
centers_1 = ones(k,3) ;
centers_2 = ones(k,3) ;
while check ~= 1
    for i = 1:k
        
        
        % Identify the center of the i th cluster
        if first == 1
           center = double(input_image(y(i),x(i),:)) ;
           centers_1(i,:) = center ;
        else
            center(1) = mean(nonzeros(clus(:,:,i).*double(input_image(:,:,1))),'all');
            center(2) = mean(nonzeros(clus(:,:,i).*double(input_image(:,:,2))),'all');
            center(3) = mean(nonzeros(clus(:,:,i).*double(input_image(:,:,3))),'all');
            centers_1(i,:) = center ;
        end
        %calc A
        page = double(ones(sz(1),sz(2)));
        page = page.*center(1);
        A = double(input_image(:,:,1))- page ;
        A = A.*A;
        
        
        % calc B
        page = double(ones(sz(1),sz(2)));
        page = page.*center(2);
        B = double(input_image(:,:,2))- page ;
        B = B.*B;
        
        % calc C
        page = double(ones(sz(1),sz(2)));
        page = page.*center(3);
        C = double(input_image(:,:,3))- page ;
        C = C.*C;
        
        % pow 2 Size
        D(:,:,i) = A+B+C;
    end
    Min = min(D,[],3) ;
    for j = 1:k
        clus(:,:,j) = (D(:,:,j) == Min); % clustering 
    end
    if(first == 0 && isequal(centers_1,centers_2))
        check = 1;
    end
    first = 0 ;
    centers_2 = centers_1 ;
end
output_image = uint8(zeros(sz(1),sz(2))) ;
for i = 1 : k
    output_image = output_image + uint8(clus(:,:,i).*i.*20) ;
end
end
function output_image = otsu(image)
var = ones(123,1) ; % between class variance  
for i = 1 : 123 % check the max of between class variance  
    [mub,muf] = meanCounter(i,image) ; 
    [Wb,Wf] = weightCounter(i,image) ;
    var(i) = Wb.*Wf.*(mub-muf).*(mub-muf) ;
end
Max = max(var);
k =  find(Max == var) ;
output_image = image >= k ;
output_image = output_image.*255 ;
end

function [Mub,Muf] = meanCounter(bound,image)
elementsofb = image < bound ;
elementsoff = image >= bound ;
numberOfelementsb = sum(elementsofb,'all') ;
numberOfelementsf = sum(elementsoff,'all') ;
Mub = sum(double(elementsofb).*double(image),'all')./numberOfelementsb ;
Muf = sum(double(elementsoff).*double(image),'all')./numberOfelementsf ;
end

function [Wb,Wf] = weightCounter(bound,image)
Wb = image < bound ;
Wb = sum(Wb,'all')./numel(Wb) ;
Wf = image >= bound ; 
Wf = sum(Wf,'all')./numel(Wf) ;
end