%% Q1.1
clear;clc;close all;
img=sobel("edgePic.jpg"); %calling the sobel function
figure('Name','edge-detected image');
imshow(img); %showing final image
title('edge-detected image using sobel');

figure('Name','edge-detected binarized image');
imshow(imbinarize(img)); %showing final binarized image
title('edge-detected binarized image using sobel');
%% Q1.2
clear;clc;
img=kirsch("edgePic.jpg"); %calling the kirsch function
figure('Name','edge-detected image');
imshow(img); %showing final image
title('edge-detected image using kirsch');

figure('Name','edge-detected binarized image');
imshow(imbinarize(img)); %showing final binarized image
title('edge-detected binarized image using kirsch');
%% functions
function edgesImg=sobel(fileName)
    Mx=[-1 0 1; -2 0 2; -1 0 1]; %x kernel
    My=[-1 -2 -1; 0 0 0; 1 2 1]; %y kernel

    img=imread(fileName); %reading the image file
    if ndims(img)==3 %if image is rgb, convert it to grayscale
        img=rgb2gray(img);
    end
    Gx=conv2(im2double(img),Mx); %x gradient
    Gy=conv2(im2double(img),My); %y gradient
    edgesImg=sqrt(Gx.^2+Gy.^2); %derivative matrix
end

function edgesImg=kirsch(fileName)
    g1=[5 5 5;-3 0 -3; -3 -3 -3]; %defining kernels 1 to 8
    g2=[5 5 -3;5 0 -3; -3 -3 -3];
    g3=rot90(g1);
    g4=rot90(g2);
    g5=rot90(g3);
    g6=rot90(g4);
    g7=rot90(g5);
    g8=rot90(g6);
    
    img=imread(fileName); %reading the image file
    if ndims(img)==3 %if image is rgb, convert it to grayscale
        img=rgb2gray(img);
    end
    h1=conv2(im2double(img),g1);
    h2=conv2(im2double(img),g2);
    h3=conv2(im2double(img),g3);
    h4=conv2(im2double(img),g4);
    h5=conv2(im2double(img),g5);
    h6=conv2(im2double(img),g6);
    h7=conv2(im2double(img),g7);
    h8=conv2(im2double(img),g8);

    edgesImg=max(h1,h2);
    edgesImg=max(edgesImg,h3);
    edgesImg=max(edgesImg,h4);
    edgesImg=max(edgesImg,h5);
    edgesImg=max(edgesImg,h6);
    edgesImg=max(edgesImg,h7);
    edgesImg=max(edgesImg,h8);

end
