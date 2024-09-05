%% Q3 part a
clear;clc;close all;
img=imread("white house.jpg");
figure('Name','main image');
imshow(img);
SP_noise=imnoise(img,'salt & pepper',0.1);
figure('Name','salt & pepper noisy image');
imshow(SP_noise);
gauss_noise=imnoise(img,'gaussian',0,0.1);
figure('Name','gaussian noisy image');
imshow(gauss_noise);
poisson_noise=imnoise(img,'poisson');
figure('Name','poisson noisy image');
imshow(poisson_noise);
speckle_noise=imnoise(img,'speckle',0.1);
figure('Name','speckle noisy image');
imshow(speckle_noise);

%% Q3 part c_3 for gaussian
clear;clc;close all;

figure;
img1=gaussianFilter("white house_S&P.jpg",9,1.5);
imshow(img1);
title("gaussian filter effect on S&P noise, size=9, sigma=1.5");

figure;
img2=gaussianFilter("white house_gaussian.jpg",7,2.2);
imshow(img2);
title("gaussian filter effect on gaussian noise, size=7, sigma=2.2");

figure;
img3=gaussianFilter("white house_poisson.jpg",5,1.2);
imshow(img3);
title("gaussian filter effect on poisson noise, size=5, sigma=1.2");

figure;
img4=gaussianFilter("poisson_noisy_pic.png",5,1);
imshow(img4);
title("gaussian filter effect on poisson noise, size=5, sigma=1");

figure;
img5=gaussianFilter("white house_speckle.jpg",7,1.2);
imshow(img5);
title("gaussian filter effect on speckle noise, size=7, sigma=1.2");
%% Q3 part c_3 for median
clear;clc;close all;

%Uncomment each code part for denoising different images
figure;
img1=medianFilter("white house_S&P.jpg",3);
imshow(img1);
title("median filter effect on S&P noise");

%figure;
% img2=medianFilter("white house_gaussian.jpg",5);
% imshow(img2);
% title("median filter effect on gaussian noise");

% figure;
% img3=medianFilter("white house_poisson.jpg",3);
% imshow(img3);
% title("median filter effect on poisson noise");

% figure;
% img4=medianFilter("poisson_noisy_pic.png",5);
% imshow(img4);
% title("median filter effect on poisson noise");

% figure;
% img5=medianFilter("white house_speckle.jpg",5);
% imshow(img5);
% title("median filter effect on speckle noise");

%% Q3 part c_4
clear;clc;close all;
origImg=imread("white house_original.jpg"); %reading original image

SP_noisy=imread("white house_S&P.jpg"); %reading salt & pepper noisy image
SP_noisy=SP_noisy(2:end,:,:); %resizing image

gaussian_noisy=imread("white house_gaussian.jpg"); %reading gaussian noisy image
gaussian_noisy=gaussian_noisy(2:end,:,:);

poisson_noisy=imread("white house_poisson.jpg"); %reading poisson noisy image
poisson_noisy=poisson_noisy(2:end,:,:);

speckle_noisy=imread("white house_speckle.jpg"); %reading speckle noisy image
speckle_noisy=speckle_noisy(2:end,:,:);

SP_denoised_median=imread("S&P noise_median_kernel3.jpg"); %reading denoised salt&pepper (median) image
SP_denoised_gaussian=imread("S&P noise_gaussian_kernel9_sigma1.5.jpg"); %reading denoised salt&pepper (gaussian) image

gaussian_denoised_median=imread("Gaussian noise_median_kernel5.jpg"); %reading denoised gaussian (median) image
gaussian_denoised_gaussian=imread("Gaussian noise_gaussian_kernel7,sigma2.2.jpg"); %reading denoised gaussian (gaussian) image

poisson_denoised_median=imread("poisson noise_median_kernel3.jpg"); %reading denoised poisson (median) image
poisson_denoised_gaussian=imread("poisson noise_gaussian_kernel5_sigma1.2.jpg"); %reading denoised poisson (gaussian) image
poisson_denoised_median=poisson_denoised_median(2:end,:,:); %resizing images
poisson_denoised_gaussian=poisson_denoised_gaussian(2:end,:,:);

speckle_denoised_median=imread("speckle noise_median_kernel5.jpg"); %reading denoised speckle (median) image
speckle_denoised_gaussian=imread("speckle noise_gaussian_kernel7_sigma1.2.jpg"); %reading denoised speckle (gaussian) image

snr_SP_noisy=SNR_cal(origImg,SP_noisy); %calculating SNR using defined function for salt&pepper noisy
snr_gaussian_noisy=SNR_cal(origImg,gaussian_noisy); %calculating SNR using defined function for gaussian noisy
snr_poisson_noisy=SNR_cal(origImg,poisson_noisy); %calculating SNR using defined function for poisson noisy
snr_speckle_noisy=SNR_cal(origImg,speckle_noisy); %calculating SNR using defined function for speckle noisy

snr_SP_median=SNR_cal(origImg,SP_denoised_median); %calculating SNR using defined function for denoised salt&pepper (median)
snr_SP_gaussian=SNR_cal(origImg,SP_denoised_gaussian); %calculating SNR using defined function for denoised salt&pepper (gaussian)
snr_gaus_median=SNR_cal(origImg,gaussian_denoised_median); %calculating SNR using defined function for denoised gaussian (median)
snr_gaus_gaussian=SNR_cal(origImg,gaussian_denoised_gaussian); %calculating SNR using defined function for denoised gaussian (gaussian)
snr_poisson_median=SNR_cal(origImg,poisson_denoised_median); %calculating SNR using defined function for denoised poisson (median)
snr_poisson_gaussian=SNR_cal(origImg,poisson_denoised_gaussian); %calculating SNR using defined function for denoised poisson (gaussian)
snr_speckle_median=SNR_cal(origImg,speckle_denoised_median); %calculating SNR using defined function for denoised speckle (median)
snr_speckle_gaussian=SNR_cal(origImg,speckle_denoised_gaussian); %calculating SNR using defined function for denoised speckle (gaussian)
%% Q3 part c_5
clear;clc;close all;

figure;
img1=meanFilter("white house_S&P.jpg",5);
imshow(img1);
title("mean filter effect on S&P noise, size=5");

figure;
img2=consFilter("white house_S&P_for consFilter.jpg",3);
imshow(img2);
title("conservative filter effect on S&P noise, size=3");



%% functions
function filteredImg=gaussianFilter(fileName,kernelSize,sigma)

    
    img=imread(fileName); %reading file
    img=im2double(img);
    imRows=size(img(:,:,1),1); %number of picture rows
    imCols=size(img(:,:,1),2); %number of picture columns

    % adding zero values for picture edges
    img2=[zeros(floor(kernelSize/2),size(img(:,:,1),2),3);img;zeros(floor(kernelSize/2),size(img(:,:,1),2),3)];
    img3=[zeros(size(img2(:,:,1),1),floor(kernelSize/2),3),img2,zeros(size(img2(:,:,1),1),floor(kernelSize/2),3)];
    img=img3;
    
    kernel=zeros(kernelSize,kernelSize); %creating kernel matrix
    if (mod(kernelSize,2)==1) %if kernel size is odd
        center=floor(kernelSize/2)+1;
        for i=1:kernelSize
            for j=1:kernelSize
                x=i-center;
                y=j-center;
                kernel(i,j)=1/(2*pi*sigma^2)*exp((-x^2-y^2)/(2*sigma^2)); %gaussian distribution
            end
        end
    else %if kernel size is even
        firstQuart=zeros(kernelSize/2,kernelSize/2); %first quarter of kernel matrix
        center=kernelSize/2;
        for i=1:kernelSize/2
            for j=1:kernelSize/2
                x=i-center;
                y=j-center;
                firstQuart(i,j)=1/(2*pi*sigma^2)*exp((-x^2-y^2)/(2*sigma^2));
            end
        end
        kernel(1:center,1:center)=firstQuart; %completing kernel using first quarter
        kernel(center+1:end,1:center)=flipud(firstQuart);
        kernel(1:center,center+1:end)=fliplr(firstQuart);
        kernel(center+1:end,center+1:end)=rot90(firstQuart,2);
    end
    mean_kernel=kernel/sum(sum(kernel)); %dividing by sum, for weighted average
    
    imRed=img(:,:,1); %red matrix of picture
    imGreen=img(:,:,2); %green matrix of picture
    imBlue=img(:,:,3); %blue matrix of picture
    
    redConv=conv2(imRed,mean_kernel); %convolution for each color
    greenConv=conv2(imGreen,mean_kernel);
    blueConv=conv2(imBlue,mean_kernel);
    
    %defining new image
    newIm=zeros(imRows,imCols,3);
    newImRed=imRed;
    newImGreen=imGreen;
    newImBlue=imBlue;
    
    
    for i=1:imRows %for all rows
        for j=1:imCols %for all cols
            newImRed(i,j)=redConv(i,j);
            newImGreen(i,j)=greenConv(i,j);
            newImBlue(i,j)=blueConv(i,j);
        end
    end
   
   start=floor(kernelSize/2)+1;
   newIm(:,:,1)=newImRed(start:end-start+1,start:end-start+1); %filtered picture red matrix
   newIm(:,:,2)=newImGreen(start:end-start+1,start:end-start+1); %filtered picture green matrix
   newIm(:,:,3)=newImBlue(start:end-start+1,start:end-start+1); %filtered picture blue matrix
   filteredImg=newIm;
end

function filteredImg=medianFilter(fileName,kernelSize)

    img=imread(fileName); %reading file
    img=im2double(img);
    
    imRows=size(img(:,:,1),1); %number of picture rows
    imCols=size(img(:,:,1),2); %number of picture columns
    imRed=img(:,:,1); %red matrix of picture
    imGreen=img(:,:,2); %green matrix of picture
    imBlue=img(:,:,3); %blue matrix of picture
    
    newIm=zeros(imRows,imCols,3);
    newImRed=newIm(:,:,1);
    newImGreen=newIm(:,:,2);
    newImBlue=newIm(:,:,3);
    
    center=floor(kernelSize/2)+1; %center of kernel matrix
    for i=center:imRows-center %for all rows
        for j=center:imCols-center %for all cols
            %creating neighbors matrix, and reshaping it to a vactor
            subMatrixRed=reshape(imRed(i-floor(kernelSize/2):i+floor(kernelSize/2),j-floor(kernelSize/2):j+floor(kernelSize/2)),1,[]);
            %calculating the median of neighbors
             newImRed(i,j)=median(subMatrixRed);
            
            %doing the same for green matrix
            subMatrixGreen=reshape(imGreen(i-floor(kernelSize/2):i+floor(kernelSize/2),j-floor(kernelSize/2):j+floor(kernelSize/2)),1,[]);
            newImGreen(i,j)=median(subMatrixGreen);
           
            %doing the same for blue matrix
            subMatrixBlue=reshape(imBlue(i-floor(kernelSize/2):i+floor(kernelSize/2),j-floor(kernelSize/2):j+floor(kernelSize/2)),1,[]);
            newImBlue(i,j)=median(subMatrixBlue);
            
        end
    end
    
   newIm(:,:,1)=newImRed; %filtered picture red matrix
   newIm(:,:,2)=newImGreen; %filtered picture green matrix
   newIm(:,:,3)=newImBlue; %filtered picture blue matrix
   filteredImg=newIm;
end

function snr = SNR_cal(origImg,filtImg)
    logArg=sum(origImg.^2,'all')/sum((origImg-filtImg).^2,'all'); %argument inside log
    snr=10*log10(logArg);
end



function filteredImg=meanFilter(fileName,kernelSize)

    
    img=imread(fileName); %reading file
    img=im2double(img);
    imRows=size(img(:,:,1),1); %number of picture rows
    imCols=size(img(:,:,1),2); %number of picture columns

    % adding zero values for picture edges
    img2=[zeros(floor(kernelSize/2),size(img(:,:,1),2),3);img;zeros(floor(kernelSize/2),size(img(:,:,1),2),3)];
    img3=[zeros(size(img2(:,:,1),1),floor(kernelSize/2),3),img2,zeros(size(img2(:,:,1),1),floor(kernelSize/2),3)];
    img=img3;
    
    kernel=ones(kernelSize,kernelSize); %creating kernel matrix
    mean_kernel=kernel/sum(sum(kernel)); %dividing by sum, for average
    
    imRed=img(:,:,1); %red matrix of picture
    imGreen=img(:,:,2); %green matrix of picture
    imBlue=img(:,:,3); %blue matrix of picture
    
    redConv=conv2(imRed,mean_kernel); %convolution for each color
    greenConv=conv2(imGreen,mean_kernel);
    blueConv=conv2(imBlue,mean_kernel);
    
    %defining new image
    newIm=zeros(imRows,imCols,3);
    newImRed=imRed;
    newImGreen=imGreen;
    newImBlue=imBlue;
    
    
    for i=1:imRows %for all rows
        for j=1:imCols %for all cols
            newImRed(i,j)=redConv(i,j);
            newImGreen(i,j)=greenConv(i,j);
            newImBlue(i,j)=blueConv(i,j);
        end
    end
   
   start=floor(kernelSize/2)+1;
   newIm(:,:,1)=newImRed(start:end-start+1,start:end-start+1); %filtered picture red matrix
   newIm(:,:,2)=newImGreen(start:end-start+1,start:end-start+1); %filtered picture green matrix
   newIm(:,:,3)=newImBlue(start:end-start+1,start:end-start+1); %filtered picture blue matrix
   filteredImg=newIm;
end

function filteredImg=consFilter(fileName,kernelSize)

    img=imread(fileName); %reading file
    img=im2double(img);
    
    imRows=size(img(:,:,1),1); %number of picture rows
    imCols=size(img(:,:,1),2); %number of picture columns
    imRed=img(:,:,1); %red matrix of picture
    imGreen=img(:,:,2); %green matrix of picture
    imBlue=img(:,:,3); %blue matrix of picture
    
    newIm=zeros(imRows,imCols,3);
    newImRed=imRed;
    newImGreen=imGreen;
    newImBlue=imBlue;
    
    center=floor(kernelSize/2)+1; %center of kernel matrix
    for i=center:imRows-center %for all rows
        for j=center:imCols-center %for all cols
            %creating neighbors matrix, and calculating it's min and max
            subMatrixRed=imRed(i-floor(kernelSize/2):i+floor(kernelSize/2),j-floor(kernelSize/2):j+floor(kernelSize/2));
            subCenter=find(subMatrixRed==imRed(i,j));
            subMatrixRed=reshape(subMatrixRed,1,[]);
            subMatrixRed(subCenter)=[]; %removing center element from matrix
            max_subMatrixRed=max(subMatrixRed);
            min_subMatrixRed=min(subMatrixRed);
            %calculating the new value
             if imRed(i,j)<min_subMatrixRed
                 newImRed(i,j)=min_subMatrixRed;
             elseif imRed(i,j)>max_subMatrixRed
                 newImRed(i,j)=max_subMatrixRed;
             end
             
            
            %doing the same for green matrix
            subMatrixGreen=reshape(imGreen(i-floor(kernelSize/2):i+floor(kernelSize/2),j-floor(kernelSize/2):j+floor(kernelSize/2)),1,[]);
            subMatrixGreen(subCenter)=[];
            max_subMatrixGreen=max(subMatrixGreen);
            min_subMatrixGreen=min(subMatrixGreen);
            %calculating the new value
             if imGreen(i,j)<min_subMatrixGreen
                 newImGreen(i,j)=min_subMatrixGreen;
             elseif imGreen(i,j)>max_subMatrixGreen
                 newImGreen(i,j)=max_subMatrixGreen;
             end
           
            %doing the same for blue matrix
            subMatrixBlue=reshape(imBlue(i-floor(kernelSize/2):i+floor(kernelSize/2),j-floor(kernelSize/2):j+floor(kernelSize/2)),1,[]);
            subMatrixBlue(subCenter)=[];
            max_subMatrixBlue=max(subMatrixBlue);
            min_subMatrixBlue=min(subMatrixBlue);
            %calculating the new value
             if imBlue(i,j)<min_subMatrixBlue
                 newImBlue(i,j)=min_subMatrixBlue;
             elseif imBlue(i,j)>max_subMatrixBlue
                 newImBlue(i,j)=max_subMatrixBlue;
             end
            
        end
    end
    
   newIm(:,:,1)=newImRed; %filtered picture red matrix
   newIm(:,:,2)=newImGreen; %filtered picture green matrix
   newIm(:,:,3)=newImBlue; %filtered picture blue matrix
   filteredImg=newIm;
end


