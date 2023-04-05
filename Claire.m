clc;    % Clear the command window.
clear all;
close all;
workspace;  % Make sure the workspace panel is showing.
format longg;
format compact;
%% Read in PLM images and combine
A = imread('D:\Everything\Tissue Related\Histology\CollContTest\10xPS-2-p1.tif');
B = imread('D:\Everything\Tissue Related\Histology\CollContTest\10xPS-2-p2.tif');
Z = imadd(A,B);
%% Read in brightfild image
L = imread('D:\Everything\Tissue Related\Histology\CollContTest\10xPS-2-unp.tif');

% attempt to threshold out ecm
%Brightfield
HSV = rgb2hsv(L); %RGB iamge to hue, saturation, and value (HSV) values of an HSV image
S = HSV(:,:,2); %values of HSV image - look here to see about thresholding this image
BF_BW = imbinarize(S,0.5); %this threshold can be changed 
figure(); subplot(2,3,1);
imshow(L);
subplot(2,3,2);
imshow(HSV);
subplot(2,3,3);
imshow(BF_BW);

%% Process Images
%Brightfield
HSV = rgb2hsv(L); %RGB iamge to hue, saturation, and value (HSV) values of an HSV image
S = HSV(:,:,2); %values of HSV image
BW = imbinarize(S,0.2);
%PLM
PLM_grey = rgb2gray(Z);
BW2 = imbinarize(PLM_grey, 0.2); %Binarizaion threshold taken from Johnston et al. 2021
figure(1);imshow(BW2);
%% Define regions of interest
figure(1); imshow(L);
roi_Media1 = images.roi.Freehand('Color','c','StripeColor','r');
draw(roi_Media1);
Mask1 = createMask(roi_Media1);
%% Calculate Collagen content
areaMask1 = nnz(Mask1 & BW);
areaCollagen1 = nnz(Mask1 & BW & BW2);%Mask1 & BW & BW2
Collagen_Content1 = areaCollagen1/areaMask1
roi_Media2 = images.roi.Freehand('Color','c','StripeColor','r');
draw(roi_Media2);
Mask2 = createMask(roi_Media2);
areaMask2 = nnz(Mask2 & BW);
areaCollagen2 = nnz(Mask2 & BW & BW2);
Collagen_Content2 = areaCollagen2/areaMask2
roi_Media3 = images.roi.Freehand('Color','c','StripeColor','r');
draw(roi_Media3);
Mask3 = createMask(roi_Media3);
areaMask3 = nnz(Mask3 & BW);
areaCollagen3 = nnz(Mask3 & BW & BW2);
Collagen_Content3 = areaCollagen3/areaMask3
Collagen_Content = (Collagen_Content1 + Collagen_Content2 + Collagen_Content3)/3