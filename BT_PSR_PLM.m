clear
clc
close all
% READ IN PLM IMAGES + MERGE
A = imread('3L-1-10X-P-0-112-3-1500.tif');
B = imread('3L-1-10X-P-45-67-3-1500.tif');
Z = imadd(A,B);
figure(1);imshow(Z);
title('Combined PLM'); 

% READ IN PRE IMAGEJ MERGED PLM
% X = imread('MAX_0-to-45.tif');
% Y = imread('MAX_0-and-45.tif');
% figure(1);imshow(Z);
% READ IN LIGHT IMAGES
L = imread('3L-1-10X-BF-3.tif'); 
figure(2); imshow(L); 
title('Brightfield');

% PROCESS LIGHT IMAGES
HSV = rgb2hsv(L);
S = HSV(:,:,2); 

% % T H R E S H O L D !! TISSUE ROI
[level EM] = graythresh(L);
BW = im2bw(S,0.15); % Select level or hard code threshold value
invert = ~BW;
tissue_area = sum(invert(:) == 0); %tissue content - count black pixels
figure(3); imshow(invert); title('Tissue Content = Black');

% PROCESS PLM IMAGES
PLM_grey = rgb2gray(Z);
% % T H R E S H O L D !! PLM ROI
BW2 = imbinarize(PLM_grey, 0.25);
collagen_area = sum(BW2(:) == 1);

content = (collagen_area)/tissue_area;
content = content*100;
figure(4); imshow(BW2); 
title(sprintf('Collagen content = %d ', content));