%% collagen content code by RJohnston + editted by BTornifoglio

clc;    % Clear the command window.
clear all;
close all;
workspace;  % Make sure the workspace panel is showing.
format longg;
format compact;

% Define a starting folder. Paste in the working directory. Pop up will
% come up prompting you to change folder if needed, if not, just click
% select folder right away.
start_path = fullfile('D:\Everything\Tissue Related\Histology\plm suss\MULTI ANG');
% Ask user to confirm or change.
topLevelFolder = uigetdir(start_path);
if topLevelFolder == 0
    return;
end

% Get list of all subfolders.
allSubFolders = genpath(topLevelFolder);
% Parse into a cell array.
remain = allSubFolders;
listOfFolderNames = {};
%%
% Loop true folders
while true
    [singleSubFolder, remain] = strtok(remain, ';');
    if isempty(singleSubFolder)
        break;
    end
    listOfFolderNames = [listOfFolderNames singleSubFolder];
end
numberOfFolders = length(listOfFolderNames)

% Process all image files in those folders.
counter = 1;
for k = 1 : numberOfFolders
    
    AA(k)=listOfFolderNames(k);
   
end

i=1;
for j = 1 :numberOfFolders
    
    cd(AA{j})
    disp RUN_CODE
    pwd
    Q=AA{j};
    
    if exist([pwd filesep 'L.tif'], 'file')==2
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
       %% % M A N U A L L Y    C H A N G E    F I L E N A M E S %
        % READ IN TWO POLARISED LIGHT IMAGES AND COMBINE
        A = imread('55-CC1-10X-P-0-112-1.tif');
        B = imread('55-CC1-10X-P-5-107-1.tif');
        two_angles = imadd(A,B);
        C = imread('55-CC1-10X-P-10-100-1.tif');
        three_angles = imadd(two_angles,C);
        D = imread('55-CC1-10X-P-15-95-1.tif');
        four_angles = imadd(three_angles,D);
        E = imread('55-CC1-10X-P-20-89-1.tif');
        five_angles = imadd(four_angles,E);
        F = imread('55-CC1-10X-P-25-86-1.tif');
        six_angles = imadd(five_angles,F);
        G = imread('55-CC1-10X-P-30-82-1.tif');
        seven_angles = imadd(six_angles,G);
        H = imread('55-CC1-10X-P-35-76-1.tif');
        eight_angles = imadd(seven_angles,H);
        I = imread('55-CC1-10X-P-40-70-1.tif');
        nine_angles = imadd(eight_angles,I);
        J = imread('55-CC1-10X-P-45-65-1.tif');
        ten_angles = imadd(nine_angles,J);
        Z = ten_angles;
        
        subplot(3,3,1); imshow(two_angles);
        subplot(3,3,2); imshow(three_angles);
        subplot(3,3,3); imshow(four_angles);
        subplot(3,3,4); imshow(five_angles);
        subplot(3,3,5); imshow(six_angles);
        subplot(3,3,6); imshow(seven_angles);
        subplot(3,3,7); imshow(eight_angles);
        subplot(3,3,8); imshow(nine_angles);
        subplot(3,3,9); imshow(ten_angles);
        
        figure(1);imshow(ten_angles);
        title('Combined PLM'); 
   %     saveas(gcf,'CombinedPLM.tif'); %save merged PLM
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % M A N U A L L Y    C H A N G E    F I L E N A M E  %  
        % READ IN LIGHT IMAGES
        L = imread('55-CC1-10X-BF-1.tif'); 
        figure(2); imshow(L); 
        title('Brightfield'); 
    %    saveas(gcf,'BF.tif');
        
        % PROCESS LIGHT IMAGES
        HSV = rgb2hsv(L);
       % figure(5); imshow(HSV); 
       % title('HSV - Tissue Content'); 
        %saveas(gcf,'HSV.tif');
        
        S=HSV(:,:,2);
       % figure(6); imshow(S); 
       % title('Saturation Channel to determine tissue content'); 
        %saveas(gcf,'SaturationChannel.tif');
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%T H R E S H O L D !!
        [level EM] = graythresh(L);
        BW = im2bw(S,0.2); % Select level or hard code threshold value
        invert = ~BW;
        tissue_area = sum(invert(:) == 0); %tissue content - count black pixels
        figure(3); imshow(invert); title('Tissue Content = Black');
    %    saveas(gcf,'TissueROI.tif');   
        
        % PROCESS PLM IMAGES
        PLM_grey = rgb2gray(Z);
        %figure(3); imshow(PLM_grey); title('Greyscale PLM');
        %saveas(gcf,'Greyscale.tif');
         
        % Adaptive method
        % BW2 = imbinarize(PLM_grey, 'adaptive','Sensitivity',0.5); adds
        % speckles
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%T H R E S H O L D !!
        BW2 = imbinarize(PLM_grey, 0.2);
        %collagen content - count white pixels
        collagen_area = sum(BW2(:) == 1);
        %collagen content = white
       % figure(4); imshow(BW2); title('Binarise PLM - collagen area');
       % saveas(gcf,'CollagenArea.tif');      

       
        %Clean up stray blobs in image smaller than 80 pixels
      %  BW2_clean = bwareaopen(BW2, 15);
      %  clean_collagen_area = sum(BW2_clean(:) == 1);

        % CALCULATE PERCENT COLLAGEN
        %%content clean misses a lot...use collagen area
        content = (collagen_area)/tissue_area;
        content = content*100;
       % content_clean = (clean_collagen_area)/tissue_area
      %  content_clean = content_clean*100;
        figure(4); imshow(BW2); 
        title(sprintf('Collagen content = %d ', content));
     %   saveas(gcf,'Coll%.tif');
       % figure(16); imshow(BW2_clean); 
       % title(sprintf('Collagen content clean = %d ', content_clean)),saveas(gcf,'CollagenContentClean_Ex.tif');
        i=i+1;
       
    end
end

% 
% figure(5)
% subplot(2,2,1)
% imshow(L); title('BF'); 
% subplot(2,2,2)
% imshow(Z);  title('Combined PLM'); 
% subplot(2,2,3)
% imshow(invert); title('Tissue ROI - threshold (0.25)');
% subplot(2,2,4)
% imshow(BW2); title(sprintf('Collagen content (thres(0.3)) = %d ', content));
% saveas(gcf,'ROI1.tif');