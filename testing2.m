clearvars
close all
clc

%% image aqc
orig_im1    = (double(imread('Pics for assignment 2/Cig_on_Orange1.JPG'))/255);
orig_im2    = (double(imread('Pics for assignment 2/Cig01.JPG'))/255);
orig_im3    = (double(imread('Pics for assignment 2/Cig03.JPG'))/255);
orig_im4    = (double(imread('Pics for assignment 2/Cig05.JPG'))/255);
orig_im5   = (double(imread('Pics for assignment 2/Cig07.JPG'))/255);
orig_im6   = (double(imread('Pics for assignment 2/Cig08.JPG'))/255);
orig_im7   = (double(imread('Pics for assignment 2/Cig09.JPG'))/255);
orig_im8   = (double(imread('Pics for assignment 2/Cig13.JPG'))/255);

hsv_im1 = rgb2hsv(orig_im1);
hsv_im2 = rgb2hsv(orig_im2);
hsv_im3 = rgb2hsv(orig_im3);
hsv_im4 = rgb2hsv(orig_im4);
hsv_im5 = rgb2hsv(orig_im5);
hsv_im6 = rgb2hsv(orig_im6);
hsv_im7 = rgb2hsv(orig_im7);
hsv_im8 = rgb2hsv(orig_im8);

%% taking R from all images
im1 = orig_im1(:,:,1); % disk = 6/6, sigma = 2 (otherwise only 2 detected)
im2 = orig_im2(:,:,1); % disk = 2/2
im3 = orig_im3(:,:,1); % disk = 6
im1 = orig_im4(:,:,1);
im1 = orig_im5(:,:,1);
im6 = orig_im6(:,:,1); %disk = 4/6, 0.9, 0.86
im7 = orig_im7(:,:,1);
im8 = orig_im8(:,:,1); %same as #6

num_cigs_condition=0;
itter=1;
sigma_delta = 0;
% while(~num_cigs_condition)
if(1)
    %% parameters
    disk_size       = 2*itter;
    disk2_size      = 1;
    sigma           = 2.68 + 0.3*sigma_delta;
    eccen_thresh    = 0.87; %0.94
    solid_thresh    = [0.1,1.1];%[0.78,1];
    major_minor     = [2.25,6.25];
    PA_ratio        = 0.30;
    rect_ratio      = 0.64;

    %% image processing
    se1 = strel('diamond',disk_size);
    se2 = strel('disk',disk2_size);

%     BW1 = edge(im1,'canny');
%     BW1_1 = imdilate(BW1,se1);
%     BW1_2 = bwmorph(BW1_1,'dia');
%     BW1_3 = 1-BW1;
% %     BW1_4 = imdilate(BW1_3,se2);
%     BW1_44 = imfill(BW1_3,'holes');
%     BW1_5 = BW1_3 - bwareaopen(BW1_3,3000);
    
        BW1 = edge(im1,'canny');
    BW1_1 = imdilate(BW1,se1);
    BW1_2 = bwmorph(BW1_1,'dia');
    BW1_3 = 1-BW1_2;
    BW1_4 = imdilate(BW1_3,se2);
    BW1_44 = imfill(BW1_2,'holes');
    BW1_5 = bwareaopen(imerode(BW1_2,se2),45);

    [B,L,N] = bwboundaries(BW1_5);
    stats=  regionprops(L, 'all');
    rect = zeros(4,length(B));

    %% finding cigs
    for k = 1:length(B)
       rect(1,k) = stats(k).Area;
       rect(2,k) = stats(k).Perimeter/stats(k).Area;
       rect(3,k) = stats(k).Eccentricity;
       rect(4,k) = stats(k).Solidity;
       rect(6,k) = stats(k).MajorAxisLength/stats(k).MinorAxisLength;
       rect(7,k) = rect(1,k)/(stats(k).MajorAxisLength*stats(k).MinorAxisLength);
       rect(5,k) = ( (stats(k).Eccentricity > eccen_thresh) &...
           rect(2,k) < PA_ratio &...
           (stats(k).Solidity > solid_thresh(1)) & (stats(k).Solidity < solid_thresh(2)) &...
           (rect(6,k) > major_minor(1)) & ( rect(6,k) < major_minor(2) ) &...
           rect(7,k) > rect_ratio & rect(1,k) < 3000);
       
    end

    %% eliminating false positives
    if(~isempty(rect))
        potential_cigs = find(rect(5,:)>0);
        if(~isempty(potential_cigs))
            A_cig = mean(rect(1,potential_cigs));
            rect(5,potential_cigs) = logical(rect(1,potential_cigs) > 0.5*A_cig) &...
                logical(rect(1,potential_cigs) < 1.75*A_cig);
        end
        num_cigs = sum(rect(5,:));
    end
    
    if itter ==2
        if(num_cigs > 5)
            sigma_delta = sigma_delta + 1;
        else
%             break;
        end
    else
        itter = itter + 1;
    end
    
    num_cigs_condition = (num_cigs > 0) & (num_cigs < 5);
end
find(rect(5,:) > 0)
num_cigs

%% plotting

figure
montage([BW1, BW1_1])
figure
montage([BW1_3, BW1_44])
figure
montage([im1, BW1_5])
hold on

for k = 1:length(B)
   boundary = B{k}; 
    if (find(k==potential_cigs))%(rect(3,k) > eccen_thresh)
       plot(boundary(:,2), boundary(:,1), 'r','LineWidth',2);
       text(stats(k).Centroid(1),stats(k).Centroid(2),...
       num2str(k), 'FontSize',20,'Color','cyan');
    end
   hold on;
end

%% testing plots

% figure;
% subplot(2,2,1), imshow(hsv_im1(:,:,1));
% subplot(2,2,2), imshow(hsv_im1(:,:,2));
% subplot(2,2,3), imshow(hsv_im1(:,:,3));
% subplot(2,2,4), imshow(hsv_im1);
% 
% figure;
% subplot(2,2,1), imshow(hsv_im2(:,:,1));
% subplot(2,2,2), imshow(hsv_im2(:,:,2));
% subplot(2,2,3), imshow(hsv_im2(:,:,3));
% subplot(2,2,4), imshow(hsv_im2);
% 
% figure;
% subplot(2,2,1), imshow(hsv_im3(:,:,1));
% subplot(2,2,2), imshow(hsv_im3(:,:,2));
% subplot(2,2,3), imshow(hsv_im3(:,:,3));
% subplot(2,2,4), imshow(hsv_im3);
% 
% figure;
% subplot(2,2,1), imshow(hsv_im4(:,:,1));
% subplot(2,2,2), imshow(hsv_im4(:,:,2));
% subplot(2,2,3), imshow(hsv_im4(:,:,3));
% subplot(2,2,4), imshow(hsv_im4);
% figure;
% subplot(2,2,1), imshow(hsv_im5(:,:,1));
% subplot(2,2,2), imshow(hsv_im5(:,:,2));
% subplot(2,2,3), imshow(hsv_im5(:,:,3));
% subplot(2,2,4), imshow(hsv_im5);
% 
% figure;
% subplot(2,2,1), imshow(hsv_im6(:,:,1));
% subplot(2,2,2), imshow(hsv_im6(:,:,2));
% subplot(2,2,3), imshow(hsv_im6(:,:,3));
% subplot(2,2,4), imshow(hsv_im6);
% 
% figure;
% subplot(2,2,1), imshow(hsv_im7(:,:,1));
% subplot(2,2,2), imshow(hsv_im7(:,:,2));
% subplot(2,2,3), imshow(hsv_im7(:,:,3));
% subplot(2,2,4), imshow(hsv_im7);
% 
% figure;
% subplot(2,2,1), imshow(hsv_im8(:,:,1));
% subplot(2,2,2), imshow(hsv_im8(:,:,2));
% subplot(2,2,3), imshow(hsv_im8(:,:,3));
% subplot(2,2,4), imshow(hsv_im8);