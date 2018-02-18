clc
close all
clear all

orig_im1    = rgb2gray(double(imread('Pics for assignment 2/Cig_on_Orange1.JPG'))/255);
orig_im1   = rgb2gray(double(imread('Pics for assignment 2/Cig01.JPG'))/255);
orig_im1    = rgb2gray(double(imread('Pics for assignment 2/Cig03.JPG'))/255);
orig_im1    = (double(imread('Pics for assignment 2/Cig05.JPG'))/255);
orig_im12   = rgb2gray(double(imread('Pics for assignment 2/Cig07.JPG'))/255);
orig_im12   = rgb2gray(double(imread('Pics for assignment 2/Cig08.JPG'))/255);
orig_im12   = rgb2gray(double(imread('Pics for assignment 2/Cig09.JPG'))/255);
orig_im12   = rgb2gray(double(imread('Pics for assignment 2/Cig13.JPG'))/255);

im1_bw = rgb2gray(orig_im1);

im1 = (im1_bw-min(im1_bw(:)));
im1 = im1/max(im1(:));

mask1 = [-1,0,0,0,0,0,0,1];

d_im1_x = conv2(im1,mask1,'same');
d_im1_y = conv2(im1,mask1','same');
g_im1 = d_im1_x.^2+d_im1_y.^2;

% imshow(im1), figure, imshow(d_im1_x), figure, imshow(g_im1)

se = strel('disk', 10);
cc = imclose(g_im1, se);

cc = (cc-min(cc(:)));
cc = cc/max(cc(:));

% dark = find(cc<0.02 | cc>0.8);
% cc(dark) = 0;


% figure, imshow(cc)

% figure;
montage([im1,cc]);

[BW thresh] = edge(im1,'canny');
% [BW thresh] = edge(im1,'log', thresh,3);
imshow(BW)
[H T R] = hough(BW);
P = houghpeaks(H,120);

lines = houghlines(BW,T,R,P);

line_xy = [lines.point1; lines.point2];

figure, imshow(im1), hold on;

for k = 1:length(lines)
   
    xy = [lines(k).point1; lines(k).point2];
    plot(xy(:,1),xy(:,2),'LineWidth',2,'Color','cyan');
    
end





% [eg thresh] = edge(im1,'log');

% eg = edge(im1,'log',thresh,2.7);

% figure, imshow(eg)

% dd_im1 = conv2(d_im1,mask1,'same');
% dd_im2 = conv2(d_im2,mask1,'same');
% dd_im3 = conv2(d_im3,mask1,'same');


% mask2 = [0,-1/4,0;
%         -1/4,1,-1/4;
%         0,-1/4,0];
% ddd_im1 = conv2(im1,mask2,'same');
% ddd_im2 = conv2(im2,mask2,'same');
% ddd_im3 = conv2(im3,mask2,'same');
% 
% ddd_im1 = ddd_im1 + abs(min(ddd_im1(:)));
% ddd_im2 = ddd_im2 + abs(min(ddd_im2(:)));
% ddd_im3 = ddd_im3 + abs(min(ddd_im3(:)));
% 
% cig_dims = [20,20];
% cig_zone = -ones(cig_dims(1),cig_dims(2));
% bg_zone = ones(cig_dims(1)/2,cig_dims(2));

% [x,y] = meshgrid(1:cig_dim(1),1:cig_dim(2));
% 
% mask = [bg_zone; cig_zone; bg_zone];
% 
% edges = conv2(im1,mask);
% ddd_im3 = ddd_im3/max(ddd_im3(:));

% figure;
% imshow(ddd_im3);

% light = find(ddd_im3>0.4);
% ddd_im3(light)=1;


% figure
% imshow(edges);

% figure;
% imshow(ddd_im3);
