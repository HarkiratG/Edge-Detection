clearvars
close all
clc

%% parameters
disk_size       = 2;
disk2_size      = 1;

eccen_thresh    = 0.95;
solid_thresh    = 0.70;
major_minor(1,:)= 2.2:0.1:6.2;
PA_ratio(1,:)   = 0.0:0.005:0.28;
rect_ratio      = 0.64;


%% image aqc
orig_im1   = (double(imread('Pics for assignment 2/Cig_on_Orange1.JPG'))/255);
orig_im2   = (double(imread('Pics for assignment 2/Cig01.JPG'))/255);
orig_im1   = (double(imread('Pics for assignment 2/Cig03.JPG'))/255);
orig_im1   = (double(imread('Pics for assignment 2/Cig05.JPG'))/255);
orig_im1   = (double(imread('Pics for assignment 2/Cig07.JPG'))/255);
orig_im6   = (double(imread('Pics for assignment 2/Cig08.JPG'))/255);
orig_im1   = (double(imread('Pics for assignment 2/Cig09.JPG'))/255);
orig_im8   = (double(imread('Pics for assignment 2/Cig13.JPG'))/255);

% % taking R from all images
im1 = orig_im1(:,:,1);

%% image processing
se1 = strel('diamond',disk_size);
se2 = strel('disk',disk2_size);

BW1 = edge(im1,'canny');
BW1_1 = imdilate(BW1,se1);
BW1_2 = imerode(BW1_1,se2);
BW1_3 = bwareaopen(BW1_2 ,45);

[B,L,N] = bwboundaries(BW1_3);
stats=  regionprops(L, 'all');
rect = zeros(8,length(B));

% major_minor(2,:) = (median(major_minor(1,:)) - min(major_minor(1,:)) - ...
% 	abs(major_minor(1,:) - median(major_minor(1,:))))/(median(major_minor(1,:)) - min(major_minor(1,:)));
% 
% PA_ratio(2,:) = (median(PA_ratio(1,:)) - min(PA_ratio(1,:)) - ...
% 	abs(PA_ratio(1,:) - median(PA_ratio(1,:))))/(median(PA_ratio(1,:)) - min(PA_ratio(1,:)));
% median_index = find(PA_ratio(1,:) == median(PA_ratio(1,:)));

%% finding cigs
for k = 1:length(B)
	rect(1,k) = stats(k).Area;
	rect(2,k) = stats(k).Perimeter/stats(k).Area;
	rect(3,k) = stats(k).Eccentricity;
	rect(4,k) = stats(k).Solidity;
	rect(5,k) = rect(1,k)/(stats(k).MajorAxisLength*stats(k).MinorAxisLength);
	rect(6,k) = stats(k).MajorAxisLength/stats(k).MinorAxisLength;
	
	%% weights
% % 	PA_ratio
	
% 	       if(~isempty(PA_ratio(2,round(rect(2,k),2) == PA_ratio(1,:))))
% 	           rect(8,k) = 3*( (PA_ratio(2,round(rect(2,k),2) == PA_ratio(1,:))) );
% 	       else
% 	           rect(8,k) = 0;
% 	       end
	if rect(2,k) < max(PA_ratio(1,:))
		rect(7,k) = 3;
	else
		rect(7,k) = 0;
	end
	
% % 	eccentricity
	if rect(3,k) > eccen_thresh
		rect(7,k) = rect(7,k) + 5;
	else
		rect(7,k) = rect(7,k) + 5*( rect(3,k) );
	end
	
% % 	solidity
	if rect(4,k) > solid_thresh
		rect(7,k) = rect(7,k) + 4*( rect(4,k) );
	end
	
% % 	rect_ratio
% 	       rect(8,k) = rect(8,k) + 2*( (rect(5,k)-0.60)/(1-0.60) );
	if rect(5,k) > 0.6
		rect(7,k) = rect(7,k) + 2*rect(5,k)/0.6;
		
% 		mjor_minor
% 		       if(~isempty(major_minor(2,round(rect(6,k),2) == major_minor(1,:))))
% 		           rect(8,k) = rect(8,k) + 4*( (major_minor(2,round(rect(6,k),1) == major_minor(1,:))) );
% 		       end
		if rect(6,k) > min(major_minor(1,:)) & rect(6,k) < max(major_minor(1,:))
			rect(7,k) = rect(7,k) + 4;
		end
		
	end
% % 	removing minute details
	rect(7,k) = rect(7,k) * (rect(1,k)>60);
	
end
rect(7,:) = rect(7,:) / max(rect(7,:));
rect(8,rect(7,:) > 0.97) = 1;

allcig(1,:) = find(rect(8,:) > 0);
allcig(2:9,:) = rect(:,rect(8,:) > 0);
allcig(:,allcig(2,:) < 0.5*max(allcig(2,:))) = 0;

allcig = allcig(:,allcig(1,:) > 0);

find(rect(8,:) > 0)
sum(rect(8,:))

%% plotting
% figure
% montage([BW1, BW1_1])
% figure
% montage([BW1_2, BW1_3])
% figure
imshow(orig_im1);
hold on

% for k = 1:length(B)
% 	boundary = B{k};
% 	if (rect(8,k) > 0.65)
% 		plot(boundary(:,2), boundary(:,1), 'r','LineWidth',2);
% 		text(stats(k).Centroid(1),stats(k).Centroid(2),...
% 			num2str(k), 'FontSize',10,'Color','cyan');
% 	end
% 	hold on;
% end

dims = size(allcig);

for k = 1:dims(2)
	boundary = B{allcig(1,k)};
	if (allcig(9,k) > 0.65)
		plot(boundary(:,2), boundary(:,1), 'r','LineWidth',2);
	end
	hold on;
end
