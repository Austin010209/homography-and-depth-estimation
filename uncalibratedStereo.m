clc;
clear;
close all;
rng(17);


I1all = imread('I21.jpg');     %bookshelf1.jpg
I2all = imread('I22.jpg');     %bookshelf2.jpg
scale = 1;
I1all = imresize(I1all, scale);  %image resize
I2all = imresize(I2all, scale);
I1 = rgb2gray(I1all);            %convert the image to gray
I2 = rgb2gray(I2all);

figure();                        %show original images
imagesc(I1)
colormap(gray);
title('left image');
figure();
imagesc(I2)
colormap(gray);
title('right image');

%use libraries
ptsI1  = detectSURFFeatures(I1);
ptsI2 = detectSURFFeatures(I2);
[featuresI1, validPtsI1] = extractFeatures(I1, ptsI1);
[featuresI2, validPtsI2] = extractFeatures(I2, ptsI2);
indexPairs = matchFeatures(featuresI1, featuresI2);
matchedI1 = validPtsI1(indexPairs(:,1));
matchedI2 = validPtsI2(indexPairs(:,2));

locs1 = matchedI1.Location;
locs2 = matchedI2.Location;




[F,inliersIndex] = estimateFundamentalMatrix(matchedI1,...
    matchedI2,'Method','RANSAC',...
    'NumTrials',50000,'DistanceThreshold',1.35);  %1.35  %best distthresh:1 or 1.5, 1.35

%epipole: Fe1=0 (for image 1, the position of camera 2)
[U,~,V] = svd(F);
e1 = V(:, end);    %null(F);
e1 = unhomo(e1);   %epipole e1
e2 = U(:, end);    %null(F')'
e2 = unhomo(e2);   %epipole e1

inliers1 = locs1(inliersIndex,:);          %leave only inliers
inliers2 = locs2(inliersIndex,:);
select = randi([1 size(inliers1,1)],8,1);  %randomly select 8 points



pts1 = inliers1(select,:);
pts2 = inliers2(select,:);
colors = de2bi((0:7)');
colors(1,:) = [0.9290 0.6940 0.1250];   %(because black is not obvious)
sizerect = size(I1,2)/100;

figure();
imagesc(I1);
colormap(gray);
title(sprintf('epipole in image 1 is (x,y) = (%d,%d)', round(e1(1)), round(e1(2))));
for i=1:8
    pt1 = pts1(i,:);
    ptx = pt1(1);
    pty = pt1(2);
    rectangle('Position',[ptx-sizerect, pty-sizerect, sizerect*2, sizerect*2], 'EdgeColor', colors(i,:), 'LineWidth', 2);
    %find end points of the epipolar line in the image
    [X0, Y0, X1, Y1] = start_end(ptx, pty, e1(1), e1(2), size(I1));
    line([X0 X1], [Y0 Y1], 'Color', colors(i,:), 'LineStyle', '--', 'LineWidth', 1);
end

figure();
imagesc(I2);
colormap(gray);
title(sprintf('epipole in image 2 is (x,y) = (%d,%d)', round(e2(1)), round(e2(2))));
for i=1:8
    pt2 = pts2(i,:);
    ptx = pt2(1);
    pty = pt2(2);
    rectangle('Position',[ptx-sizerect, pty-sizerect, sizerect*2, sizerect*2], 'EdgeColor', colors(i,:), 'LineWidth', 2);
    %find end points of the epipolar line in the image
    [X0, Y0, X1, Y1] = start_end(ptx, pty, e2(1), e2(2), size(I2));
    line([X0 X1], [Y0 Y1], 'Color', colors(i,:), 'LineStyle', '--', 'LineWidth', 1);
end


%next part: estimate disparity
[t1, t2] = estimateUncalibratedRectification(F, inliers1,...
    inliers2, size(I2));
[I1Rect,I2Rect] = rectifyStereoImages(I1,I2,t1,t2);

disparityRange = [-56 56];  %original [-56 56]
disparityMap = disparitySGM(I1Rect,I2Rect,'DisparityRange',disparityRange, 'UniquenessThreshold',0);
figure();
imagesc(disparityMap)  

title('disparity map')
colormap default
colorbar   




function [newpt] = unhomo(pt)
%unhomo
%   calculates the regular coordinate

if(size(pt)==[3,1])
    newpt = [pt(1); pt(2)]/pt(3);
else
    newpt = [pt(1) pt(2)]/pt(3);
end

end



function [X0, Y0, X1, Y1] = start_end(x0, y0, x1, y1, Size)
%start_end 
%   calculate end points of epipolar lines




Size = [Size(2), Size(1)];

m = (y1-y0)/(x1-x0);
if(y1==y0)
    m=0.001;
end
if(x1==x0)
    m=999;
end
b = y0-m*x0;
y_x = @(x) (m*x+b);
x_y = @(y) (y-b)/m;

%x: 0,Size(1)
%y: 0,Size(2)
taken = 0;
y_1 = y_x(0);
y_2 = y_x(Size(1));
if(y_1>=0 && y_1<=Size(2))
    X0 = 0;
    Y0 = y_1;
    taken = 1;
end
if(y_2>=0 && y_2<=Size(2))
    if taken==1
        X1 = Size(1);
        Y1 = y_2;
        taken = 2;
    else
        X0 = Size(1);
        Y0 = y_2;
        taken = 1;
    end
end

if taken==2
    return
end





x_1 = x_y(0);
x_2 = x_y(Size(2));
if(x_1 >= 0 && x_1 <= Size(1))
    if(taken==1)
        X1 = x_1;
        Y1 = 0;
        taken = 2;
    else
        X0 = x_1;
        Y0 = 0;
        taken = 1;
    end
end
if taken==2
    return
end





if(x_2 >= 0 && x_2 <= Size(1))
    if(taken == 1)
        X1 = x_2;
        Y1 = Size(2);
        taken = 2;
    else
        X0 = x_2;
        Y0 = Size(2);
        taken = 1;
    end
end
if taken==2
    return
end

end

function [newpt] = homo(pt)
%homo: calculates homogeneous coordinate

if(size(pt)==[2,1])
    newpt = [pt; 1];
else
    newpt = [pt, 1];
end

end

