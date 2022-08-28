clc;
clear;
close all;
rng(5);                    %setting random seed
I1 = imread('I1.jpg');    %read the first image
I2 = imread('I2.jpg');    %read the second image
scale = 1;
I1 = imresize(I1, scale);
I2 = imresize(I2, scale);

I1 = rgb2gray(I1); %convert the image to gray
I2 = rgb2gray(I2);
size1 = size(I1);
size2 = size(I2);

%use libraries
ptsI1  = detectSURFFeatures(I1);   %detectSURFFeatures
ptsI2 = detectSURFFeatures(I2);    %detectSURFFeatures
[featuresI1, validPtsI1] = extractFeatures(I1, ptsI1);   %extract features
[featuresI2, validPtsI2] = extractFeatures(I2, ptsI2);   %extract features
indexPairs = matchFeatures(featuresI1, featuresI2);  %get matching features
matchedI1 = validPtsI1(indexPairs(:,1));   %matching features
matchedI2 = validPtsI2(indexPairs(:,2));   %matching features

Locs1 = matchedI1.Location;        %we only care about location
Locs2 = matchedI2.Location;
locs1 = [Locs1(:,2), Locs1(:,1)];  %note that the x and y in Locs1 and 2 are reversed
locs2 = [Locs2(:,2), Locs2(:,1)];  %because we want to work in matlab coordinate

sz_matching = size(locs1,1);       %size of matching pairs
tol1 = sqrt(size1(1)^2 + size1(2)^2)/100;  %tolerance of an exact fit model(in terms of samson distance)
tol2 = tol1/3;                     %tolerance of an least square model
recordC = 0;                       %maximum current size of consensus set
recordH = zeros(3);                %the H corresponding to the largest consensus set


for j = 1:1000
    indices = randi([1 sz_matching],4,1);   %select randomly 4 points
    curH = H4(locs1, locs2, indices);       %calcualte an exact fit
    C1 = get_consensus_set(locs1, locs2, curH, tol1); %get consensus set for this H

    if(size(C1,1) > 0.5*sz_matching)        %if the size of consensus set is large enough
        [C, M1, M2inv] = normalization(C1); %normalize data
        H_norm = H_multiple(C);             %compute least square solution of H
        H = M2inv * H_norm * M1;            %get the H
        C2 = get_consensus_set(locs1, locs2, H, tol2); %get consus set for this new H under a tighter restriction
        if(size(C2, 1) > recordC)           %if size of consensus set is current largest
            recordC = size(C2, 1);          %set this to be the largest
            recordH = H;                    %record this best H
        end
    end
end



%next part  (to stitch the image)
%the extreme value of position of transformed image is originated from corners of image
coords = [1,1; 1,size2(2); size2(1),1; size2(1),size2(2)];
M = round(unhomo(recordH\homo(coords)')); %get transformed coordinates (M is 2*4, first row x, 2nd row y)

%min and max for the transformed image
minx = min(M(1,:));
maxx = max(M(1,:));
miny = min(M(2,:));
maxy = max(M(2,:));
%get size of combined two images by the min and max position of the image
Sizex = max(size1(1), maxx) - min(minx, 1) +1;     %max - min of image, x
Sizey = max(size1(2), maxy) - min(miny, 1) +1;     %max - min of image, y

Fullim = uint8(zeros(Sizex, Sizey, 3));  %construct full image
im1x = max(1, -minx+2);   %x position of (1,1) of image1
im1y = max(1, -miny+2);   %y position of (1,1) of image1
Fullim(im1x : im1x+size1(1)-1, im1y : im1y+size1(2)-1, 1) = I1; %draw I1


for i = 1 : Sizex
    for j = 1 : Sizey  %for every pixel in Fullim
        %H is for original position, so we need to convert back
        oripos = [i - im1x+1; j - im1y+1];
        pos = round(unhomo(recordH*homo(oripos))); %position in image2
        px = pos(1);
        py = pos(2);
        if(px>0 && px<=size2(1) && py>0 && py<=size2(2))   %if image position is valid
            Fullim(i, j, 2) = I2(px,py);  %take the G channel of im2
            Fullim(i, j, 3) = I2(px,py);  %take the B channel of im2
            % Fullim(i, j, 2:3) = reshape(uint8(ones(1,2))*I2(px,py), [1,1,2]);
        end
    end
end
figure();
imagesc(Fullim);   %show result






%functions
function [newpt] = unhomo(pt)
%unhomo Convert pt from homogeneous coordinate to 2D coordinate
%   assume pt is a 3 vector or a matrix with 3 rows or columns

if(size(pt,1) == 3)                       %if pt has 3 rows
    newpt = [pt(1,:); pt(2,:)]./pt(3,:);  %suppress the last row
else
    newpt = [pt(1,:); pt(2,:)]./pt(3,:);  %suppress the last column
end

end



function [temp, M1, M2inv] = normalization(C)
%normalization:   It normalize a matrix with specific form
%   assume C is N*4
%   where first two columns are in one block(coordiantes in im1)
%   and last two columns are in another block(coordinates in im2)

N = size(C,1);       %number of pairs
means = mean(C);
temp = C - means;    %subtract mean first

std1 = norm(temp(:,1:2), 'fro')/sqrt(2*N);  %calculate standard deviation of points in im1 by norm function
std2 = norm(temp(:,3:4), 'fro')/sqrt(2*N);  %calculate standard deviation of points in im2 by norm function
temp(:, 1:2) = temp(:, 1:2)/std1;           %divide two blocks by their stds
temp(:, 3:4) = temp(:, 3:4)/std2;

%now normalization is done; since we need to convert it back sometime, we
%need to construct our transformation in matrix
scale1 = diag([1/std1, 1/std1, 1]);         %rescaling by standard deviation
trans1 = eye(3);
trans1(1,3) = -means(1);
trans1(2,3) = -means(2);                    %translation from the mean of x and y
M1 = scale1*trans1;

scale2 = diag([1/std2, 1/std2, 1]);         %similarly
trans2 = eye(3);
trans2(1,3) = -means(3);
trans2(2,3) = -means(4);
M2 = scale2*trans2;

M2inv = inv(M2);                             %since we only need M2 inverse later
end


function [newpt] = homo(pt)
%homo Convert pt from 2D coordinate to homogeneous coordinate
%   assume pt is a 2 vector or a matrix with 2 rows or columns

if(size(pt,1) == 2)                      %if pt has 2 rows
    newpt = [pt; ones( 1, size(pt,2) )]; %add ones in 3rd row
else
    newpt = [pt, ones( size(pt,1), 1 )]; %add ones in 3rd column
end

end


function [H] = H_multiple(C)
%H_multiple: calcualte the homography matrix given multiple(more, or much more, than 4 pairs)
%   Input: C: N*4 matrix
%   first column: x of matching point in im1
%   second column: y of matching point in im1
%   third column: x of matching point in im2
%   fourth column: y of matching point in im2
%
%   Output: H is a 3*3 homography matrix

N = size(C,1);         %size of matching pairs
M = zeros(2*N,9);      %initialize data matrix
for i = 1:N            %construct data matrix each 2 lines
    x1 = C(i,1);       %x of i th matching point in im1
    y1 = C(i,2);       %y of i th matching point in im1
    x1_1 = C(i,3);     %x of i th matching point in im2
    y1_1 = C(i,4);     %y of i th matching point in im2
    M(2*i-1, :) = [x1, y1, 1, 0, 0, 0, -x1_1*x1, -x1_1*y1, -x1_1];    %put data in
    M(2*i, :) = [0, 0, 0, x1, y1, 1, -y1_1*x1, -y1_1*y1, -y1_1];
end
[~,~,V] = svd(M);
null_space_vector = V(:, end);  %last column of V

H = reshape(null_space_vector, 3, 3)';   %return H

end



function [H] = H4(loc1, loc2, index)
%H4: it calculates the exact homography from 4 matching pairs
%   Input is essentially 4 pairs of matched points
%   in form of all matching points from I1 and I2, and point selected
%   Output is 3*3 homography matrix

M = zeros(8,9);              %initialize data matrix
for i = 1:4                  %construct data matrix every two rows
    pt1 = loc1(index(i), :); %load point 1
    pt2 = loc2(index(i), :); %load point 2
    x1 = pt1(1);             %x of point 1
    y1 = pt1(2);             %y of point 1
    x1_1 = pt2(1);           %x of point 2
    y1_1 = pt2(2);           %y of point 2
    M(2*i-1, :) = [x1, y1, 1, 0, 0, 0, -x1_1*x1, -x1_1*y1, -x1_1];   %put data in
    M(2*i, :) = [0, 0, 0, x1, y1, 1, -y1_1*x1, -y1_1*y1, -y1_1];
end
%find null space of M
[~,~,V] = svd(M);
null_space_vector = V(:, end);   %or we can do: null(M)

H = reshape(null_space_vector, 3, 3)';  %return H
end


function C = get_consensus_set(loc1, loc2, H, tol)
%get_consensus_set:
%   get consensus set of matching points and homography matrix H, unders
%   some tolerance, tol
%   Input:
%   loc1 is the matching points in im1
%   loc2 is the matching points in im2
%   H is the homography matrix computed elsewhere
%   tol is the sampson distance tolerance in transforming points
%  
%   Output: consensus set, a 4 column matrix
%   first two columns are inliers in im1 (x and y)
%   last two columns are inliers in im2 (x and y)

C = zeros(size(loc1,1), 4);  %initialize C
curpo = 1;                   %to add points in order(and potentially last line can use curpos to filter zeros out)
for i = 1:size(loc1,1)
    pt = loc1(i,:);          %extract point from im1
    match_p = unhomo( H * homo(pt') );    %calculate matching point in im2 by H
    error = norm(match_p - loc2(i,:)');   %calculate error
    if(error < tol)                       %if error is small
        C(curpo, :) = [pt, loc2(i,:)];    %add this point and its matching point to consensus set
        curpo = curpo + 1;                %iterate
    end
end
C = C(C(:,1) ~= 0,:);        %suppress zeros tail

end


