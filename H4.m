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

