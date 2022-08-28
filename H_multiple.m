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

