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

