function [newpt] = unhomo(pt)
%unhomo Convert pt from homogeneous coordinate to 2D coordinate
%   assume pt is a 3 vector or a matrix with 3 rows or columns

if(size(pt,1) == 3)                       %if pt has 3 rows
    newpt = [pt(1,:); pt(2,:)]./pt(3,:);  %suppress the last row
else
    newpt = [pt(1,:); pt(2,:)]./pt(3,:);  %suppress the last column
end

end

