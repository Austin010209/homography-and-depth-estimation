function [newpt] = homo(pt)
%homo Convert pt from 2D coordinate to homogeneous coordinate
%   assume pt is a 2 vector or a matrix with 2 rows or columns

if(size(pt,1) == 2)                      %if pt has 2 rows
    newpt = [pt; ones( 1, size(pt,2) )]; %add ones in 3rd row
else
    newpt = [pt, ones( size(pt,1), 1 )]; %add ones in 3rd column
end

end
