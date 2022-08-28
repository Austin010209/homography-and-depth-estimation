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

