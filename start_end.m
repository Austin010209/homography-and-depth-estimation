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

