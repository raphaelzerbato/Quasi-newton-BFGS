function lol=f_2(varargin)
x = varargin{1};
lol = log(1+3*(x(1,2)-(x(1,1)^3-x(1,1)))^2+(x(1,1)-4/3)^2);
end

