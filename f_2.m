function lol=f_2(varargin)
x = varargin{1};
lol = (1-x(1,1))^2+100*(x(2,1)-x(1,1)^2)^2;
end

%(1-x(2,1))^2+100*(x(1,2)-x(1,1)^2)^2