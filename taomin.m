function e = taomin(y,x,v,phi)
%UNTITLED5 Summary of this function goes here
%   Detailed explanation goes here

e = mean(sum( ( -y(2:end) + x(1)*y(1:end-1) + x(2)*v(1:end-1) + x(3)*phi(1:end-1) ).^2));

end

