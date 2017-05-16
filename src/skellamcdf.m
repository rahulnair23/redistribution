function [out,x] = skellamcdf(lambda1,lambda2,range)
%SKELLAMCDF returns the CDF 
%   Detailed explanation goes here
if nargin==3
    x = range;
else
x = -1000:1000;
end

out = zeros(size(x));
for i = 1:length(x)
    out(i) = skellampmf(x(i),lambda1,lambda2);
end
out = cumsum(out);
end %function

    