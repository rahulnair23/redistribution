function  out  = skellaminvCDF(x,y,alpha)
%SKELLAMINVCDF Returns the inverse of CDF 
%   y-cummulative probabilities
%   x-values
%   alpha-desired probability
%
% out-returns the value at alpha value (since it is discrete, returns min)


if alpha<0.5
    %start from lower tail
    for i = 1:length(x)
        if y(i) > alpha
            out = x(i-1);
            break
        end
    end
else
    %start from upper tail
    for i = length(x):-1:1
        if y(i) < alpha
            out = x(i+1);
            break
        end
    end
end
        