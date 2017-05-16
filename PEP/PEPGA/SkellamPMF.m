function [ out ] = SkellamPMF(x,lambda1,lambda2)
%SKELLAMPMF Returns the probability of a skellam distributed variate y
% taking a value of x
%   Computes the prob(y = x) where x has to be integer
%   lambda1 and lambda2 are the Poisson rates of the underlying processes
%   x can be vector or scalar, but the lambda's should be scalar

%make sure its a valid process
if lambda1==0 && lambda2==0
    error('skellamPMF:noProcess','No distribution.');
end

if isscalar(x)
    % error checking
    if floor(x)-x ~= 0
        error('x must be integer');
    end
    if lambda1==0 
            out = - poisspdf(x,lambda2);
    elseif lambda2==0
           out = poisspdf(x,lambda1);
    else
        out = getvalue(x,lambda1,lambda2);
    end
else
    out = zeros(size(x));
    for i = 1:length(out)
        if floor(x(i))-x(i) ~= 0
            error('x must be integer');
        end 
    end
    if lambda1==0
            out = - poisspdf(x,lambda2);
    elseif lambda2==0
            out = poisspdf(x,lambda1);
    else
        for i = 1:length(out)
             out(i) = getvalue(x(i),lambda1,lambda2);
        end
    end
end

end %function

function pout=getvalue(x,lambda1,lambda2)
    
         pout = exp(-lambda1-lambda2)*...
            ((lambda1/lambda2)^(x/2))*...
            besseli(x,2*sqrt(lambda1*lambda2));
end %function

