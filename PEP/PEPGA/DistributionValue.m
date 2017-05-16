function val = DistributionValue(state)
%DISTRIBUTIONVALUE Computes the cummulative distribution function for a
%given state.
%   State is a vector of values (v1,v2,....vn)
%   Returns F(x1<=v1,x2<=v2...xn<=vn)

persistent jointCDF;
global nodes;
global LB;
global UB;
    
%val=1;
%     for i =1:nodes
%         val = val * jointCDF{i}(state(i)-LB(i)+1);
%     end %

%test the uniform distribution
val = prod((state-LB+1)./(UB-LB));

end

