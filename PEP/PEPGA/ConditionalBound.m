function out = ConditionalBound(state)
%CONDITIONALBOUND Computes the conditional bound, to check if (state) is a
%PEP
%   See Prekopa

out = -999*ones(size(state));

    for i = 1:nodes
        for j = LB(i):UB(i)
            temp = state(i);
            state(i) = j;
            if DistributionValue(state)>alpha
                out(i) = j;
                state(i) = temp;
                break;
            end %if
        end %j
    end %i
end
