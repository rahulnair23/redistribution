function v = MyJointDistributionFunction(state,params,type)
%MYJOINTDISTRIBUTIONFUNCTION Computes the CDF at state given distribution
%parameters params

%unpack
lambda1=params{1};
lambda2 = params{2};
minState = -50;
v=1;
    for i = 1:length(state)
        if type(i)
            v=v*sum(SkellamPMF(minState:state(i),lambda1(i),lambda2(i)));
        else
            v=v*(1-sum(SkellamPMF(minState:state(i),lambda1(i),lambda2(i))));
        end
    end
end

