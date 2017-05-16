function S=Main(s,e,l1,l2,p)
%MAIN Base of operations for PEP generation

if nargin==0
    %test prob
    s=ones(1,28);
    e=5*ones(1,28);
    l1 = [1 2 1 2 1 2 1 2 1 1 1 2 2 2];
    l2 = [2 1 2 2 2 1 2 2 1 2 2 2 2 2];
    l1 = [l1 l1];
    l2=[l2 l2];
    p=0.9;
end

alpha = p;
LB=s;
UB=e;
lambda1 = l1;
lambda2 = l2;
nodes=length(LB);

% Correct the bounds 
% while slow_DistVal(LB)<=0.001
%     LB=LB+1;
% end
% while slow_DistVal(UB)<=0.99
%     UB=UB+1;
% end

%Create a look up table
jointCDF={};
InitDistribution;
tic
%Call recursive PEP generator
S = PEPGenerator(LB,UB);
toc
    function v=slow_DistVal(state)
        F=zeros(size(state));
        for i = 1:nodes
            vals = -50:state;
            F(i) = sum(SkellamPMF(vals,lambda1(i),lambda2(i)));
        end
        v=prod(F);
    end %function
        
    function InitDistribution
        %INITDISTRIBUTION:Precomputes the marginals for each state
        smallest=-100;
        jointCDF = cell(nodes,1);
        for i =1:nodes
            temp=zeros(UB(i)-LB(i),1);
            k=1;
            for j = LB(i):UB(i)
                if k==1
                    z=smallest;
                    while z <=LB(i)
                        temp(k)=temp(k)+SkellamPMF(z,lambda1(i),lambda2(i));
                        z=z+1;
                    end
                else
                    temp(k) = temp(k-1)+SkellamPMF(j,lambda1(i),lambda2(i));
                end %if
                k=k+1;
            end %for j
            jointCDF{i}=temp;
        end %i
     end %function InitDistr
 
    function val = DistributionValue(state)
    %DISTRIBUTIONVALUE Computes the cummulative distribution function for a
    %given state.
    %   State is a vector of values (v1,v2,....vn)
    %   Returns F(x1<=v1,x2<=v2...xn<=vn)

        val=1;
        for i =1:nodes
            u=jointCDF{i}(state(i)-LB(i)+1);
            val = val * u;
        end %

    %test the uniform distribution
    %val = prod((state-LB+1)./(UB-LB));

    end %Computes F(.)

    function out = PEPGenerator(startv,endv)
    %PEPGENERATOR Generates a set of p-level efficient points within the
    %hyperrectangle defined by start,end
    %   Recursive
    %   startv: The lower corner of the hyperrec e.g. (0,0,...0)
    %   endv  : The upper corner of e.g. (1,1,..1)
    %
    %
    %relay('Corners: (%2.0f, %2.0f) (%2.0f, %2.0f)\n',startv(1),startv(2),endv(1),endv(2));

    minSize = 1000; %Minimum size of the neighborhood till enumeration is used

    %The hyperrectangle is "above" the p-frontier
    if DistributionValue(startv)>alpha
        %Check to see if the start point is PEP
        if all(ConditionalBound(startv)==startv)
            out = startv;
        else 
            out = [];
        end
        return
    end

    %The hyperrectangle is "below the p-frontier
    if DistributionValue(endv)<=alpha
        %ignore
        out=[];
        return
    end

    %The p-frontier passes this zone
    if prod(endv-startv)>=minSize
        %recurse
        out=[];
        w = findsplit(startv,endv);
        %relay('Split: %2.0f %2.0f\n',w(1),w(2));
        VerticesStart = NDVertices(startv,w);
        VerticesEnd = NDVertices(w,endv);
        for i = 1:length(VerticesStart)
            out = vertcat(out,...
                PEPGenerator(VerticesStart(i,:),VerticesEnd(i,:)));
        end
    else
        %do a local search
        out = PEPEnumerate(startv,endv);
    end
    end %Recursive

    function w=findsplit(s,e)
    %start from the middle
    w=floor(s+(e-s)/2);
%         while DistributionValue(w)<=alpha
%             w=min(w+1,UB);
%         end
%         w=w-1;
    end %function

    function o = PEPEnumerate(low,high)
    %PEPENUMERATE Summary of this function goes here
    %   Detailed explanation goes here

    o=[];
    end

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

end

