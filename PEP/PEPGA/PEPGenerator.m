function out = PEPGenerator(startv,endv,p)
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
if DistributionValue(startv)>p
    %Check to see if the start point is PEP
    if all(ConditionalBound(startv)==startv)
        out = startv;
    else 
        out = [];
    end
    return
end

%The hyperrectangle is "below the p-frontier
if DistributionValue(endv)<=p
    %ignore
    out=[];
    return
end

%The p-frontier passes this zone
if prod(endv-startv)>=minSize
    %recurse
    out=[];
    w = findsplit(startv,endv,p);
    %relay('Split: %2.0f %2.0f\n',w(1),w(2));
    VerticesStart = NDVertices(startv,w);
    VerticesEnd = NDVertices(w,endv);
    for i = 1:length(VerticesStart)
        out = vertcat(out,...
            PEPGenerator(VerticesStart(i,:),VerticesEnd(i,:),p));
    end
else
    %do a local search
    out = PEPEnumerate(startv,endv,p);
end
end

function w=findsplit(s,e,p)
%start from the middle
w=ceil(s+(e-s)/2);
    while DistributionValue(w)<=p
        w=w+1;
    end
end %function