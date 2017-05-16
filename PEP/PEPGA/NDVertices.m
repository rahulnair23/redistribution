function out = NDVertices(l,u)
%NDVERTICES Generates all the vertices of the hyperrectangle defined by the
%corner vectors l, and u.

if nargin<2
    error('NDVertices:InsufficientInputs','Incorrect number of inputs');
end
if length(l)~=length(u)
    error('NDVertices:IncorrectSize','Inputs of different size.');
end
%Ensure all integer

n=length(l);

v = ('1'==dec2bin(0:2^n-1));
    out =zeros(size(v));
    for i = 1:n
        out(v(:,i)==0,i)=l(i);
        out(v(:,i)==1,i)=u(i);
    end
end %function

