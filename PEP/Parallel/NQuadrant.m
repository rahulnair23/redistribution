 function [lower,upper] = NQuadrant(start,split,endv,n)
    %NQUADRANT Returns the n'th quadrant of a hyper rectangle defined by the
    %corner points start and endv. Split is the point that induces the
    %quandrants. 
    %   The n'th quadrant is defined the corner points lower and upper

    %v = ('1'==dec2bin(n,length(start))); %too slow
    v=mydec2bin(n,length(start));
    
    %set defaults as corner points
    lower=start;
    upper=split;
    
    %update the vertices
    lower(v==1)=split(v==1);
    upper(v==1)=endv(v==1);
 end
    
 %% Clone of dec2bin
function out= mydec2bin(d,n)
out=uint8(rem(floor(d(:)*pow2(1-n:0)),2));
end