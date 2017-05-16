function w=findsplit(s,e,Prefs)

        if GetCount(s,e)>Prefs.minsize
            %Get all dimensions that cannot be split further
            up=~(s==e|e==s+1|s==e+1);
            w=s;

            %split dimensions that can
            w(up)=floor(s(up)+(e(up)-s(up))/2);

            %check to make sure that split doesn't give the same node
            if all(s==w)||all(e==w)
                w=[];
            end
        else
            %return nothing so the Enumeration algorithm kicks in
            w=[];
        end
    end %function
    
%% Helper func
function c=GetCount(st,en)
        v=abs(en-st)+1;
        v(v==0)=1;
        c=prod(v);
end
