function PEP = Prekopa(low,high)
    %PREKOPA Implements Prek\'{o}pa's PeP enumeration algorithm for the local
    %region defined by low,high
    %relay('Prekopa algorithm: Search space %4.0f\n',prod(high-low+1));
        
    %Call recursive PEP generator
        PEP=[]; %set of local PEP's
        fx=0;  %internal counter showes the number of fixed dimensions
        val = []; %fixed values

        Recursive(fx,val);

        %Add all the PEP's to global list
        for itr=1:size(localPEP,1)
            curr=localPEP(itr,:);
            if ~SHash.containsKey(key(curr))
                   SHash.put(key(curr),curr);
                   fprintf(1,'*');
            end %if
        end % for
        
    function Recursive(fixed,value)
        %relay('*');
        %set iterator
        k_itr=0;
        
        %step 1: Deterimine the bounds
        z = zeros(nodes,1);
        for i = 1:nodes
            if i<=fixed
                z(i)=value(i);
            else
                for j = low(i):high(i)
                    
                    if i==1
                        currstate=[j; high(2:end)];
                    elseif i==nodes
                        currstate=[z(1:end-1); j];
                    else 
                        currstate=[z(1:i-1); j; high(i+1:end)];
                    end
                    
                    if DistributionValue(currstate)>=Prefs.p
                        z(i)=j;
                        break
                    end
                    
                end %for j
            end %if fixed
        end %for i
        
        %Step 2:Add to PEP list
        if ~Dominated(z)
            localPEP(end+1,:)=z(:);
        end
        
        %Step 3: Recurse
        k_itr=k_itr+1;
        FNFV = fixed+1;
        if FNFV<=nodes
            while z(FNFV)+k_itr<UB(FNFV)
                Recursive(FNFV,[value; z(FNFV)+k_itr]);
               k=k+1;
            end %while
        end %if fixed+1<=nodes   
    end %function

    function flag = Dominated(state)
    %DOMINAted Checks to see if the state is dominated by any already in
    %set of PEP's (out)
        flag = false;
        for i = 1:size(localPEP,1)
            if all(localPEP(i,:)<=state')
                flag = true;
                return;
            end
        end
    end %function
    end %function Prekopa
    
    function val = DistributionValue(s,varargin)
    %DISTRIBUTIONVALUE Computes the g(u,v) value by look up. 
    %   g(u,v)=P(u_i<=\xi_i<=v_i, \forall i). Since we assume independence,
    %   g(u,v)=\prod(F(v_i)-F(u_i))
    %   Returns the value of the distribution
    
        if Prefs.precompute
            val = 1;
            half = nodes/2;
            for u=1:half;
                val = val * (marginals(s(u)-LB(u)+1,u)-...
                    marginals(s(half+u)-LB(half+u)+1,half+u));
            end
        else
            val = Prefs.jointDistribution(s,Prefs.parameters);
        end
%      if nargin==1
%         plot3(s(1),s(2),val,'Marker','o','MarkerSize',9,'MarkerFaceColor','r','MarkerEdgeColor','r');
%         drawnow;
%         SearchMatrix(s(1),s(2))=true;
%     end
    end %Computes F(.) by lookup

