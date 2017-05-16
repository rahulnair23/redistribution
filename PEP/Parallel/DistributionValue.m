  function val = DistributionValue(s,Prefs)
    %DISTRIBUTIONVALUE Computes the g(u,v) value by look up. 
    %   g(u,v)=P(u_i<=\xi_i<=v_i, \forall i). Since we assume independence,
    %   g(u,v)=\prod(F(v_i)-F(u_i))
    %   Returns the value of the distribution
    
    val = 1;
    half = length(s)/2;
    for u=1:half;
        val = val * (Prefs.marginals(s(u)-Prefs.LB(u)+1,u)-...
            Prefs.marginals(s(half+u)-Prefs.LB(half+u)+1,half+u));
    end
%     
    end %Computes F(.) by lookup