  function result = ConditionalBound(state,Prefs)
    %CONDITIONALBOUND Computes the conditional bound, l(v). Returns true if
    %the bound is met l(v)=v; and false otherwise.
    %   See Beraldi and Rusz..2001

    if strcmpi(Prefs.algorithm,'backwardenumeration')
        result=true;
    else
        if DistributionValue(state,Prefs)>=Prefs.p
            %assume it to be PEP
            result = true;
        else
            %It is not PEP, so reject
            result = false;
            return
        end
    end

    %The actual conditional bound check
    %NOTE: This doesn't actually compute the entire conditional bound l(v),
    %for that see earlier versions. This just checks one dimension lower to
    %ensure that there is no smaller vector that satisfies the PEP
    %condition, in which case the current point is not PEP.
        
        for c = 1:Prefs.nodes
            
            %store current state
            temp = state(c);

            %increment state
            if Prefs.type(c)
                %check for a smaller value
                if state(c)>Prefs.LB(c)
                    state(c) = state(c)-1;
                else
                    continue
                end
                
            else
                
                %check for a larger value
                if state(c)<Prefs.UB(c)
                    state(c) = state(c)+1;
                else
                    %skip this dimension
                    continue
                end
            end

            %do the check
            if DistributionValue(state,Prefs)>=Prefs.p
                result = false;
                break
            end

            %correct state for next iteration
            state(c) = temp;
            
        end %for   
    end % function
