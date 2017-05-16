 function PEPGenerator(startv,endv,Prefs)
    %PEPGENERATOR Generates a set of p-level efficient points within the
    %hyperrectangle defined by start,end
    %   Recursive
    %   startv: The lower corner of the hyperrec e.g. (0,0,...0)
    %   endv  : The upper corner of e.g. (1,1,..1)
    
   %The hyperrectangle is "above" the p-frontier
     if DistributionValue(startv,Prefs)>Prefs.p;
         return
     end

    %The hyperrectangle is "below the p-frontier
    if DistributionValue(endv,Prefs)<=Prefs.p;
        %ignore
        return
    end

    %The p-frontier passes this zone
    %recurse
    
        w = findsplit(startv,endv,Prefs);
        if ~isempty(w)
            %fprintf(1,'Split\n');
            parfor i = 0:2^Prefs.nodes-1
                %generate the n'th quadrant and recurse
                [l,u]=NQuadrant(startv,w,endv,i);
                PEPGenerator(l,u,Prefs);

                %No can't do in PARFOR
                %check for time
                %if toc(myRunTime)>Prefs.maxtime
                    %we've run out of time
                %    break
                %end
            end
        else
            fprintf(1,'Enumerate\n');
            %do a local search
            switch lower(Prefs.algorithm)
                case 'completeenumeration'
                    CompleteEnumeration(startv,endv);
                case 'backwardenumeration'
                    PEPEnumerate(startv,endv,Prefs);
               case 'prekopa'
                   error('PEPRD:NotTested','Not tested');
                    %Prekopa(startv,endv);
                otherwise
                    error('PEPRD:Prefs:Algorithm','Check Prefs.algorithm.');
            end %switch
        end %if ~isempty(w)
    end %Recursive
