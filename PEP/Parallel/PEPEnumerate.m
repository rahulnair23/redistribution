 function PEPEnumerate(low,high,Prefs)
    %PEPENUMERATE Does the backward enumeration for space defined by corner
    %point low, to high.
    
    currGen=java.util.Hashtable;
    currGen.put(key(high),high);
    
        while ~currGen.isEmpty
        
            %populate the next generation
            nextGen = java.util.Hashtable(Prefs.nodes*currGen.size);
            myNodes = currGen.elements;
            while myNodes.hasMoreElements
                state = myNodes.next;
                
                for i = 1:Prefs.nodes
                    if Prefs.type(i)
                        %move down
                        if state(i)>low(i)
                            state(i)=state(i)-1;
                            tempkey=key(state);
                            if DistributionValue(state,Prefs)>=Prefs.p
                                nextGen.put(tempkey,state);
                            end
                            state(i)=state(i)+1;
                        end %if
                    else
                        %move up
                        if state(i)<low(i)
                            state(i) = state(i)+1;
                            tempkey=key(state);
                            if DistributionValue(state,Prefs)>=Prefs.p
                                nextGen.put(tempkey,state);
                            end
                            state(i)=state(i)-1;
                        end %if
                    end %if Prefs.type(i)
                end %for nodes
            end
            
            %Compute conditional bounds (and record PEP)
            deleteList = java.util.Hashtable; 
            allNodes = nextGen.elements;
            while allNodes.hasMoreElements
                state=allNodes.next;
                
                if ConditionalBound(state,Prefs)
                   %record PEP and delete from nextGen
                    deleteList.put(key(state),state);
                    WritePEP(state);
                end
            end %while

            %update nextGen
            allDNodes = deleteList.elements;
            while allDNodes.hasMoreElements
                state=allDNodes.next;
                nextGen.remove(key(state));
            end
            
            fprintf(1,'Generation Size: %5.0f.\n',nextGen.size);
            currGen = nextGen;
            %force gc
            %java.lang.System.gc;
            
        end %while          
    end %function
