 function WritePEP(state,val)
    %WRITEPEP Writes the PEP's to file tab separated and sets S
    %if ~ismember_myversion(state(:),S)
    
    
    %if ~Prefs.SHash.containsKey(key(state))
%         if size(Prefs.S,2)<Prefs.PEPCounter
%             %increase capacity
%             S(:,end+capacity)=0;
%         end
%         S(:,PEPCounter)=state(:);
%         PEPCounter=PEPCounter+1;
%         
%         %add to hashtable
%         SHash.put(key(state),1);

        pepFile = fopen('PEP.txt','at');
            fprintf(pepFile,'%3.0f\t',state);
            if nargin==2
                fprintf(pepFile,'(%0.3f)\n',val);
            else
                fprintf(pepFile,'\n');
            end
        fclose(pepFile);
    %end
end %futnction