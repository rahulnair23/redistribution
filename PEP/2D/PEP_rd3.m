function output=PEP_rd2(LB,UB,userPrefs)
%MAIN Procedure for generating PEP's using a recursive decomposition
%algorithm.
%   Searches a hyperrectangle defined by corner points LB and UB to get
%   points v and v', such that g(v,v')>=p and there exists no w and
%   w' such that g(w,w')>=p, where w<=v and w'>=v'. Uses a recursive decomposition to
%   disregard areas of the lattice that are guaranteed not to have a PEP.
%   
%INPUTS:
%LB     : Integer vector of the smallest node in the lattice
%UB     : Vector of the largest node in the lattice
%Prefs  : A structure containing the following fields (defaults)
%.p      : Desired probability level 
%.type   : Logical vector indicating the if r.v. is from an equation of the
%form P(Tx \ge \xi)\ge p or if it is P(Tx \le \xi)\ge p
%.jointDistribution : 
%A function handle that calculates the joint distribution
%value using parameters defined in data type Prefs.parameters (@SkellamPMF)
%.parameters : 
%Stores the parameters used in the joint distribution function
%.isPMF : Logical scalar, if the distribution function handle computes the
%cummulative function or just mass (default: 1)
%.maxtime : Maximum run time allowed in seconds (10800) 
%.minsize : Size at which the enumeration algorithm begins (min(10^5,2^n)4)
%.precompute: If the marginals are independent, the distribution value can
%be precomputed to save considerable time (default: 1)
%.marginalDistribution
%.algorithm: Determines the algorithm for local search. Possible values
%BackwardEnumeration, (Prekopa).
%.verbose: Boolean, status messages displayed (false)/true
%. WritePEPtoFile (true)/false record text file
%
% OUTPUTS: out a structure containing
%.S     : a matrix of PEP's, with each row representing a PEP.
%.time  : run time in seconds
%.reason : termination reason
%
% Rahul Nair, 4/3/09

%% Set up
global Prefs;           %clone of userPrefs
global nodes;           %number of stations
PEPCounter=1;           %the PEP counter
S = [];                 %Set of PEP's (each column is a PEP)
capacity = 10^5;        %Capacity of PEP set to start, incremented when exceeded
SHash = java.util.Hashtable(10^5);

if nargin==0
    %demo
    warning('PEPRD:Demo','Running demo problem only.');
    nodes = 2;
    domain = 30;
    LB=ones(nodes,1);
    UB = domain *ones(nodes,1);
    userPrefs.isPMF=false;
    userPrefs.verbose = false;
    userPrefs.algorithm = 'prekopa';
    userPrefs.type=true(1,nodes);
    userPrefs.jointDistribution = @JointDis;
    userPrefs.parameters = 13;
    userPrefs.precompute = false;
    userPrefs.WritePEPtoFile = false;
    userPrefs.p = 0.5;
    userPrefs.minsize = 10^50;
    %Prefs all at defaults
end

myRunTime=tic;
%store the visit
SearchMatrix = logical(sparse(domain,domain));


nodes=length(LB);
PEP_rd_prefs(userPrefs);
S=zeros(nodes,capacity); 
x = 1:30;
y = 1:30;
for ii=1:length(x);
    for jj = 1:length(y);
        z(ii,jj) = JointDis([ii jj],Prefs.parameters);
    end
end

surf(x,y,z);
hold all;
%relay('Count: %5.0f\n', prod(UB-LB+1));
%Create a look up table
if Prefs.precompute
    TightenedBounds;
    marginals=[];
    InitDistribution;
end

%Set up correct directions
startv=LB;
startv(~Prefs.type)=UB(~Prefs.type);
endv=LB;
endv(Prefs.type)=UB(Prefs.type);

%Call recursive PEP generator
PEPGenerator(startv,endv);

output.runtime = toc(myRunTime);
%Trim excess zeros from S
S=S(:,1:PEPCounter-1);
S=S';
output.S =S;

%plot the visit
%spy(SearchMatrix);

%% NESTED: Tighten Bounds
    function TightenedBounds
    %TIGHTENBOUNDS Resets the LB and UB vectors such that the minimum point
    %is guaranteed not to cut off PEP's.
        for i = 1:nodes
            t=UB;
            t(~Prefs.type)=LB(~Prefs.type);
            
            if Prefs.type(i)
                %work upwards
                for j = LB(i):UB(i)
                    t(i)=j;
                        if Prefs.jointDistribution(t,...
                                Prefs.parameters,Prefs.type)>=Prefs.p
                            LB(i)=j;
                            break
                        end %if
                end%for j
            else
                %work downwards
                for j=UB(i):-1:LB(i)
                    t(i)=j;
                        if Prefs.jointDistribution(t,...
                                Prefs.parameters,Prefs.type)>=Prefs.p
                            UB(i)=j;
                            break
                        end %if
                end %for j
            end
        end %for i
    end %Resets LB 

%% NESTED: Initialize distribution (look up matrix)
    function InitDistribution
        %INITDISTRIBUTION:Precomputes the marginals for each state
        %creates a matrix 'marginals' of size (max(UB-LB),nodes) with the
        %cummulative density for each dimension. 
        %marginals(k,i)=F_i(x<=LB(i)+k)
        smallest=-100;
        marginals = zeros(max(UB-LB),nodes);
        for i =1:nodes
            k=1;
            for j = LB(i):UB(i)
                if Prefs.isPMF
                    if k==1
                        z=smallest;
                        while z <=LB(i)
                            marginals(k,i)=marginals(k,i)+...
                                Prefs.marginalDistribution(z,Prefs.parameters{1}(i),...
                                                            Prefs.parameters{2}(i));
                                %TODO: Parameters are specific to current
                                %problem, need to change
                            z=z+1;
                        end
                    else
                        marginals(k,i)=marginals(k-1,i)+...
                                Prefs.marginalDistribution(j,Prefs.parameters{1}(i),...
                                                            Prefs.parameters{2}(i));
                    end %if
                    k=k+1;
                else
                    %its a CDF already
                    warning('PEPRD:InitDistribution:NotTested','Feature not tested yet.');
                    marginals(k,i) = Prefs.marginalDistribution(z,Prefs.parameters);
                    k=k+1;
                end %if isPMF
            end %for j
        end %i
     end % Creates look up table

%% NESTED: Computes the distrbution value
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
     if nargin==1
        plot3(s(1),s(2),val,'Marker','o','MarkerSize',9,'MarkerFaceColor','r','MarkerEdgeColor','r');
        drawnow;
        SearchMatrix(s(1),s(2))=true;
    end
    end %Computes F(.) by lookup

%% NESTED: Main recursive procedure
    function PEPGenerator(startv,endv)
    %PEPGENERATOR Generates a set of p-level efficient points within the
    %hyperrectangle defined by start,end
    %   Recursive
    %   startv: The lower corner of the hyperrec e.g. (0,0,...0)
    %   endv  : The upper corner of e.g. (1,1,..1)
    
   %The hyperrectangle is "above" the p-frontier
    if DistributionValue(startv)>Prefs.p
        return
    end

    %The hyperrectangle is "below the p-frontier
    if DistributionValue(endv)<=Prefs.p
        %ignore
        return
    end

    %The p-frontier passes this zone
    %recurse
    
        w = findsplit(startv,endv);
        if ~isempty(w)
            %fprintf(1,'Split\n');
            for i = 2^nodes-1:-1:0
                %generate the n'th quadrant and recurse
                [l,u]=NQuadrant(startv,w,endv,i);
                PEPGenerator(l,u);

                %check for time
                if toc(myRunTime)>Prefs.maxtime
                    %we've run out of time
                    break
                end
            end
        else
            %fprintf(1,'Enum\n');
            %do a local search
            switch lower(Prefs.algorithm)
                case 'completeenumeration'
                    CompleteEnumeration(startv,endv);
                case 'backwardenumeration'
                    PEPEnumerate(startv,endv);
                case 'prekopa'
                %    error('PEPRD:NotTested','Not tested');
                    Prekopa(startv,endv);
                otherwise
                    error('PEPRD:Prefs:Algorithm','Check Prefs.algorithm.');
            end %switch
        end %if ~isempty(w)
    end %Recursive

%% NESTED: Complete Enumeration
    function CompleteEnumeration(s,e)
        
        %Generate all nodes
        allnodes = allVertices(s,e);
        for b =1:size(allnodes,1)
            curr_state = allnodes(b,:);
            %check for PEP
            if ConditionalBound(curr_state')
                WritePEP(curr_state,DistributionValue(curr_state'));
            end
        end %for b
    end %function

%% NESTED: Splitting function
    function w=findsplit(s,e)
        
        if prod(s-e+1)>Prefs.minsize
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

%% NESTED: Backward Enumeration algorithm
    function PEPEnumerate(low,high)
    %PEPENUMERATE Does the backward enumeration for space defined by corner
    %point low, to high.
    currGen=java.util.Hashtable;
    currGen.put(key(high),high);
    
        while ~currGen.isEmpty
        
            %populate the next generation
            nextGen = java.util.Hashtable(nodes*currGen.size);
            myNodes = currGen.elements;
            while myNodes.hasMoreElements
                state = myNodes.next;
                
                for i = 1:nodes
                    if state(i)>low(i)
                        state(i)=state(i)-1;
                        tempkey=key(state);
                        if DistributionValue(state)>=Prefs.p
                            nextGen.put(tempkey,state);
                        end
                        state(i)=state(i)+1;
                    end %if
                end %for
            end
            
            %Compute conditional bounds (and record PEP)
            deleteList = java.util.Hashtable; 
            allNodes = nextGen.elements;
            while allNodes.hasMoreElements
                state=allNodes.next;
    
               if ConditionalBound(state)
                   %record PEP and delete from nextGen
                    deleteList.put(key(state),state);
                    WritePEP(state,DistributionValue(state));
                end
            end %while

            %update nextGen
            allDNodes = deleteList.elements;
            while allDNodes.hasMoreElements
                state=allDNodes.next;
                nextGen.remove(key(state));
            end
            %fprintf(1,'Generation Size: %6.0f.\n',nextGen.size);
            
            currGen = nextGen;
            %force gc
            %java.lang.System.gc;
            
        end %while          
    end %function

%% NESTED: Prekopa's Enumeration Algorithm    
    function Prekopa(low,high)
    %PREKOPA Implements Prek\'{o}pa's PeP enumeration algorithm for the local
    %region defined by low,high
    %relay('Prekopa algorithm: Search space %4.0f\n',prod(high-low+1));
        %Call recursive PEP generator
        localPEP=[]; %set of local PEP's
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
               k_itr=k_itr+1;
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

%% NESTED: Conditional Bound Calculation
    function result = ConditionalBound(state)
    %CONDITIONALBOUND Computes the conditional bound, l(v). Returns true if
    %the bound is met l(v)=v; and false otherwise.
    %   See Beraldi and Rusz..2001

    
        if DistributionValue(state,true)>=Prefs.p
            %assume it to be PEP
            result = true;
        else
            %It is not PEP, so reject
            result = false;
            return
        end
    
    %The actual conditional bound check
    %NOTE: This doesn't actually compute the entire conditional bound l(v),
    %for that see earlier versions. This just checks one dimension lower to
    %ensure that there is no smaller vector that satisfies the PEP
    %condition, in which case the current point is not PEP.
        
        for c = 1:nodes
            
            %store current state
            temp = state(c);

            %increment state
            if Prefs.type(c)
                %check for a smaller value
                if state(c)>LB(c)
                    state(c) = state(c)-1;
                else
                    continue
                end
                
            else
                
                %check for a larger value
                if state(c)<UB(c)
                    state(c) = state(c)+1;
                else
                    %skip this dimension
                    continue
                end
            end

            %do the check
            if DistributionValue(state,true)>=Prefs.p
                result = false;
                break
            end

            %correct state for next iteration
            state(c) = temp;
            
        end %for   
    end % function

%% NESTED: Quadrant generation 
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

%% NESTED Write PEP's
    function WritePEP(state,val)
    %WRITEPEP Writes the PEP's to file tab separated and sets S
    %if ~ismember_myversion(state(:),S)
    if ~SHash.containsKey(key(state))
        if size(S,2)<PEPCounter
            %increase capacity
            S(:,end+capacity)=0;
        end
        S(:,PEPCounter)=state(:);
        PEPCounter=PEPCounter+1;
        %fprintf(1,'%5.0f\n',PEPCounter);
        plot3(state(1),state(2),...
            Prefs.jointDistribution(state,Prefs.parameters),...
            'Marker','o','MarkerSize',12,'MarkerFaceColor','b','MarkerEdgeColor','b');
        drawnow;
        %add to hashtable
        SHash.put(key(state),1);
        if Prefs.WritePEPtoFile
                pepFile = fopen('PEP.txt','at');
                    fprintf(pepFile,'%3.0f\t',state);
                    if nargin==2
                        fprintf(pepFile,'(%0.3f)\n',val);
                    else
                        fprintf(pepFile,'\n');
                    end
                fclose(pepFile);
        end
    end
end %futnction
end


%% Handle Preferences
function []= PEP_rd_prefs(userPrefs)
global nodes Prefs;
%set up defaults
Prefs.p = 0.9;
Prefs.type = ones(nodes,1);
Prefs.jointDistribution = @MyJointDistributionFunction;
Prefs.marginalDistribution = @SkellamPMF;
Prefs.parameters = {ones(nodes,1),2*ones(nodes,1)};
Prefs.isPMF=true;
Prefs.maxtime = 10800;
Prefs.minsize = 10^4;
Prefs.precompute = true;
Prefs.algorithm = 'backwardenumeration';
Prefs.verbose = false;
Prefs.WritePEPtoFile = true;

%override defaults
 Fields = fieldnames(Prefs);
    for i = 1:length(Fields)
        if isfield(userPrefs, Fields{i})
            Prefs.(Fields{i}) = userPrefs.(Fields{i});
        end % isfield
    end % i
end %function
%% Hash functions
function out = key(item)
%get row vectors
if size(item,1)~=1
    item = item';
end
%Alt 1:
out = char(100+item);
%Alt 2:
% out=0;
% t=length(item);
% for i = 1:t
%     out = out+10^(t-i)*item(i);
% end
end

%% Vertex
function v = aVertex(start,endv,n)
%AVERTEX Generate i'th vertex of hypercube defined by s,e
    ve = ('1'==dec2bin(n,length(start)));
    v=start;
    v(ve==1)=endv(ve==1);
end

%% Reporting
function [] = relay(TextItem,varargin)
% RELAY (TextItem,VariableItem) : This function outputs to both screen and logfile
% 	Essentially, this command duplicates the work of two fprintf commands, and
% 	checks to see that the disk file is being written correctly.
% INPUTS:
%	TextItem		String	Text to print, including formats of variables as per fprintf command
%	varargin		Any		(Optional) The associated variable to print, as per fprintf command

global Prefs

    if Prefs.verbose
        LogFile1 = fopen('Runlog.txt','at');

            Num_Bytes_Written = 0;
            if nargin > 1	% For printing output that includes a variable
                try
                    fprintf(1,TextItem,varargin{:});
                    Num_Bytes_Written = fprintf (LogFile1,TextItem, varargin{:});
                end
            else
                try
                    Num_Bytes_Written = fprintf (LogFile1,TextItem);
                    fprintf(1,TextItem);
                end
            end
        fclose(LogFile1);
    end
end %func

%% Clone of dec2bin
function out= mydec2bin(d,n)
out=uint8(rem(floor(d(:)*pow2(1-n:0)),2));
end
%% Clone of ismember
function out=ismember_myversion(state,S)
A=bsxfun(@eq,state,S);
out=any(sum(A)==length(state));
end
%% All Vertices
function out = allVertices(l,u)
%NDVERTICES Generates all the vertices of the hyperrectangle defined by the
%corner vectors l, and u.

n=length(l);

    v = mydec2bin(0:2^n-1,n);
    out =zeros(size(v));
    for i = 1:n
        out(v(:,i)==0,i)=l(i);
        out(v(:,i)==1,i)=u(i);
    end
 %remove duplicates
 if any(l==u)
    out = unique(out,'rows');
 end
end %function