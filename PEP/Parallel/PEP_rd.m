function output=PEP_rd(uLB,uUB,userPrefs)
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

% OUTPUTS: out a structure containing
%.S     : a matrix of PEP's, with each row representing a PEP.
%.time  : run time in seconds
%.reason : termination reason
%
% Rahul Nair, 4/3/09

%% Set up
global Prefs;           %clone of userPrefs
%global PEPCounter;           %the PEP counter
%global S;                 %Set of PEP's (each column is a PEP)
%global SHash;           %Hashtable to store the PEP's
%global marginals;

%start parallel solver
matlabpool('open',4,'FileDependencies',{'D:\rahul\SpeedTest\ver4'});
%PEPCounter = 1;
%S = [];
%SHash = java.util.Hashtable(10^5);

if nargin==0
    %demo
    warning('PEPRD:Demo','Running demo problem only.');
    userPrefs.nodes = 6;
    uLB=-10*ones(userPrefs.nodes,1);
    uUB = 10 *ones(userPrefs.nodes,1);
    userPrefs.isPMF=true;
    userPrefs.verbose = false;
    userPrefs.algorithm = 'backwardenumeration';
    userPrefs.minsize=(2^userPrefs.nodes)+1;
    userPrefs.type=logical([ones(1,userPrefs.nodes/2) zeros(1,userPrefs.nodes/2)]);
    userPrefs.p = 0.9;
    userPrefs.minsize = 1000;
    userPrefs.LB = uLB;
    userPrefs.UB = uUB;

    %Prefs all at defaults
end

LB=uLB;
UB=uUB;

myRunTime=tic;

PEP_rd_prefs(userPrefs);

TightenedBounds;


%relay('Count: %5.0f\n', prod(UB-LB+1));
%Create a look up table
if Prefs.precompute
    InitDistribution;
end

%Set up correct directions
%startv=LB;
%startv(~Prefs.type)=UB(~Prefs.type);
%endv=LB;
%endv(Prefs.type)=UB(Prefs.type);
startv=LB;
endv = UB;
%Call recursive PEP generator
PEPGenerator(startv,endv,Prefs);

matlabpool close

output.runtime = toc(myRunTime);

%TODO:Retrive all PEP's


%% NESTED: Tighten Bounds
    function TightenedBounds
    %TIGHTENBOUNDS Resets the LB and UB vectors such that the minimum point
    %is guaranteed not to cut off PEP's.
        for i = 1:Prefs.nodes
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
        Prefs.LB = LB;
        Prefs.UB = UB;
    end %Resets LB 

%% NESTED: Initialize distribution (look up matrix)
    function InitDistribution
        %INITDISTRIBUTION:Precomputes the marginals for each state
        %creates a matrix 'marginals' of size (max(UB-LB),nodes) with the
        %cummulative density for each dimension. 
        %marginals(k,i)=F_i(x<=LB(i)+k)
        smallest=-100;
        marginals = zeros(max(UB-LB),Prefs.nodes);
        for i =1:Prefs.nodes
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
        Prefs.marginals = marginals;
     end % Creates look up table

%% NESTED: Computes the distrbution value
   %moved outside

%% NESTED: Main recursive procedure
   
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
    %moved outside

%% NESTED: Get Count
    %Moved outside with splitting function
    
%% NESTED: Backward Enumeration algorithm
   %Moved outside
   
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

%% NESTED: Conditional Bound Calculation
  
%% NESTED: Quadrant generation 
   %Moved outside for Parallel Toolbox

%% NESTED Write PEP's
   %moved outside
end


%% Handle Preferences
function []= PEP_rd_prefs(userPrefs)
global Prefs;
capacity=10^5;

%set up defaults
Prefs.p = 0.9;
Prefs.nodes = 10;
Prefs.type = ones(Prefs.nodes,1);
Prefs.jointDistribution = @MyJointDistributionFunction;
Prefs.marginalDistribution = @SkellamPMF;
Prefs.parameters = {ones(Prefs.nodes,1),2*ones(Prefs.nodes,1)};
Prefs.isPMF=true;
Prefs.maxtime = 10800;
Prefs.minsize = 500;
Prefs.precompute = true;
Prefs.algorithm = 'backwardenumeration';
Prefs.verbose = false;
Prefs.LB = ones(Prefs.nodes,1);
Prefs.UB = 10*ones(Prefs.nodes,1);
Prefs.PEPCounter = 1;
Prefs.S =zeros(Prefs.nodes,capacity); 
Prefs.SHash = java.util.Hashtable(capacity);

%override defaults
 Fields = fieldnames(Prefs);
    for i = 1:length(Fields)
        if isfield(userPrefs, Fields{i})
            Prefs.(Fields{i}) = userPrefs.(Fields{i});
        end % isfield
    end % i
end %function

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