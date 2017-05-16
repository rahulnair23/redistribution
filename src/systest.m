function out = systest (userData,userPrefs)
%SYSTEST Simulates the redistribution for multiple periods
%
%INPUTS (all optional)
%Needs following inputs default in (brackets)
%   userData : Structure containing the following fields
%   .periods    : Number of time periods to simulate (10)
%   .nodes      : Number of stations in the problem (3)
%   .lambda1    : Matrix of Poisson rates for vehicle checkouts
%                 (periods,nodes) unless cyclic is true in which case it
%                 can be of the size (states,nodes)
%   .lambda2    : Matrix of Poisson rates for vehicle returns (same as
%                 lambda1 above)
%   .cyclic     : logical to indicate if network cycles between states, the
%                 simulation cycles through each state till 'periods' is
%                 reached
%   .states     : Number of states
%   .alpha      : Desired system reliability (0.9)
%   .inventory  : Vector of inventory at stations, (1,nodes)
%   .capacity   : Vector of capacities at stations (1,nodes)
%   .cost       : Matrix of relocation costs (nodes,nodes)
%
%   userPrefs : Structure containing fields
%   .verbose    : logical if text output printed (true)/false
%   .static     : Logical, true if solving with expected value (false)
%   .seed       : Random number generator seed

%OUTPUT
%Returns a structure out that has
%   .data       :input data for debugging
%   .periods    :number of time periods for simulation
%   .cost       :vector with reloc costs for each period
%   .alpha      :Reliability achieved for time period (vector)
%   .balkVeh    :Number of vehicle requests balked
%   .balkSpace  :Number of space requests balked
%
% RN Jan09

%% globals 
global Prefs
global data

%% Set in/out put
%inputs
if nargin == 0
    userData = [];
    userPrefs = [];
elseif nargin == 1
    userPrefs = [];
end
set_inputs(userData,userPrefs);
periods = data.periods;
states = data.states;
nodes = data.nodes;
rand('seed',Prefs.seed); % runs are repeatable

LPdata = data; %clone for LPsolve;

%outputs
out.data = data;
out.periods = periods;
out.obj_val = zeros(1,periods);

if Prefs.static
     relay('\n\n\nWARNING: Running static model.\n\n\n');
end
%% Main

for p = 1:periods
    
    if data.cyclic
        %pull up correct state
        s = mod(p,states);
        if s==0
            s = states;
        end
    else
        s = p;
    end
    
    %print basic stats
    if Prefs.verbose
        relay('--------------------------\n');
        relay('Time Period %2.0f (State: %1.0f)\n',p,s);
        relay('Available Inventory: ');
        relay(' %3.0f ',LPdata.inventory);
        relay(1,'\n');
    end
   
    % set up data
    LPdata.lambda1 = data.lambda1(s,:);
    LPdata.lambda2 = data.lambda2(s,:);
    
    %set the problem version
    if Prefs.static
        up.static = true;
    else
        up.static = false;
    end

    %solve the LP
    LPout = relocation(LPdata);
            
    if isempty(LPout.fval)
       errflag = true;
    else
        errflag = false;
    end
    
    %intepret results
    if ~errflag
        curr_inventory = LPdata.inventory;
        for i =1:nodes
            for j = 1:nodes
                if LPout.y(i,j)>0
                    relay('Relocate %d units from %d to %d.\n',LPout.y(i,j),i,j);
                    curr_inventory(i) = curr_inventory(i) - LPout.y(i,j);
                    curr_inventory(j) = curr_inventory(j) + LPout.y(i,j);
                    %compute relocation costs
                    out.obj_val(p) = out.obj_val(p) + data.cost(i,j);
                end
            end
        end
        out.alpha(p) = LPout.alpha; %record true reliability achieved
    else
        %should not happen, but just in case, notify
        relay('Error solving this period. Assuming no redistribution.LPSOLVE msg: \n');
        relay(LPout.status);
        relay('\n');
    end
    
    relay('Relocated Inventory: ');
    relay(' %3.0f ',curr_inventory);
    relay(1,'\n');
    
    %generate random demands
    rnd_demand = zeros(1,nodes);
    for i = 2:nodes %no demand at the super node
        %generate cdf
        [y,x] = skellamcdf(LPdata.lambda1(i),LPdata.lambda2(i));
        %get a random variate
        rnd_demand(i) = skellaminvCDF(x,y,rand);
    end
    
    relay('      Random demand: ');
    relay(' %3.0f ',rnd_demand);
    relay(1,'\n');
    %adjust current inventory to reflect random event
    curr_inventory = curr_inventory - rnd_demand;
    %record disappointments
    out.balkVeh(p) = -1*sum(curr_inventory(curr_inventory<0));
    out.balkSpace(p) = sum(curr_inventory(curr_inventory>data.capacity)-...
                    data.capacity(curr_inventory>data.capacity));
    
    %correct imbalances
    curr_inventory(curr_inventory<0)=0;
    curr_inventory(curr_inventory>data.capacity)=...
        data.capacity(curr_inventory>data.capacity); %vehicles
    
    %reset for next period
    LPdata.inventory = curr_inventory;
end
        
end %function

%% Handle Preferences
function set_inputs(userData,userPrefs)
global Prefs data

%Default data
% data.nodes = 3;
% data.states = 2;
% data.cyclic = true;
% data.periods = 10;
% data.lambda1 = [0 5 2; 0 2 5]; 
% data.lambda2 = [0 2 5; 0 5 2];
% data.inventory = [20 3 3]; 
% data.capacity = [100 8 8];
% data.cost = [0.0 100 100;
%              100 0.0 1.0;
%              200 1.0 0.0];
% data.alpha = 0.9;

data.nodes = 8;
data.states = 2;
data.cyclic = true;
data.periods = 10;
data.lambda1 = [0 5 2 5 2 5 2 5;
                0 2 5 2 5 2 5 2]; 
data.lambda2 = [0 2 5 2 5 2 5 2; 
                0 5 2 2 5 2 5 2];
data.inventory = [20 3 3 3 3 3 3 3]; 
data.capacity = [100 8 8 8 8 8 8 8];
data.cost =  [0  100.0000  100.0000  100.0000  100.0000  100.0000  100.0000  100.0000;
  100.0000         0    0.9157    0.7577    0.0462    0.3816    0.7547    0.3404;
  100.0000    0.1576         0    0.7431    0.0971    0.7655    0.2760    0.5853;
  100.0000    0.9706    0.9595         0    0.8235    0.7952    0.6797    0.2238;
  100.0000    0.9572    0.6557    0.6555         0    0.1869    0.6551    0.7513;
  100.0000    0.4854    0.0357    0.1712    0.3171         0    0.1626    0.2551;
  100.0000    0.8003    0.8491    0.7060    0.9502    0.4456         0    0.5060;
  100.0000    0.1419    0.9340    0.0318    0.0344    0.6463    0.4984         0];
data.alpha = 0.9;


Prefs.verbose = true;
Prefs.static = false;
Prefs.seed = 0;

% Override data fields with user values
Fields = fieldnames(data);
for i = 1:length(Fields)
	if isfield(userData, Fields{i})
		data.(Fields{i}) = userData.(Fields{i});
	end % isfield
end 
% Override prefs fields with user values
Fields = fieldnames(Prefs);
for i = 1:length(Fields)
	if isfield(userPrefs, Fields{i})
		Prefs.(Fields{i}) = userPrefs.(Fields{i});
	end % isfield
end 

end %function

%% relay
function [] = relay(TextItem,varargin)
global Prefs
    if Prefs.verbose
        if nargin > 1	
            % For printing output that includes a variable
            fprintf (TextItem, varargin{:});
        else
            fprintf (TextItem);
        end
    end % Prefs.verbose
end %relay