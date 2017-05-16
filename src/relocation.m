%RELOCATION Solves the fleet relocation problem
%   Uses the LPSOLVE package
% INPUTS
%  data: A structure containing
%   
%   lambda1 : vector of Poisson rates for user arrivals at station
%   lambda2 : vector of Poisson rates for returns at station
%   cost : matrix of relocation costs between two stations
%   capacity : vector of station capacities (docks)
%   inventory : vector of current inventories at station
%   alpha : The reliability factor 0<alpha<1
%
% userPrefs: A structure containing
%   static: Boolean, indicating whether solve using the mean,
%           can be used to compute the Value of stochastic solution  
%           (false)/true
% 
% OUTPUTS
% out : A structure containing
%   x : logical matrix of relocation variables
%   y : integer matrix of how many vehicles to relocate
%   fval : objective function value (relocation cost)
%   status : LPSOLVE status
%   alpha : True reliability achieved
%
% RN

function [out] = relocation(LPdata)
%% Error checking
nodes = LPdata.nodes;

%% Determine the reliability allocation
% Apportioned equally
nn = nodes-1; %exclude the super sink
alpha_j = ((nn-1+LPdata.alpha)/nn)*ones(nodes,1);
%note: alpha_j(1) is not needed. provided for index consistency later


%% Generate inverses
% -Skellam distributed unless one process non-existent, in which case it is
% Poisson.
% -Capacity infeasibility handled, but reducing the reliability by
% step_size, till the infeasibility is resolved.

step_size = 0.0005;
Finv1 = zeros(nodes,1);
Finv2 = zeros(nodes,1);
for i = 2:nodes %the first is always the depot
    
    %edit: 2/25/2009
    %Verify that both Checkotuts and Returns are performed
     if (LPdata.lambda1(i)==0||LPdata.lambda2(i)==0)
         %The aggregate flow reduces to a poisson process
         if LPdata.lambda1(i)==0
            %No Checkouts
            relay('Port %2.0f has no checkout process.\n',i);
            Finv1(i) = 0;
            Finv2(i) = - poissinv((1+alpha_j(i))/2,LPdata.lambda2(i));
            %Check capacity constraint
            while -Finv2(i)> LPdata.capacity(i)
                alpha_j(i) = alpha_j(i)-step_size;
                Finv2(i) = - poissinv((1+alpha_j(i))/2,LPdata.lambda2(i));               
            end %while
            relay('Port %2.0f: alpha_j reduced to %1.2f.\n',i,alpha_j(i));
         else
             %No Returns
             relay('Port %2.0f has no return process.\n',i);
             Finv1(i) = poissinv((1+alpha_j(i))/2,LPdata.lambda1(i));
             Finv2(i) = 0;
             
             %Check capacity constraint
             while Finv1(i)>LPdata.capacity(i)
                 alpha_j(i) = alpha_j(i)-step_size;
                 Finv1(i) = poissinv((1+alpha_j(i))/2,LPdata.lambda1(i));
             end %while
              relay('Port %2.0f: alpha_j reduced to %1.2f.\n',i,alpha_j(i));
         end
     else
         %generate Skellam variables
        [prob,range] = skellamcdf(LPdata.lambda1(i),LPdata.lambda2(i),-100:100);
        if any(isnan(prob))
            relay('Can''t compute the Skellam CDF for port %2.0f. Terminating.\n',i);
            return;
        end
        Finv1(i) = skellaminvCDF(range,prob,((1+alpha_j(i))/2));
        Finv2(i) = skellaminvCDF(range,prob,((1-alpha_j(i))/2));

        %check capacity infeasibilities
%              while max(Finv1(i),0)-min(Finv2(i),0)>LPdata.capacity(i)
%                  alpha_j(i) = alpha_j(i) - step_size;
%                  Finv1(i) = skellaminvCDF(range,prob,((1+alpha_j(i))/2));
%                  Finv2(i) = skellaminvCDF(range,prob,((1-alpha_j(i))/2));   
% %             relay('Port: %2.0f alpha_j : %0.3f (%2.0f,%2.0f)\n',...
% %                     i,alpha_j(i),...
% %                     max(Finv1(i),0)-min(Finv2(i),0),LPdata.capacity(i));
%              end %while
        %Check capacity infeasibilities (ALT)
            while max(Finv1(i),0)-min(Finv2(i),0)>LPdata.capacity(i)
                %failure rates
                t1 = (1+alpha_j(i))/2;
                t2 = (1-alpha_j(i))/2;
                switch true
                    case max(Finv1(i),-Finv2(i))==Finv1(i)
                        t1= t1-step_size;
                        Finv1(i) = skellaminvCDF(range,prob,t1);                        
                    case max(Finv1(i),-Finv2(i))==-Finv2(i)
                        t2=t2+step_size;
                        Finv2(i) = skellaminvCDF(range,prob,t2);                        
                end %switch
                alpha_j(i) = t1-t2;
%                 relay('Port: %2.0f t1 : %0.3f t2 : %0.3f alpha_j : %0.3f (%2.0f,%2.0f)\n',...
%                     i,t1,t2,alpha_j(i),...
%                     max(Finv1(i),0)-min(Finv2(i),0),LPdata.capacity(i));
            end %while
        
            relay('Port %2.0f: alpha_j reduced to %1.3f.\n',i,alpha_j(i));
     end %if 
end

%% Update reliability
% may be updated later if there are supply infeasibilities
out.alpha = 1 - sum(ones(nodes-1,1)-alpha_j(2:nodes));

relay('System reliability (after Capacity infeasibilities): %0.2f.\n',out.alpha);

%% Check for supply infeasibilities
% \sum max(Finv1,0) < \sum v_i
% \sum -min(Finv2,0) < \sum (c_i-v_i)

if (sum(max(Finv1(2:end),0))>sum(LPdata.inventory(2:end)) || ...
        -sum(min(Finv2(2:end),0))>sum(LPdata.capacity(2:end)-LPdata.inventory(2:end)))
    %There will be supply infeasibilities
    relay('There are supply infeasibilities. Recovering partial redistribution.\n');
else
    relay('No supply infeasibilities this period.\n');
end %if 

%% build the constraint matrix
%% Constraint 1: Vehicle demand constraint 
% v_j + \sum_i y_{ij} - \sum_i y_{ji} \ge Finv1
A1 = [];
for j = 1:nodes
    temp = zeros(nodes);
    temp(:,j) = 1;
    temp(j,:) = -1; %edit 6 Feb 09 (account for relcations out)
    temp(j,j) = 0;
    A1 = vertcat(A1,reshape(temp,1,nodes*nodes));
end
A1 = horzcat(zeros(nodes,nodes*nodes),A1); %pad for x's
b1 = Finv1(:) - LPdata.inventory(:);
e1 = ones(nodes,1);

%% Constraint 2: Spaces constraint 
% -\left( c_j - v_j + \sum_i y_{ji} - \sum_i y_{ij} \right) \le Finv2
A2 = [];
for j = 1:nodes
    temp = zeros(nodes);
    temp(j,:) = -1;
    temp(:,j) = 1; %edit 6 Feb 09 (account for relocations into station)
    temp(j,j) = 0;
    A2 = vertcat(A2,reshape(temp,1,nodes*nodes));
end
A2 = horzcat(zeros(nodes,nodes*nodes),A2); %pad for x's
b2 = Finv2(:) + LPdata.capacity(:) - LPdata.inventory(:);
e2 = -1*ones(nodes,1);

%% Constraint 3: Relate dv's
A3 = [];
for i = 1:nodes
    for j = 1:nodes
        temp_x = zeros(nodes);
        temp_y = zeros(nodes);
        temp_x(i,j) = -LPdata.capacity(i);
        temp_y(i,j) = 1;
        A3 = vertcat(A3,horzcat(reshape(temp_x,1,nodes*nodes),...
                             reshape(temp_y,1,nodes*nodes)));
    end
end
b3 = zeros(nodes*nodes,1);
e3 = -1*ones(nodes*nodes,1);

%% Constraint 4 : Cannot redistribute more than available
A4 = [];
for j = 1:nodes
    temp = zeros(nodes);
    temp(j,:) = 1;
    temp(j,j) = 0;
    A4 = vertcat(A4,reshape(temp,1,nodes*nodes));
end
A4 = horzcat(zeros(nodes,nodes*nodes),A4); %pad for x's
e4 = -1*ones(nodes,1);
b4 = LPdata.inventory(:);

%% Constraint 5: Cannot move if no spaces exist
A5 = [];
for j = 1:nodes
    temp_y = zeros(nodes);
    %temp_y(j,:) = -1; %outgoing
    temp_y(:,j) = 1; %incoming
    temp_y(j,j) = 0; %self loop
    A5 = vertcat(A5,reshape(temp_y,1,nodes*nodes));
end 
A5 = horzcat(zeros(nodes,nodes*nodes),A5);
e5 = -1*ones(nodes,1);
b5 = LPdata.capacity(:) - LPdata.inventory(:);

%% Other inputs
% Binary and Integer constraints
bint = 1:nodes*nodes; %the x's
xint = (nodes*nodes)+1:2*nodes*nodes; %the y's

% Variable bounds
y_max = sum(sum(LPdata.capacity)); 
vlb = zeros(2*nodes*nodes,1);
vub = [ones(nodes*nodes,1);y_max*ones(nodes*nodes,1)];
col_names = set_column_names(nodes);

% Objective
Obj = [reshape(LPdata.cost,1,nodes*nodes), 0.1 *ones(1,nodes*nodes)];

%% feed to LPSOLVE
A = vertcat(A1,A2,A3,A4,A5);
e = vertcat(e1,e2,e3,e4,e5);
b = vertcat(b1,b2,b3,b4,b5);
[fval,xy,duals,stat,var_names] = lp_solve(Obj,A,b,e,vlb,vub,xint,bint,'',[],[],col_names);
out.status = lp_solve_status(stat);
    if ~isempty(xy)
        temp = interpret(xy,var_names);
        out.x = temp.x;
        out.y = temp.y;
        out.fval = fval;
    else
        out.fval = [];
    end

%% interpret 
    function t = interpret(xy,var_names)
        x = zeros(nodes);
        y = zeros(nodes);
        
        for ii = 1:length(var_names)
            %process name
            [var,r,c] = process_name(var_names{ii});
            if var=='x'
                x(r,c) = xy(ii);
            else
                y(r,c) = xy(ii);
            end
        end
        t.x = logical(x);
        t.y = y;
        
        %TODO: see if the super node was used.
    end %function
end %function

%% Set column names for debug
function colnames = set_column_names(n)
x = cell(n,n);
y = cell(n,n);

for i = 1:n
    for j = 1:n
        x(i,j) = {['x(' num2str(i) ',' num2str(j) ,')']};
        y(i,j) = {['y(' num2str(i) ',' num2str(j) ,')']};
    end
end
colnames = horzcat(reshape(x,1,n*n),reshape(y,1,n*n));
end %function

%% retrieve variable info
function [var,r,c] = process_name(name)
if name(1) == 'x'
    var = 'x';
else
    var = 'y';
end
A = sscanf(name(2:end),'(%d,%d)');
r = A(1);
c = A(2);
end %function

%% LP_SOLVE status
function msg = lp_solve_status(status)
switch status
    case -2
        msg =  'Out of memory';
    case 0
        msg = 'Optimal solution found';
    case 1
        msg = 'Suboptimal: Optimality not guaranteed';
    case 2
        msg = 'The problem is infeasible';
    case 3
        msg = 'The model is unbounded';
    case 4
        msg = 'The model is degenerative';
    case 5
        msg = 'Numerical failure encountered';
    case 6
        msg = 'User abort. See put_abortfunc';
    case 7
        msg = 'A timeout occured';
    case 9
        msg = 'The model could be solved by presolve. This can only happen if presolve is active via set_presolve';
    case 10
        msg = 'The Branch & Bound  routine failed';
    case 11
        msg = 'The B&B was stopped because of a break-at-first (see set_break_at_first) or a break-at-value (see set_break_at_value)';
    case 12
        msg = 'A feasible B&B solution was found';
    case 13
        msg = 'No feasible B&B solution found';
    otherwise
        msg = 'Termination code unknown.';
end %switch
end