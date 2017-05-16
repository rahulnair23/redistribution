%RELOCATION_STATIC Solves the static fleet relocation problem
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
% OUTPUTS
% out : A structure containing
%   x : logical matrix of relocation variables
%   y : integer matrix of how many vehicles to relocate
%   fval : objective function value (relocation cost)
%   status : LPSOLVE status
%   alpha : (For consistency) = 1
%
% Rahul Nair
% 1/1/09
function [out] = relocation_static(LPdata)

nodes = length(LPdata.lambda1);
out.alpha = 1;
%Get Skellam mean values
xi_bar = LPdata.lambda1-LPdata.lambda2;

%% build the constraint matrix
%% Constraint 1: Vehicle demand constraint 
% v_j + \sum_i y_{ij} \ge Finv1
A1 = [];
for j = 1:nodes
    temp = zeros(nodes);
    temp(:,j) = 1;
    temp(j,j) = 0;
    A1 = vertcat(A1,reshape(temp,1,nodes*nodes));
end
A1 = horzcat(zeros(nodes,nodes*nodes),A1); %pad for x's
b1 = xi_bar(:) - LPdata.inventory(:);
e1 = ones(nodes,1);

%% Constraint 2: Spaces constraint 
A2 = [];
for j = 1:nodes
    temp = zeros(nodes);
    temp(j,:) = -1;
    temp(j,j) = 0;
    A2 = vertcat(A2,reshape(temp,1,nodes*nodes));
end
A2 = horzcat(zeros(nodes,nodes*nodes),A2); %pad for x's
b2 = xi_bar(:) + LPdata.capacity(:) - LPdata.inventory(:);
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

%% Constraint 4(ALT) : Cannot redistribute more than available
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

%% Binary and Integer constraints
bint = 1:nodes*nodes; %the x's
xint = (nodes*nodes)+1:2*nodes*nodes; %the y's

%% Variable bounds
y_max = sum(sum(LPdata.capacity)); 
vlb = zeros(2*nodes*nodes,1);
vub = [ones(nodes*nodes,1);y_max*ones(nodes*nodes,1)];
col_names = set_column_names(nodes);
%% Objective
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