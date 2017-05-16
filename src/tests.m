

%% Load cplex libraries
NET.addAssembly('C:\ILOG\CPLEX101\bin\x86_win32\ILOG.CPLEX.dll');
NET.addAssembly('C:\ILOG\CPLEX101\bin\x86_win32\ILOG.Concert.dll');

%Define variables
myCplex = ILOG.CPLEX.Cplex;
x = myCplex.NumVar(0,10);
y = myCplex.NumVar(0,10);

%set up a obj variable
o = myCplex.LinearNumExpr;
o.AddTerm(5,x);
o.AddTerm(3,y);
myCplex.AddMaximize(o);

%set up constraint
c1 = myCplex.LinearNumExpr;
c1.AddTerm(1,x);
c1.AddTerm(3,y);
myCplex.AddLe(c1,20);

if myCplex.Solve
    %get the values
    fprintf(1,'Output \n x : %3.2f\n y : %3.2f\n',...
        myCplex.GetValue(x),myCplex.GetValue(y));
end
    
myCplex.End;

%% arrays
myCplex = ILOG.CPLEX.Cplex;

%create a string array in NET
names = NET.createArray('System.String',5);
for i = 0:4
    names.Set(i,['x' num2str(i)]);
end
x = myCplex.BoolVarArray(5,names);
obj_coeffs = [5 3 5 3 5];
obj_c = NET.convertArray(obj_coeffs,'System.Double',5);

%set up objective
obj = myCplex.LinearNumExpr;
for i = 0:4
    obj.AddTerm(obj_coeffs(i+1),x.Get(i));
end
myCplex.AddMaximize(obj);

%set up constraint

if myCplex.Solve
    %get the values
    fprintf(1,'Output:\n');
        for i =0:4
            fprintf(1,'x%1.0f : %2.2f\n',i,myCplex.GetValue(x.Get(2)));
        end
end
    
myCplex.End;
clear;






