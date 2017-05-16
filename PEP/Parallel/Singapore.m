function out=Singapore
%SINGAPORE Generates PEP for singapore
load data.mat;

nodes = 28;
lambda1 = [data.lambda1(3,2:end) data.lambda1(3,2:end)];
lambda2 = [data.lambda2(3,2:end) data.lambda2(3,2:end)];
startLB=[5*ones(nodes/2,1); -6*ones(nodes/2,1)];
startUB=[7*ones(nodes/2,1); -5*ones(nodes/2,1)];
type = [true(1,nodes/2) false(1,nodes/2)];

%Set up preferences
userPrefs.nodes = nodes;
userPrefs.parameters = {lambda1,lambda2};
userPrefs.verbose = true;
userPrefs.maxtime =  604800;
userPrefs.algorithm = 'backwardenumeration';
userPrefs.type = type;
userPrefs.minsize = 10^12;%6*10^7;
userPrefs.p = 0.9;

%solve
%PEPEnumerate(startLB,startUB,userPrefs);
out=PEP_rd(startLB,startUB,userPrefs);

end %function
     