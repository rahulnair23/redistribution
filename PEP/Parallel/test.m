%%parallel comp
matlabpool open
parfor i=1:10
    WritePEP([i i],10);
end
%% Test duplicity
k = java.util.Hashtable;
k.put('1',1);
k.put('1',2);
k.put('1',3);
k.put('2',3);


%% .NET Dictionary test
currGen = NET.createGeneric('System.Collections.Generic.Dictionary',...
  {'System.String','System.String'},10^7);
for j=1:10
    for i = 1:100
        currGen.Add(char(i),'');
    end
    nextGen = NET.createGeneric('System.Collections.Generic.Dictionary',...
  {'System.String','System.String'},10^7);
    currGen=nextGen;
end
    
    

%% More java tests
tic
currGen = java.util.Hashtable(2^20);
for i =1:10
    
    %do something with currGen
    for j = 1:10000
        currGen.put(char(j),j);
    end
    
    %create next gen from currGen
    nextGen = java.util.Hashtable(2^20);
    
    %Update generation
    currGen=0;
    currGen=nextGen;
end
toc

%% minus vs bsxfun

A=1:20;
tic
for i = 1:10000
    k=bsxfun(@minus,A(1:10),A(11:end));
end
toc
tic
for i = 1:10000
    k=A(1:10)-A(11:end);
end
toc
%% prod vs. bsxfun and loop
A=1:20;
tic
for i =1:500000
    k=prod(A);
end
toc

tic
for i = 1:500000
    k=cumprod(A);
    k=k(end);
end
toc

tic
for i = 1:500000
    for j = 1:20
        k=k*A(j);
    end
end
toc
%% are java objects value or handle classes
k1=java.util.Hashtable;
k1.put('1',1);
k1.put('2',2);
k2=k1.clone;
k2.remove('2');

k1.size
%So java classes are handle classes


%% Test assignment
len =10^6;
tic
k = zeros(len,30);
toc
tic
k(1:len,1:30)=0;
toc
%% Demonstrates BSXFUN
len=100000;
    nodes=28;
    itr=100;
    A = 0.5*ones(nodes,1);
    B = rand(nodes,len);
    
    tic;
        for i=1:itr
            k=repmat(A,1,len)>B;
        end
    toc

    tic;
        for i=1:itr
            k=bsxfun(@gt,A,B);
        end
    toc

    