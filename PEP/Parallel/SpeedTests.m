function SpeedTests 

%p = [0.8 0.9 0.95];
%n = [4 6 8 10 12 16];
p=0.95;
n=[8 10 12 16];

relay('SysTime \t n \t p \t #PEP \t Time(be) \t Algo \n'); 
    for j=1:length(n)

     %generate inputs
     nodes = n(j);
     userPrefs.parameters = {ones(nodes,1), 2*ones(nodes,1)};
     userPrefs.type = logical([ones(1,nodes/2) zeros(1,nodes/2)]);
     userPrefs.verbose = true;
     userPrefs.maxtime = 10800;
     userPrefs.algorithm = 'completeenumeration';
     
     LB = -10*ones(nodes,1);
     UB = 10*ones(nodes,1);

          for i =1:length(p)
              userPrefs.p = p(i);
             try
                 %run 
                     %s=tic;
                        out = PEP_rd(LB,UB,userPrefs);
                     %dt = toc(s);
                  relay('%s \t %2.0f \t %0.2f \t %5.0f \t %5.3f \t RD',...
                        datestr(now),n(j),p(i),size(out.S,1),out.runtime);
                 relay('\n');
             catch ME
                relay('%s \t %2.0f \t %0.2f \t RD \t',datestr(now),n(j),p(i));
                relay(ME.message);
                relay('\n');
             end
         end 
     end 
 end 