function SpeedTests 

p=[0.5 0.6 0.7 0.8 0.9];
n=[4 5 8 10 15 20];
domain = [5 10 15 20 30 50];

relay('SysTime \t n \t p \t Domain \t #PEP \t Time(be) \t Algo \n'); 
    for j=1:length(n)

        for k = 1:length(domain)
            
         %generate inputs
         nodes = n(j);
         userPrefs.type = true(1,nodes);
         userPrefs.verbose = true;
         userPrefs.maxtime = 10800;
         userPrefs.algorithm = 'backwardenumeration';
         userPrefs.jointDistribution = @JointDis;
         userPrefs.parameters = domain(k);
         userPrefs.precompute = false;
         userPrefs.WritePEPtoFile = false;

         LB = ones(nodes,1);
         UB = domain(k)*ones(nodes,1);

              for i =1:length(p)
                 userPrefs.p = p(i);
                 %Backward enumeration
                 try
                     %run 
                    userPrefs.minsize = 10^50; %for enumerate
                    out = PEP_rd(LB,UB,userPrefs);
                    relay('%s \t %2.0f \t %0.2f \t %5.0f \t %5.0f \t %5.3f \t BE',...
                            datestr(now),n(j),p(i),domain(k),size(out.S,1),out.runtime);
                     relay('\n');
                  catch ME
                     relay('%s \t %2.0f \t %0.2f \t RD \t %5.0f \t',datestr(now),n(j),p(i),domain(k));
                     relay(ME.message);
                     relay('\n');
                 end
                 
                 %splitting algorithm
                 try
                     %run 
                    userPrefs.minsize = min(1000,2^nodes+1); 
                    out = PEP_rd(LB,UB,userPrefs);
                    relay('%s \t %2.0f \t %0.2f \t %5.0f \t %5.0f \t %5.3f \t RD',...
                            datestr(now),n(j),p(i),domain(k),size(out.S,1),out.runtime);
                     relay('\n');
                  catch ME
                      relay('%s \t %2.0f \t %0.2f \t RD \t %5.0f \t',datestr(now),n(j),p(i),domain(k));
                     relay(ME.message);
                     relay('\n');
                 end
                 
              end %for i
        end %for k
     end %for j
 end 