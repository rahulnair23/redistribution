fid = fopen('Runlog.txt');
C = textscan(fid,'%n %n %n %n %n %s','HeaderLines',1);
n = C{1};
p = C{2};
dom = C{3};
PEP_num = C{4};
run_time = C{5};
algo = [C{6}];
isBE = strcmp(algo,'BE');
fail = run_time == -1;

plot(n(~fail & isBE),run_time(~fail & isBE),'xg',...
    n(fail & isBE),0,'xr',...
    n(~fail & ~isBE),run_time(~fail & ~ isBE),'+b',...
    n(fail & ~isBE),0,'or');    


%% PLots
x = 1:30;
y = 1:30;
for i=1:length(x);
    for j = 1:length(y);
        z(i,j) = JointDis([i j],8);
    end
end
surf(x,y,z);