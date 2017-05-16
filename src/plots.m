%% Plot one solution
a=zeros(15,15);
a(1,9)=3;
a(2,11)=2;
a(5,13)=1;
a(11,14)=1;
a(9,15)=2;
a=a(2:end,2:end);
load Singapore.mat;
hold on;
plot(LatLong(:,2),LatLong(:,1),'ob','MarkerFaceColor','b');
box('on');

for i=1:14
    for j = 1:14
        if a(i,j)>0
             line(LatLong([i j],2)',LatLong([i j],1)','LineWidth',a(i,j));
        end
    end
end

%% Determine lambda2 and do the Chi-squared test
periods = [ 0 9 12 18 24]; 
load checkouts.mat;
load returns.mat;
load Singapore.mat;
x = 0:15;
SelectedPorts = 1:14;
lambda = zeros(length(SelectedPorts),length(periods)-1);
chitest = false(length(SelectedPorts),length(periods)-1); %booleans
out = zeros(length(x),length(periods)-1,length(SelectedPorts));
for p = SelectedPorts
        for t = 1:length(periods)-1
            b = unique(RetDayID(Retport==p));
            temp = zeros(size(b));
            for i  = 1:length(temp)
                temp(i) = sum(Retport==p & RetDayID==b(i) & RetHour>=periods(t) & RetHour<periods(t+1));
            end
            %fit poisson num
            lambda(p==SelectedPorts,t) = poissfit(temp);
            %store the test
            chitest(p==SelectedPorts,t) = ...
                chi2gof(temp,...
                'cdf',{@poisscdf,lambda(p==SelectedPorts,t)},...
                'nbins',5);
        end
end %p
save('lambda2.mat','lambda');
[r,c] = find(chitest);
if ~isempty(r)
    fprintf('These failed the Chi squared test\n');
   for i = 1:length(r)
        fprintf('Port: %3.0f, Time Period: %3.0f \n',...
            SelectedPorts(r(i)),c(i));
   end
end %if

%% Determine lambda1 and do the Chi-squared test
periods = [ 0 9 12 18 24]; 
load checkouts.mat;
load Singapore.mat;
x = 0:15;
SelectedPorts = 1:14;
lambda = zeros(length(SelectedPorts),length(periods)-1);
chitest = false(length(SelectedPorts),length(periods)-1); %booleans
out = zeros(length(x),length(periods)-1,length(SelectedPorts));
for p = SelectedPorts
        for t = 1:length(periods)-1
            b = unique(dayID(port==p));
            temp = zeros(size(b));
            for i  = 1:length(temp)
                temp(i) = sum(port==p & dayID==b(i) & hour>=periods(t) & hour<periods(t+1));
            end
            %fit poisson num
            lambda(p==SelectedPorts,t) = poissfit(temp);
            %store the test
            chitest(p==SelectedPorts,t) = ...
                chi2gof(temp,...
                'cdf',{@poisscdf,lambda(p==SelectedPorts,t)},...
                'nbins',5);
        end
end %p
save('lambda1.mat','lambda');
[r,c] = find(chitest);
if ~isempty(r)
    fprintf('These failed the Chi squared test\n');
   for i = 1:length(r)
        fprintf('Port: %3.0f, Time Period: %3.0f \n',...
            SelectedPorts(r(i)),c(i));
   end
end %if

%% Distribution of demand imbalance by time period
load checkouts.mat;
load Singapore.mat;
load returns.mat;

SelectedPorts =1:14;
periods = [ 0 9 12 18 24]; 
NumPorts = length(SelectedPorts);
x = -20:20;

for p = SelectedPorts
    %initialize
    out = zeros(length(x),length(periods)-1);
    %get the days the port operated for
    c = sort(unique(dayID(port==p)));
    operated(p) = c(end)-c(1);
    
    for t = 1:length(periods)-1
        temp = zeros(operated(p),1);
        for i = 1:length(c)
            checkouts = sum(port==p & dayID==c(i) & hour>=periods(t) & hour<periods(t+1) ); 
            returns = sum(Retport==p & RetDayID==c(i) & hour>=periods(t) & hour<periods(t+1)); 
            temp(i) = checkouts-returns;
        end %for
        out(:,t) = hist(temp,x);
    end %t
    %plot the distribution
    figure;
    hold all;
    plot(x,out./repmat(sum(out),length(x),1),'-x');
      for i = 1:length(periods)-1
            leg{i} = [num2str(periods(i)) '<= t < ' num2str(periods(i+1))];
      end %i
    legend(leg);
    box('on');
    xlim([-10 10]);
    title(['Distribution of Net Flow by time of day for Port (' num2str(p) ') ' names{p}]);
    xlabel('Net Flow in Vehicles');
    ylabel('Probability');
end %p

%% Distribution of Daily Imbalance
SelectedPorts =2:14;
NumPorts = length(SelectedPorts);
x = -20:20;
out = zeros(length(x),NumPorts);

load checkouts.mat;
load Singapore.mat;
load returns.mat;

for p = SelectedPorts
    c = sort(unique(dayID(port==p)));
    operated(p) = c(end)-c(1);
    temp = zeros(operated(p),1);
    for i = 1:length(c)
        checkouts = sum(port==p & dayID==c(i)); %number of checkouts 
        returns = sum(Retport==p & RetDayID==c(i)); 
        temp(i) = checkouts-returns;
    end %for
    out(:,SelectedPorts==p) = hist(temp,x);
end %p

%plot
plot(x,out./repmat(operated(SelectedPorts),length(x),1),'-x');
legend(names(SelectedPorts));
xlabel('Aggregate Flow (Checkouts-Returns)');
ylabel('Probability');
title('Distribution of Daily imbalance');
xlim([-10 10]);
ylim([0 0.6]);

%% Fleet imbalance by day for each station
SelectedPorts = 1:14;
NumPorts = length(SelectedPorts);
out = zeros(NumPorts,3);

load checkouts.mat;
load Singapore.mat;
load returns.mat;
out = zeros(NumPorts,3); %Numdays with +checkouts,0,-checkouts
for p = SelectedPorts
    c = unique(dayID(port==p));
    for i = 1:length(c)
        checkouts = sum(port==p & dayID==c(i)); %number of checkouts 
        returns = sum(Retport==p & RetDayID==c(i)); 
        switch true
            case checkouts-returns > 0
                out(SelectedPorts==p,1) = out(SelectedPorts==p,1)+1;
            case checkouts-returns ==0
                out(SelectedPorts==p,2) = out(SelectedPorts==p,2)+1;
            case checkouts-returns<0
                out(SelectedPorts==p,3) = out(SelectedPorts==p,3)-1;
        end %switch
    end %for
end %p

%plot
b1 = bar(out(:,2),'FaceColor','g');
set(b1,'BarWidth',1);
hold on;
bar(out(:,1),'b');
bar(out(:,3),'r');
legend('Balanced','More checkouts','More returns');
xlabel('Station');
set(gca,'XTickLabel',SelectedPorts);
ylabel('Number of Days');
title('Daily Fleet Imbalance');
%make new tick labels
th=text(1:14,zeros(14,1),names,'HorizontalAlignment','left',...
    'rotation',90,...
    'BackgroundColor','w',...
    'FontSize',8);

%% Plot maps of capacity and fleet position
load Singapore.mat;
h=bar([capacity inventory],'EdgeColor','none');
colormap winter;
legend('capacity','base inventory');
set(gca,'xticklabel',{});
xlabel('Stations');
ylabel('Number of Vehicles');

%make new tick labels
th=text(1:14,0.5*ones(14,1),names,'HorizontalAlignment','left',...
    'rotation',90,...
    'BackgroundColor','w',...
    'FontSize',9);

%% Plot map of volumes
load Singapore.mat;
S = shaperead('H:\Rahul\Research\Sharing\data\maps\singaporeWGS1984');
mapshow(S);
hold on;
plot(LatLong(:,2),LatLong(:,1),'ob','MarkerFaceColor','b');
box('on');
%text(LatLong(:,2),LatLong(:,1),names);


%Compute the total flows
load checkouts.mat;
load returns.mat;
NumPorts = 14;
out = zeros(NumPorts,NumPorts);

for i = 1:length(port)
    out(port(i),Retport(i)) = out(port(i),Retport(i)) + 1;
end
%Normalize to get Per day values.
out = out./daysOpenOD;

%output for processing
for i = 1:14
    for j = 1:14
        if out(i,j)>0
    fprintf(1,'%s,%s,%3.4f\n',names{i},names{j},out(i,j));
        end
    end
end

for i = 1:NumPorts
    for j = 1:NumPorts
        if i~=j
            if out(i,j)>=0.01
            line(LatLong([i j],2)',LatLong([i j],1)','LineWidth',10*out(i,j));
            end
        end
    end
end
for i=1:14
    text(LatLong(i,2)+0.003,LatLong(i,1),num2str(i));
end
xlim([103.8 104.02]);
ylim([1.24 1.4]);

%% Distribution of all checkouts for each station
load checkouts.mat;
load Singapore.mat;
x = [0:20];
selectedPorts = 3:8;
out = zeros(length(x),length(selectedPorts));
for p = selectedPorts
    temp = zeros(daysOpen(p),1);
    
        for i  = 1:length(temp)
            temp(i) = sum(port==p & dayID==i);
        end
     out(:,selectedPorts==p) = hist(temp,x);
end 
plot(x,out./repmat(daysOpen(selectedPorts)',size(out,1),1));
legend(names(selectedPorts));
title('Distribution of Checkouts per Day');
xlabel('Number of Checkouts (k)');
ylabel('P(Checkouts = k)');

%% Distribution by time periods for specific station

periods = [ 0 9 12 18 24]; 
load checkouts.mat;
load Singapore.mat;
x = 0:20;
SelectedPorts = 1:14;
lambda = zeros(length(SelectedPorts),length(periods)-1);

for p = SelectedPorts
    out = zeros(length(x),length(periods)-1);
        for t = 1:length(periods)-1
            b = unique(dayID(port==p));
            temp = zeros(size(b));
            for i  = 1:length(temp)
                temp(i) = sum(port==p & dayID==b(i) & hour>=periods(t) & hour<periods(t+1));
            end
            %fit poisson num
            lambda(p,t) = poissfit(temp);
            out(:,t) = hist(temp,x);
        end
    figure;
    hold all;
    box('on');
    plot(x,out./daysOpen(p));
    title(['Distribution of demand by time of day for Port (' num2str(p) ') ' names{p}]);
    xlabel('Number of Checkouts');
    ylabel('Probability');
        for i = 1:length(periods)-1
            leg{i} = [num2str(periods(i)) '<= t < ' num2str(periods(i+1))];
        end %i
    legend(leg);
    
    %add theoretical curve
    for t = 1:length(periods)-1
        %calc the actual value
        y = poisspdf(x,lambda(p,t));
        %plot
        plot(x,y,':k');
    end %t
end %p


%% Trip duration and distance (probability plots)
load checkouts.mat;
probplot('exponential',usetime/60);
box('on');
title ('Distribution of trip duration');
xlabel('Trip duration in hours');
set(gca,'XTick',[0:24:100 200:100:500]);
grid('on');

figure;
probplot('exponential',usedistance);
box('on');
title('Distribution of trip distance');
xlabel('Trip length in kilometers');

%% Average checkouts by hour of day for each port
load checkouts.mat;
load Singapore.mat;
NumPorts = 14;
out = zeros(24,NumPorts);
for p = 1:NumPorts
    for h = 0:23
        out(h+1,p) = sum(hour==h & port==p)/ daysOpen(p);
    end
end
plot(out);
legend(names,'location','Northwest');
title('Average number of checkouts by hour');
xlabel('Hour of day');
set(gca,'XTick',[0:3:24],'XTickLabel',[0:3:24]);

%% Chi2gof test
x = poissrnd(5,1000,1);
lambda = mean(x);
k = chi2gof(x,'cdf',{@poisscdf,15});
if k
    fprintf(1,'Reject\n');
else
    fprintf(1,'Accept\n');
end

%% Set up data inputs
load Singapore.mat;
load lambda.mat;

data.nodes = 15; %includes dummy node
data.states = 4;
%append the dummy node
lambda1 = horzcat(zeros(data.states,1),lambda1');
lambda2 = horzcat(zeros(data.states,1),lambda2');
data.lambda1=lambda1;
data.lambda2=lambda2;
cost = 1000*ones(data.nodes);
cost(2:end,2:end) = duration/100;
data.cost = cost;
data.capacity = [100;capacity];
data.inventory = floor(data.capacity/2);
data.alpha = 0.9; %desired level
data.periods = 30;
data.cyclic = true;

save('data.mat','data');



