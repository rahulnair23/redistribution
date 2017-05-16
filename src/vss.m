%% Compute the VSS 
clear;
userdata.alpha = 0.9;
Veh = zeros(5,2);
Spaces = zeros(5,2);
seeds = [100:50:10000];
userPrefs.static = false;
for i = 1:20
    clear global;
    userPrefs.seed = seeds(i);
    out = systest(userdata,userPrefs);
    Veh(i,1) = sum(out.balkVeh);
    Spaces(i,1) = sum(out.balkSpace(:));
end

userPrefs.static = true;
for i = 1:20
    clear global;
    userPrefs.seed  = seeds(i);
    out = systest(userdata,userPrefs);
    Veh(i,2) = sum(out.balkVeh(:));
    Spaces(i,2) = sum(out.balkSpace(:));
end
bar(Veh);
figure;
bar(Spaces);
