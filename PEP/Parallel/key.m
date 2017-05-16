%% Hash functions
function out = key(item)
%get row vectors
if size(item,1)~=1
    item = item';
end
%Alt 1:
out = char(100+item);
%Alt 2:
% out=0;
% t=length(item);
% for i = 1:t
%     out = out+10^(t-i)*item(i);
% end
end