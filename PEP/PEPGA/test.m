

k = java.util.ArrayDeque;
%add
a = rand(4,10^5);
for i = 1:10^5
    k.addLast(a(:,i));
end
tic

while ~k.isEmpty
        
            temp = k.removeFirst;
        
        fprintf(1,'%2.0f\n',k.size);
end
toc

