function skellamCDFplot(x,y,col)
%SKELLAMPLOT Summary of this function goes here
%   Detailed explanation goes here

if nargin == 2
    col = 'r';
end
%determine axis
for i=1:length(x)
    if y(i)>0.0000001
        x_min = x(i-1);
        break
    end
end
x_max = -x_min;

xx = x(x>x_min & x<x_max);
yy = y(x>x_min & x<x_max);

hold all;
for i = 1:length(xx)-1
    line([xx(i) xx(i+1)],[yy(i) yy(i)],'LineStyle','-','Color',col);
    plot(xx(i),yy(i),'or','MarkerEdgeColor',col,'MarkerFaceColor',col,'MarkerSize',4);
    plot(xx(i+1),yy(i),'o','MarkerSize',2,'MarkerEdgeColor',col);
end

end %function
