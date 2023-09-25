function plotbound(xy,du,color)
% plot the values of the boundary on the unit circle dependent on
% the angle
% 'xy' is a 2xn matrix containing the coordinates on the unit circle
% 'du' is a 1xn matrix with the corresponding boundary values

if nargin<3
    color = {'k','r','b','g','y','c','k','k','k','k'};
end
n=length(du);
a=zeros(n,1);
for i=1:n
    x=xy(1,i);
    y=xy(2,i);
    if x>0 && y>=0
        a(i)=atan(y/x);
    elseif x<=0 && y>=0
        a(i)=pi./2 + atan(-x/y);
    elseif x<=0 && y<0
        a(i)=pi+atan(y/x);
    elseif x>0 && y<0
        a(i)=(3*pi/2) + atan(-x/y);
    end
end
[A,index]=sort(a);
axes('XTickLabel',{'0','pi/2','pi','3*pi/2','2*pi'},...
    'XTick',[0 1.571 3.142 4.712 6.283]);
xlim([0 6.283]);
box('on');
Neig = size(du,2);
for j = 1:Neig
    hold on
    plot(A,du(index,j),color{j});
    hold off
end
