close all
rand('state',0)
sdof = fem.sdof;
x = femf.gcoord(:,1);
y = femf.gcoord(:,2);
sigmatt = eval(sigmat);
sigmatt(1:sdof) = sigmatt(1:sdof) + (2*rand(sdof,1)-1)*.2;
pdeplot(p1,e1,t1,'xydata',sigmatt,'contour','off','colorbar','off'), colormap('jet')
colorbar
axis off, axis image
s=linspace(0,2*pi,100);
hold on
plot(0.3*cos(s)+.5,0.3*sin(s)+0.45,'k','linewidth',2)
hold off
set(colorbar,'FontSize',16), drawnow
print(1, '-dpng',['exp6-ex'])