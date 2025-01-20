% prepare 3D plot
delete(ax)
figure(1)
tiledlayout(1,1)
nexttile
ax = gca;
set(gcf,'Color','w')
% plot measurements
p1 = plot3(meas(1,:),meas(2,:),meas(3,:),'or','MarkerFaceColor','r');
axis equal; hold on
% plot reference
us = 0:0.2:1; thetas = linspace(0,2*pi,1000);
for i = 1:length(us)
    ps = [a_ref*cos(thetas); b_ref*sin(thetas); us(i)*h_ref*ones(1,1000)];
    ps = pos_ref + R_ref*ps;
    p2 = plot3(ps(1,:), ps(2,:), ps(3,:) , 'k', 'LineWidth', 2);
end
% plot estimate UKF
syms theta u
pos = X(1:3,k,j); or = X(4,k,j); 
a = c1(X(5,k,j),0,'lower'); b = c1(X(6,k,j),0,'lower'); h = c1(X(7,k,j),0,'lower');
x = a*cos(theta);
y = b*sin(theta);
z = u*h/2;
p3 = fsurf(x,y,z,[0 2*pi -1 1],'FaceAlpha',0.6,'FaceColor','#4DBEEE');
% translate surf plot
t = hgtransform('Parent',ax); set(p3,'Parent',t);
m = makehgtform('translate',pos); set(t,'Matrix',m)
% rotate surf plot
t = hgtransform('Parent',t); set(p3,'Parent',t);
Rz = makehgtform('zrotate',or); set(t,'Matrix',Rz);
% axes
set(gca,'TickLength',[0 0])
set(gca,'TickLabelInterpreter','latex')
set(gca,'FontSize',12)
xlabel('$x$','Interpreter','latex')
ylabel('$y$','Interpreter','latex')
zlabel('$z$','Interpreter','latex')
% legend
leg = legend([p1,p2,p3],{'Measurement','Reference','Estimate'});
leg.Interpreter = "latex";
leg.Orientation = "horizontal";
leg.NumColumns = 3;
leg.Box = "off";
leg.Layout.Tile = "south";
leg.FontSize = 12;
% draw
drawnow
hold off

% prepare 2D plot
figure(2)
tiledlayout(1,1)
nexttile
set(gcf,'Color','w')
% plot measurements
plot(meas(1,:),meas(2,:),'or','MarkerFaceColor','r');
axis equal; hold on
% plot reference
thetas = linspace(0,2*pi,1000);
ps = [a_ref*cos(thetas); b_ref*sin(thetas)];
ps = pos_ref(1:2) + R_ref(1:2,1:2)*ps;
plot(ps(1,:), ps(2,:), 'k', 'LineWidth', 2);
% plot estimate
pos = X(1:2,k,j); or = X(4,k,j); 
a = c1(X(5,k,j),0,'lower'); b = c1(X(6,k,j),0,'lower');
ps = [a*cos(thetas); b*sin(thetas)];
R_rot = [cos(or) -sin(or); sin(or) cos(or)];
ps = pos + R_rot(1:2,1:2)*ps;
plot(ps(1,:), ps(2,:), 'LineWidth', 2, 'Color', '#4DBEEE');
% axes
set(gca,'TickLength',[0 0])
set(gca,'TickLabelInterpreter','latex')
set(gca,'FontSize',12)
xlabel('$x$','Interpreter','latex')
ylabel('$y$','Interpreter','latex')
set(gca,'Box','off')
xlim([-1.2 5])

% legend
leg = legend('Measurement','Reference','Estimate');
leg.Interpreter = "latex";
leg.Orientation = "horizontal";
leg.NumColumns = 3;
leg.Box = "off";
leg.Layout.Tile = "south";
leg.FontSize = 12;
% draw
drawnow
hold off