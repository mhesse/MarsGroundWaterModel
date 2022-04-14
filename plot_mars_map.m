function [h] = plot_mars_map(Xc,Yc,F,Grid,topo_contours)
% author: Marc Hesse
% date: Dec 10, 2021
% description: takes whatever we have from computations and orients it the
% way Martians are used to see it.
Phi_half = [Yc(Grid.Ny/2+1:end,:);Yc(1:Grid.Ny/2,:)+2*pi]-2*pi;
Theta_half = [Xc(Grid.Ny/2+1:end,:);Xc(1:Grid.Ny/2,:)];
F_half = [F(Grid.Ny/2+1:end,:);F(1:Grid.Ny/2,:)];

if max(F(:)) ~= 0
    h = surf(Theta_half,Phi_half,F_half); shading interp, hold on
end
theta = topo_contours.dichotomy.topo.theta;
phi   = topo_contours.dichotomy.topo.phi;
elev = 4*ones(size(topo_contours.dichotomy.topo.phi));
phi_half = [phi(phi<=pi)];
theta_half = [theta(phi<=pi)];
elev_half = [elev(phi<=pi)];

%% Dichotomy boundary
phi_2nd_half = phi(phi>pi)-2*pi;
theta_2nd_half = theta(phi>pi);
elev_2nd_half = elev(phi>pi);
if max(F(:)) == 0
    h = plot3(theta_half,phi_half,elev_half,'k-','linewidth',2), hold on
end
plot3(theta_half,phi_half,elev_half,'k-','linewidth',2), hold on
plot3(theta_2nd_half,phi_2nd_half,elev_2nd_half,'k-','linewidth',2), hold on


%% Hellas
if isfield(topo_contours.hellas,'topo')
    elev_hellas = 4*ones(size(topo_contours.hellas.topo.theta));
    plot3(topo_contours.hellas.topo.theta,topo_contours.hellas.topo.phi-2*pi,elev_hellas,'k-','linewidth',2)
end
%% Argyre
if isfield(topo_contours.argyre,'topo')
    elev_argyre = 4*ones(size(topo_contours.argyre.topo.theta));
    plot3(topo_contours.argyre.topo.theta,topo_contours.argyre.topo.phi,elev_argyre,'k-','linewidth',2)
end
plot([0 pi pi 0],[-pi -pi pi pi],'k-')
axis equal
view([-90 90])
xlim([0 pi]), ylim([-pi pi])
set(gca,'Xtick',[0 pi/4 pi/2 3*pi/4 pi],'Ytick',[-pi,-pi/2,0,pi/2,pi],'Xticklabel',[-90 -45 0 45 90],'Yticklabel',[-180 -90 0 90 180])

xlabel 'lattitude [\circ]'
ylabel 'longitude [\circ]'