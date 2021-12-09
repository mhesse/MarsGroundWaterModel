function [] = plot_spherical_shell(hD,scale,Grid)
% Plot sphere
% [Theta_s,Phi_s] = meshgrid(linspace(0,pi,100),linspace(0,2*pi,100));
% 
% Xs = R*sin(Theta_s).*cos(Phi_s);
% Ys = R*sin(Theta_s).*sin(Phi_s);
% Zs = R*cos(Theta_s);

% solution
R = 1;
[Theta,Phi] = meshgrid([0;Grid.xc;pi],[Grid.yc;Grid.yc(1)]);
Hd = reshape(hD,Grid.Ny,Grid.Nx);
Hd = [Hd;Hd(1,:)];
Hd = [Hd(:,1) Hd Hd(:,end)];
Xh = (R+scale*Hd).*sin(pi-Theta).*cos(Phi);
Yh = (R+scale*Hd).*sin(pi-Theta).*sin(Phi);
Zh = (R+scale*Hd).*cos(pi-Theta);

% surf(Xs,Ys,Zs), hold on

s = surf(Xh,Yh,Zh);
s.CData = Hd;
shading interp
axis equal
end

