function [target_zones] = comp_target_zones(Xc,Yc,Grid,D,G,M,B,BC,N,topo_contours,dof,const,hellas,argyre)
% authors: Marc Hesse, Eric Hiatt
% date: Mar 4, 2022
% Description: This function uses the solution of the poisson equation to
% compute a smoothed topographic contour. This contour is then used
% together with the topographic contour to generate a polygon to determine
% the 'landing region' of a streamline. This is used in determining the
% drainage basins of the topographic features.

% Approximate locations of the saddle points between craters and shoreline
% determined for the thrre different shorelines and high levels in the
% craters.

if topo_contours.dichotomy.z == const.sea.Arabia.elev;
    x_saddle_argyre = 1.121; y_saddle_argyre = 0.7016;
    x_saddle_hellas = 1.539; y_saddle_hellas = 4.828;
    x_saddle_3 = 0.5131;     y_saddle_3 = 6.189;
    x_saddle_4 = 0.5131;     y_saddle_4 = 0.07722;
elseif topo_contours.dichotomy.z == const.sea.Meridiani.elev;
    x_saddle_argyre = 0.9529; y_saddle_argyre = 0.6388;
    x_saddle_hellas = 1.372; y_saddle_hellas = 4.744;
    x_saddle_3 = 0.5131;     y_saddle_3 = 6.189;
    x_saddle_4 = 0.5131;     y_saddle_4 = 0.07722;
elseif topo_contours.dichotomy.z == const.sea.Deuteronilus.elev;
    x_saddle_argyre = 1.435; y_saddle_argyre = 0.4712;
    x_saddle_hellas = 1.707; y_saddle_hellas = 4.66;
    x_saddle_3 = 0.5131;     y_saddle_3 = 6.189;
    x_saddle_4 = 0.5131;     y_saddle_4 = 0.07722;
end

% solve Laplacian with homogeneous Dir BC's
L = -D*G; fs = ones(Grid.N,1);
h = solve_lbvp(L,fs,B,BC.g*0,N);
h_fun = @(pos) interp2(Xc,Yc,reshape(h,Grid.Ny,Grid.Nx),pos(1),pos(2));

%% Find saddle points
% Compute norm of gradient
dh = G*h;
[dhx_c,dhy_c] = comp_cell_center_velocity(dh,Xc,Yc,Grid);
Gh_norm = sqrt(dhx_c.^2+dhy_c.^2);
Gh_norm(dof.inactive) = nan;
Gh_fun = @(pos) interp2(Xc,Yc,reshape(Gh_norm,Grid.Ny,Grid.Nx),pos(1),pos(2));

% Find minima starting from manually selected points
pos_saddle_argyre = fminsearch(Gh_fun,[x_saddle_argyre,y_saddle_argyre]);
h_saddle_argyre = h_fun(pos_saddle_argyre);

pos_saddle_hellas = fminsearch(Gh_fun,[x_saddle_hellas,y_saddle_hellas]);
h_saddle_hellas = h_fun(pos_saddle_hellas);

% Saddle 3 
pos_saddle_3 = fminsearch(Gh_fun,[x_saddle_3,y_saddle_3]);
h_saddle_3 = h_fun(pos_saddle_3);

% Saddle 4 
pos_saddle_4 = fminsearch(Gh_fun,[x_saddle_4,y_saddle_4]);
h_saddle_4 = h_fun(pos_saddle_4);



% Set contour for taget zones for 90% of the lower saddle point
elev = .95*min([h_saddle_argyre,h_saddle_hellas,h_saddle_3,h_saddle_4]);

[cont,n_cont] = get_contours(Grid.xc,Grid.yc,reshape(h,Grid.Ny,Grid.Nx),elev);

if n_cont == 3
    for i=1:n_cont
        if inpolygon(hellas.theta_min,hellas.phi_min,cont(i).x,cont(i).y)
            target_zones.hellas.theta = cont(i).x;
            target_zones.hellas.phi   = cont(i).y;
            target_zones.hellas.n     = length(cont(i).x);
        elseif inpolygon(argyre.theta_min,argyre.phi_min,cont(i).x,cont(i).y)
            target_zones.argyre.theta = cont(i).x;
            target_zones.argyre.phi   = cont(i).y;
            target_zones.argyre.n     = length(cont(i).x);
        else
            target_zones.dichotomy.n     = length(cont(i).x);
            theta = cont(i).x;
            phi = cont(i).y;
            if phi(1) > pi
                theta = flipud(theta);
                phi   = flipud(phi);
            end
            target_zones.dichotomy.theta = [theta;pi;pi;theta(1)];
            target_zones.dichotomy.phi   = [phi;2*pi;0;0];
        end
    end
else
    error('number of taget zones is not 3!\n')
end

figure('name','Target zones','position',[10 10 1.25*800 1.25*800 ])
subplot 131
contourf(Xc,Yc,reshape(h,Grid.Ny,Grid.Nx),1e2,'LineColor','none'); hold on
plot(topo_contours.dichotomy.topo.theta,topo_contours.dichotomy.topo.phi,'k-','linewidth',2), hold on
plot(topo_contours.hellas.topo.theta,topo_contours.hellas.topo.phi,'k-','linewidth',2)
plot(topo_contours.argyre.topo.theta,topo_contours.argyre.topo.phi,'k-','linewidth',2)
colorbar 

subplot 132
% contour(Xc,Yc,reshape(h,Grid.Ny,Grid.Nx),elev*[1 1],'LineColor','r'); hold on
plot(target_zones.hellas.theta,target_zones.hellas.phi,'r-'), hold on
plot(target_zones.argyre.theta,target_zones.argyre.phi,'r-')
plot(target_zones.dichotomy.theta,target_zones.dichotomy.phi,'r-')
plot(target_zones.dichotomy.theta(1),target_zones.dichotomy.phi(1),'ro')
plot(x_saddle_argyre,y_saddle_argyre,'r*')
plot(pos_saddle_argyre(1),pos_saddle_argyre(2),'bo')
plot(x_saddle_hellas,y_saddle_hellas,'r*')
plot(pos_saddle_hellas(1),pos_saddle_hellas(2),'bo')
plot(x_saddle_3,y_saddle_3,'r*')
plot(pos_saddle_3(1),pos_saddle_3(2),'bo')
plot(x_saddle_4,y_saddle_4,'r*')
plot(pos_saddle_4(1),pos_saddle_4(2),'bo')
plot(topo_contours.dichotomy.topo.theta,topo_contours.dichotomy.topo.phi,'k-','linewidth',2), hold on
plot(topo_contours.hellas.topo.theta,topo_contours.hellas.topo.phi,'k-','linewidth',2)
plot(topo_contours.argyre.topo.theta,topo_contours.argyre.topo.phi,'k-','linewidth',2)

subplot 133
contourf(Xc,Yc,reshape(Gh_norm,Grid.Ny,Grid.Nx),1e2); hold on
