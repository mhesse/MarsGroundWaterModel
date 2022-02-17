clear, close all, clc

% Read in the file.
load ../MarsTopoProcessing/Mars_topography.mat
load ../MarsTopoProcessing/hellas_topo.mat
load ../MarsTopoProcessing/argyre_topo.mat
load ../MarsTopoProcessing/dichotomy_topo.mat
disp 'Topo data loaded.'

% el = flipud(el);


%% Load Ed Kites recharge data
cd ../MOLA_Steady
load bestfitfanera.mat
load bestfitvnera_v2.mat
recharge = bestfitfanera(1:2:end,1:2:end);
% recharge = bestfitvnera(1:2:end,1:2:end);
rech_kite  = flipud(recharge');
cd ../MarsDrainageBasins

%% Extract contours
% is we change these 
topo_contours.dichotomy.z = -3800;
topo_contours.hellas.z = -3100; %-5800; %
topo_contours.argyre.z = -2500;

topo_contours = get_topo_contours(topo_contours,theta,phi,mars_topo,dichotomy,hellas,argyre);

%% Generating the computational domain
Grid.xmin = 0; Grid.xmax = pi;   Grid.Nx = 25*6;
Grid.ymin = 0; Grid.ymax = 2*pi; Grid.Ny = 50*6;
Grid = build_grid(Grid);

% For complex domain computations we need catesian Dref
[Dref,Gref,C,I,M] = build_ops(Grid);

Grid.geom = 'spherical_shell_theta_phi';
Grid.R_shell = 1;
Grid = build_grid(Grid);

% Operators
[D,G,C,I,M] = build_ops(Grid);
[Xc,Yc] = meshgrid(Grid.xc,Grid.yc);

[dof,dof_f,X_f,Y_f] = get_mars_dofs(Grid,Dref,Xc,Yc,topo_contours);
% load mars_dofs.mat
check_feature_identification(topo_contours,dof,X_f,Y_f,Xc,Yc);

%% Try to solve Laplace - as is
% D = Dref; G = Gref;
G(dof.dir,:) = 0;
L = -D*G;
% fs = ones(Grid.N,1);
rech_interp = interp2(Theta,Phi,rech_kite,Xc,Yc);
fs = 1*rech_interp(:);

% only impose Dir BC's on Dichtonomy boundary
% dof_dir = setdiff(dof_bnd_in,[Grid.dof_xmin;Grid.dof_ymin;Grid.dof_ymax]);

% Boundary conditions
BC.dof_dir = [dof.dichotomy_bnd.in;...
              dof.lowlands;...
              dof.hellas_bnd.in;...
              dof.hellas;...
              dof.argyre_bnd.in;...
              dof.argyre];

% heads are sort of non-dimensionalized: top = 1km base -9 km
% hellas (-5.8km) is 3.2km above base so h' = 3.2/10 = 0.32
BC.g = [0.52 *ones(length([dof.dichotomy_bnd.in;dof.lowlands]),1);...
        0.59*ones(length([dof.hellas_bnd.in;dof.hellas]),1);...
        0.65*ones(length([dof.argyre_bnd.in;dof.argyre]),1)];

[B,N,fn] = build_bnd(BC,Grid,I);

u = solve_lbvp(L,fs,B,BC.g,N);
q = -G*u;
[Vx_c,Vy_c] = comp_cell_center_velocity(q,Xc,Yc,Grid);
S = stream2(Xc,Yc,Vx_c,Vy_c,Xc(dof.active),Yc(dof.active));
Ns = length(S);
save('flow_field.mat','Xc','Yc','Vx_c','Vy_c','dof','dof_f','X_f','Y_f','S','Ns','topo_contours','Grid','u','fs')

%%
figure('position',[10 10 1.25*800 1.25*800 ])
subplot 121
h = contourf(Xc,Yc,reshape(u,Grid.Ny,Grid.Nx)*10-9,30,'LineColor','none'); hold on
plot(topo_contours.dichotomy.topo.theta,topo_contours.dichotomy.topo.phi,'k-','linewidth',2), hold on
plot(topo_contours.hellas.topo.theta,topo_contours.hellas.topo.phi,'k-','linewidth',2)
plot(topo_contours.argyre.topo.theta,topo_contours.argyre.topo.phi,'k-','linewidth',2)
axis equal
xlim([0 pi]), ylim([0 2*pi])
colorbar
subplot 122
h = contourf(Xc,Yc,reshape(fs,Grid.Ny,Grid.Nx),30,'LineColor','none'); hold on
plot(topo_contours.dichotomy.topo.theta,topo_contours.dichotomy.topo.phi,'k-','linewidth',2), hold on
plot(topo_contours.hellas.topo.theta,topo_contours.hellas.topo.phi,'k-','linewidth',2)
plot(topo_contours.argyre.topo.theta,topo_contours.argyre.topo.phi,'k-','linewidth',2)
axis equal
xlim([0 pi]), ylim([0 2*pi])
colorbar
%%
figure('position',[10 10 1.25*500 1.25*500 ])
plot_spherical_shell(u,.1,Grid)
view(90,-10)
set(gcf, 'InvertHardCopy', 'off'); 

% set(gca,'visible','off')
% saveas(gcf,'WaterOnSphere.png');
