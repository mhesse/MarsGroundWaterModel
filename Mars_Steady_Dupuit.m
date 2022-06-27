% file: Mars_Steady_Dupuit.m
% author: Marc Hesse
% date: Dec 28, 2021
clc, clear, close all

sim_type.heterogeneity_vertical = 'power_law';
% Aquifer constants
const.aq.dmax  = 10e3;   % [m] max. aquifer depth
const.aq.z_bot = -9e3;   % [m] elevation of the base of the aquifer
const.aq.m_exp = 2.5;    % [-] exponent in porosity power-law with height
const.aq.r_exp = 3;      % [-] exponent in porosity-hyd. cond. power-law 
const.aq.phi_ref = 0.3;  % [-] reference porosity for K(phi) power-law
const.aq.K_ref   = 1e-6; % [m/s] reference conductivity for K(phi) power-law
const.aq.phi_s = 0.5;    % [-] actual surface porosity

%% Recharge
r_m_yr = 5e-6;    % 5e-6             % recharge in m/Earth year
% Rest of the constants
const = mars_constants(const);
r_m_s = r_m_yr/const.gen.yr2s; % recharge in m/s

%% Hyd. conductivity
% k = 1e-14
K = 1e-7; % [m/s] 
k = K*const.gen.mu/(const.gen.gMars*const.gen.rho);

% Aquifer properties
[aq_prop,const] = aquifer_properties(const,sim_type);

%% Load topography and precomputed contours
load ../MarsTopoProcessing/Mars_topography.mat
load ../MarsTopoProcessing/hellas_topo.mat
load ../MarsTopoProcessing/argyre_topo.mat
load ../MarsTopoProcessing/dichotomy_topo.mat
load ../MarsTopoProcessing/Mars_1d_topo.mat
disp 'Topo data loaded.'

%% Determine geometry of computational domain 
% (this should be the only place where shorelines are set!)
topo_contours.dichotomy.z = const.sea.Arabia.elev;
% topo_contours.dichotomy.z = const.sea.Deuteronilus.elev;
% topo_contours.dichotomy.z = const.sea.Meridiani.elev;
% topo_contours.hellas.z = -3100; 
topo_contours.hellas.z = -2500; %
% topo_contours.hellas.z = hellas.set.elevations(end-10)
topo_contours.argyre.z = -2500;
% topo_contours.argyre.z = argyre.set.elevations(end-10)
% topo_contours.hellas.z = const.sea.Arabia.elev; %
% topo_contours.argyre.z = const.sea.Arabia.elev;

topo_contours = get_topo_contours(topo_contours,theta,phi,mars_topo,dichotomy,hellas,argyre);

%% Dimless boundary conditions
topo_contours.dichotomy.h = topo_contours.dichotomy.z - const.aq.z_bot;
topo_contours.hellas.h    = topo_contours.hellas.z    - const.aq.z_bot;
topo_contours.argyre.h    = topo_contours.argyre.z    - const.aq.z_bot;

% All shorelines are scaled to the northern ocean (dichotomy)
hc = topo_contours.dichotomy.h;
topo_contours.dichotomy.hD = topo_contours.dichotomy.h/hc; % Dichotomy
topo_contours.hellas.hD    = topo_contours.hellas.h/hc;    % Hellas
topo_contours.argyre.hD    = topo_contours.argyre.h/hc;    % Argyre

%% Dimensionless governing parameters
xc = const.gen.R;
Pi = r_m_s*xc^2/(K*hc^2)

%% Grid and ops
Grid.xmin = 0; Grid.xmax = pi;    Grid.Nx = 25*6; %20;
Grid.ymin = 0; Grid.ymax = 2*pi;  Grid.Ny = 50*6; %50;
Grid = build_grid(Grid);
Dref = build_ops(Grid); % cartesian ops for complex domain comps

Grid.geom = 'spherical_shell_theta_phi';
Grid.R_shell = 1;
Grid = build_grid(Grid);
[Xc,Yc] = meshgrid(Grid.xc,Grid.yc);

% Linear Operators
[D,G,C,I,M] = build_ops(Grid);
L = -D*G;
fs = Pi*ones(Grid.N,1);

% Mean operator (Leave for now)
Nx = Grid.Nx; Ny = Grid.Ny;
Mx = spdiags([ones(Nx,1) ones(Nx,1)]/2,[-1 0],Nx+1,Nx); % 1D mean-matrix in x-dir
My = spdiags([ones(Ny,1) ones(Ny,1)]/2,[-1 0],Ny+1,Ny); % 1D mean-matrix in y-dir
My(1,Ny) = 1/2; My(Ny+1,1) = 1/2;                       % periodicity in y-dir
Mx = kron(Mx,speye(Ny));                                % 2D mean-matrix in x-dir
My = kron(speye(Nx),My);                                % 2D mean-matrix in y-dir
M = [Mx;My];                                            % 2D mean-matrix

% Grid geometry
[dof,dof_f,X_f,Y_f] = get_mars_dofs(Grid,Dref,Xc,Yc,topo_contours);
load mars_dofs.mat
% check_feature_identification(topo_contours,dof,X_f,Y_f,Xc,Yc);

% %% Plot coss-section
% plot(topo1d.theta_deg,topo1d.topo_mean/1e3,'color',.7*[1 1 1],'linewidth',1.5), hold on
% plot([0 180],const.aq.z_bot*[1 1],'k-')
% plot(const.sea.Arabia.theta_bnd_deg,const.sea.Arabia.elev/1e3,'ro')
% % topo_contours.dichotomy.z
% ylabel('elevation: z [km]')
% ylim([-10 5])


%% 2D Residual and Jacobian
H = @(h) spdiags(M*h,0,Grid.Nf,Grid.Nf);  % diagonal matrix with hD ave on faces
dH = @(h) spdiags(G*h,0,Grid.Nf,Grid.Nf); % diagonal matrix with G*hD 
res = @(h) D*(H(h)*G*h) + fs;             % residual vactor
Jac = @(h) D*(H(h)*G+dH(h)*M);            % Jacobian matrix


%% Boundary conditions
BC.dof_dir = [dof.dichotomy_bnd.in;...
              dof.lowlands;...
              dof.hellas_bnd.in;...
              dof.hellas;...
              dof.argyre_bnd.in;...
              dof.argyre];
BC.dof_f_dir = zeros(size(BC.dof_dir)); % not needed
BC.g = [topo_contours.dichotomy.hD*ones(size([dof.dichotomy_bnd.in;dof.lowlands]));...
        topo_contours.hellas.hD*ones(size([dof.hellas_bnd.in;dof.hellas]));...
        topo_contours.argyre.hD*ones(size([dof.argyre_bnd.in;dof.argyre]))];
BC.dof_neu = [];
BC.dof_f_neu = [];
BC.qb = [];
[B,N,fn] = build_bnd(BC,Grid,I);

% Initial guess
fprintf('\nSolve Poisson as initial guess:\n')

hDini = solve_lbvp(L,fs+fn,B,BC.g,N);

figure('name','Initial guess','position',[10 10 1.25*800 1.25*800 ])
subplot 121
h = contourf(Xc,Yc,reshape(hDini,Grid.Ny,Grid.Nx)*10-9,30,'LineColor','none'); hold on
plot(topo_contours.dichotomy.topo.theta,topo_contours.dichotomy.topo.phi,'k-','linewidth',2), hold on
plot(topo_contours.hellas.topo.theta,topo_contours.hellas.topo.phi,'k-','linewidth',2)
plot(topo_contours.argyre.topo.theta,topo_contours.argyre.topo.phi,'k-','linewidth',2)
axis equal
xlim([0 pi]), ylim([0 2*pi])
colorbar

if std(fs)>1e-10 % only if fs is not constant
    subplot 122
    h = contourf(Xc,Yc,reshape(fs,Grid.Ny,Grid.Nx),30,'LineColor','none'); hold on
    plot(topo_contours.dichotomy.topo.theta,topo_contours.dichotomy.topo.phi,'k-','linewidth',2), hold on
    plot(topo_contours.hellas.topo.theta,topo_contours.hellas.topo.phi,'k-','linewidth',2)
    plot(topo_contours.argyre.topo.theta,topo_contours.argyre.topo.phi,'k-','linewidth',2)
    axis equal
    xlim([0 pi]), ylim([0 2*pi])
    colorbar
end

% %% Compute the target zones for drainage basins
target_zones = comp_target_zones(Xc,Yc,Grid,D,G,M,B,BC,N,topo_contours,dof,const,hellas,argyre);

%% Newton iteration
fprintf('\nNewton iteration for Dupuit-Boussinesq:\n')
tol = 1e-9; % convergence tolerance
nmax = 10;   % maximum number of iterations

% Initial guess
hD = hDini; %ones(Grid.N,1);  
hD(BC.dof_dir) = BC.g; % satisfy Dir BC so dhD = 0 on bnd! 
BC.g = BC.g*0;
N_act = length(dof.active);
nres = norm(res(hD)); ndhD = 1; n = 0;
while (nres > tol || ndhD > tol) && n < nmax
    dhD = solve_lbvp(Jac(hD),-res(hD),B,BC.g,N); 
    hD = hD + dhD;
    nres = norm(N'*res(hD))/N_act; ndhD = norm(N'*dhD)/N_act;
    n = n+1;
    fprintf('it = %d: nres = %3.2e  ndhD = %3.2e\n',n,nres,ndhD)
    if n == 1; ndhD = 0; end % to allow exit on first iteration
    nres_Newton(n) = nres; ndhD_Newton(n) = ndhD;
end
fprintf('Newton converged after %d iterations.\n',n)

% Matrix shape for plotting
hDm = reshape(hD,Grid.Ny,Grid.Nx); % [-]  dimensionless head as matrix 
hm  = hDm*hc/1e3;                  % [km] dimensional head as matrix
zm = hm+const.aq.z_bot/1e3;        % [km] elevation of water table
fprintf('Head has been determined.\n\n')

%% Compute fluxes and streamlines
qD = comp_flux_shell(G,1,hD,Grid);
[Vx_c,Vy_c] = comp_cell_center_velocity(qD,Xc,Yc,Grid);

% Compute streamlines on extended domain (periodicity)
[Xc_ext,Yc_ext,Vx_ext] = extend_field_spherical(Xc,Yc,Vx_c,Grid);
[Xc_ext,Yc_ext,Vy_ext] = extend_field_spherical(Xc,Yc,Vy_c,Grid);
fprintf('Computing streamlines\n\n')
S = stream2(Xc_ext,Yc_ext,Vx_ext,Vy_ext,Xc(dof.active),Yc(dof.active));
Ns = length(S);

fprintf('Saving fields.\n')
save('flow_field.mat','Xc','Yc','qD','Vx_c','Vy_c','dof','dof_f','X_f','Y_f','S','Ns','topo_contours','Grid','hD','hDm','hm','zm','fs','target_zones','Xc_ext','Yc_ext','Vx_ext','Vy_ext')

%%
figure('name','Steady Dupuit-Boussinesq','position',[10 10 1.25*800 1.25*800 ])
subplot 121
contourf(rad2deg(Xc),rad2deg(Yc),hDm,50,'LineColor','none'); hold on
plot(rad2deg(topo_contours.dichotomy.topo.theta),rad2deg(topo_contours.dichotomy.topo.phi),'k-','linewidth',2), hold on
plot(rad2deg(topo_contours.hellas.topo.theta),rad2deg(topo_contours.hellas.topo.phi),'k-','linewidth',2)
plot(rad2deg(topo_contours.argyre.topo.theta),rad2deg(topo_contours.argyre.topo.phi),'k-','linewidth',2)
% plot(rad2deg(X_faces),rad2deg(Y_faces),'m-')
colorbar('location','northoutside')
axis equal
xlim([0 180]), ylim([0 360])

title 'dimensionless'

subplot 122
contourf(rad2deg(Xc),rad2deg(Yc),zm,50,'LineColor','none'); hold on
plot(rad2deg(topo_contours.dichotomy.topo.theta),rad2deg(topo_contours.dichotomy.topo.phi),'k-','linewidth',1), hold on
plot(rad2deg(topo_contours.hellas.topo.theta),rad2deg(topo_contours.hellas.topo.phi),'k-','linewidth',1)
plot(rad2deg(topo_contours.argyre.topo.theta),rad2deg(topo_contours.argyre.topo.phi),'k-','linewidth',1)

colorbar('location','northoutside')
axis equal
xlim([0 180]), ylim([0 360])

title 'dimensional'

%%
figure('position',[10 10 1.25*500 1.25*500 ])
plot_spherical_shell(hD,.1,Grid,topo_contours)
view(90,-10)
set(gcf, 'InvertHardCopy', 'off'); 
% set(gcf,'Color',[0 0 0]);
set(gca,'xtick',[],'ytick',[],'ztick',[],'XColor','w','YColor','w','ZColor','w')
set(gca,'Color',[1 1 1]);
grid off
axis tight

saveas(gcf,'WaterOnSphere.png');
