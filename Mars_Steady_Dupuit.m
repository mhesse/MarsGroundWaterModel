% file: HomogeneousUniformRecharge2D.m
% author: Marc Hesse
% date: Dec 28, 2021

clc, clear, close all

sim_type.heterogeneity_vertical = 'power_law';
% Aquifer constants
const.aq.dmax  = 10e3; % [m] max. aquifer depth
const.aq.m_exp = 2.5;  % [-] exponent in porosity power-law with height
const.aq.r_exp = 3;    % [-] exponent in porosity-hyd. cond. power-law 
const.aq.phi_ref = 0.3;  % [-] reference porosity for K(phi) power-law
const.aq.K_ref   = 1e-6; % [m/s] reference conductivity for K(phi) power-law
const.aq.phi_s = 0.5; % [-] actual surface porosity

%% Recharge
r_m_yr = 5e-6;                 % recharge in m/Earth year
% Rest of the constants
const = mars_constants(const);
r_m_s = r_m_yr/const.gen.yr2s; % recharge in m/s

%% Hyd. conductivity
% k = 1e-14
K = 1e-7; % [m/s] 
k = K*const.gen.mu/(const.gen.gMars*const.gen.rho);

% Aquifer properties
[aq_prop,const] = aquifer_properties(const,sim_type);

%% Scales and dimless. params
hc = const.sea.Deuteronilus.z;
xc = const.gen.R;
Pi = r_m_s*xc^2/(K*hc^2)

%% Geometry
topo_contours.dichotomy.z = const.sea.Arabia.elev;
topo_contours.hellas.z = -3100; 
topo_contours.hellas.z = -5800; %
% topo_contours.argyre.z = -2500;
% topo_contours.hellas.z = const.sea.Arabia.elev; %
topo_contours.argyre.z = const.sea.Arabia.elev;
load ../MarsTopoProcessing/Mars_topography.mat
load ../MarsTopoProcessing/hellas_topo.mat
load ../MarsTopoProcessing/argyre_topo.mat
load ../MarsTopoProcessing/dichotomy_topo.mat
disp 'Topo data loaded.'

topo_contours = get_topo_contours(topo_contours,theta,phi,mars_topo,dichotomy,hellas,argyre);

%% Dimless boundary conditions
h_d = topo_contours.dichotomy.z - (-const.aq.dmax);
h_h = topo_contours.hellas.z - (-const.aq.dmax);
h_a = topo_contours.argyre.z - (-const.aq.dmax);
hc = h_d;
hD_d = h_d/hc;  % Dichotomy
hD_h = h_h/hc; % Hellas
hD_a = h_a/hc; % Argyre

%% Grid and ops
Grid.xmin = 0; Grid.xmax = pi;    Grid.Nx = 25*6; %20;
Grid.ymin = 0; Grid.ymax = 2*pi;  Grid.Ny = 50*6; %50;
Grid = build_grid(Grid);
[Dref,Gref,C,I,M] = build_ops(Grid); % cartesian ops for complex domain comps

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
% load mars_dofs.mat
check_feature_identification(topo_contours,dof,X_f,Y_f,Xc,Yc);

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
BC.g = [hD_d*ones(size([dof.dichotomy_bnd.in;dof.lowlands]));...
        hD_h*ones(size([dof.hellas_bnd.in;dof.hellas]));...
        hD_a*ones(size([dof.argyre_bnd.in;dof.argyre]))];
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
subplot 122
h = contourf(Xc,Yc,reshape(fs,Grid.Ny,Grid.Nx),30,'LineColor','none'); hold on
plot(topo_contours.dichotomy.topo.theta,topo_contours.dichotomy.topo.phi,'k-','linewidth',2), hold on
plot(topo_contours.hellas.topo.theta,topo_contours.hellas.topo.phi,'k-','linewidth',2)
plot(topo_contours.argyre.topo.theta,topo_contours.argyre.topo.phi,'k-','linewidth',2)
axis equal
xlim([0 pi]), ylim([0 2*pi])
colorbar

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

% matrix shape for plotting
hDm = reshape(hD,Grid.Ny,Grid.Nx);
hm = (hDm*hc-const.aq.dmax)/1e3;
q = -G*hD;
[Vx_c,Vy_c] = comp_cell_center_velocity(q,Xc,Yc,Grid);
S = stream2(Xc,Yc,Vx_c,Vy_c,Xc(dof.active),Yc(dof.active));
Ns = length(S);
save('flow_field_with_crater_5800.mat','Xc','Yc','Vx_c','Vy_c','dof','dof_f','X_f','Y_f','S','Ns','topo_contours','Grid','hD','hDm','hm','fs')

%%
figure('name','Steady Dupuit-Boussinesq','position',[10 10 1.25*800 1.25*800 ])
subplot 121
contourf(rad2deg(Xc),rad2deg(Yc),hDm,50,'LineColor','none'); hold on
plot(rad2deg(topo_contours.dichotomy.topo.theta),rad2deg(topo_contours.dichotomy.topo.phi),'k-','linewidth',2), hold on
plot(rad2deg(topo_contours.hellas.topo.theta),rad2deg(topo_contours.hellas.topo.phi),'k-','linewidth',2)
plot(rad2deg(topo_contours.argyre.topo.theta),rad2deg(topo_contours.argyre.topo.phi),'k-','linewidth',2)

colorbar('location','northoutside')
axis equal
xlim([0 180]), ylim([0 360])

title 'dimensionless.'

subplot 122
contourf(rad2deg(Xc),rad2deg(Yc),hm,50,'LineColor','none'); hold on
plot(rad2deg(topo_contours.dichotomy.topo.theta),rad2deg(topo_contours.dichotomy.topo.phi),'k-','linewidth',1), hold on
plot(rad2deg(topo_contours.hellas.topo.theta),rad2deg(topo_contours.hellas.topo.phi),'k-','linewidth',1)
plot(rad2deg(topo_contours.argyre.topo.theta),rad2deg(topo_contours.argyre.topo.phi),'k-','linewidth',1)

colorbar('location','northoutside')
axis equal
xlim([0 180]), ylim([0 360])

title 'dimensional'