function [const,prop,scales] = mars_constants_properties_scales()
% date: Nov 26, 2021
% author: Marc Hesse
% Description: Define the properties and scales used for the Mars
% groundwater model.

%% Constants
% General
const.yr2s = 60^2*24*365.25; %   seconds per year
const.R = 3389508;      % [m]      Mars' mean radius 
const.gEarth = 9.81;    % [m/s^2]  Earth gravity
const.gMars = 3.711;    % [m/s^2]  Mars gravity
const.gMoon = 1.62;     % [m/s^2]  Moon gravity
const.rho = 1e3;        % [kg/m^3] desity of water 
const.mu = 1e-3;        % [Pa s]   water viscosity

% Topograhy
const.topo.high.elev =  1e3;         % [m] mean highland elevation
const.topo.low.elev  = -4e3;         % [m] mean lowland elevation
const.sea.Deuteronilus.elev = -3792; % [m] elevation of Deuteronilus shoreline
const.sea.Arabia.elev = -2090;       % [m] elevation of Arabia shoreline
const.sea.Meridiani.elev = 0;        % [m] elevation of Meridiani shoreline

% Aquifer properties
const.theta_bnd = acos(-1/3);                   % [rad] southern lattitude of dichotomy bnd
const.theta_bnd_deg = rad2deg(const.theta_bnd); % [deg] southern lattitude of dichotomy bnd
const.K_s_min = 1e-6;   % [m/s] min. hyd. conductivity on surface.
const.phi_s_min = 0.3;  % [-] min. porosity on surface
const.phi_s_max = 0.5;  % [-] max. porosity on surface
const.dmax = 10e3;      % [m] max. aquifer depth
const.d0 = 2.8375e3;    % [m] porosity decay depth
const.m_best = const.dmax/(2*const.d0*log(2)); % [-] porosity exponent fit to Clifford 1993
const.m_exp = 2.5; % [-] exponent in porosity power-law with height
const.r_exp = 3;   % [-] exponent in porosity-hyd. cond. power-law 
const.n_exp = const.m_exp*const.r_exp; % [-] exponent in hyd. cond. power-law with height

%% Properties (functions)
% reference distributions
prop.ref.clifford1993.phi = @(d,phi_s,d0) phi_s*exp(-d/d0); % [-] porosity
prop.ref.MI.k0 = 10^(-12.65); % [m^2] ref. permeability at 1km depth Manning and Ingebritsen (1999)
prop.ref.MI.alpha = 3.2;
prop.ref.MI.k = @(d) prop.ref.MI.k0_MI*d.^(-prop.ref.MI.alpha);
prop.ref.MI.K = @(d) prop.ref.MI.k_MI(d)*const.rho*const.gMars/const.mu;

% Porosity as function of depth
prop.phi0_best = @(phi_s,m) phi_s/(const.dmax^m);  % [-] porosity pre-factor fit to Clifford 1993
prop.phi = @(z,phi0,m_exp) phi0*z.^m_exp;
prop.phibar = @(h,phi_s,m) (prop.phi0_best(phi_s,m)*h.^m)/(m+1);

% Hydraulic conductivity as function of porosity
prop.C = @(m_exp,n_exp) const.K_s_min/const.phi_s_min^(n_exp/m_exp); % pre-factor in permeability-porosity relation
prop.Khat = @(phi,m_exp,n_exp) prop.C(m_exp,n_exp)*phi^(n_exp/m_exp);     % hyd. conductivity as function of porosity      
prop.K0 = @(phi_s,m_exp,n_exp,dmax) prop.C(m_exp,n_exp)*phi_s^(n_exp/m_exp)/dmax^n_exp;
prop.Kbar = @(h,phi_s,m_exp,n_exp,dmax) (prop.K0(phi_s,m_exp,n_exp,dmax)*h.^n_exp)/(n_exp+1); % [m] mean hyd. conductivity

% Hydraulic conductivity as function of depth

% Elevation to height above base of aquifer
prop.z_aq = @(elev) elev-(const.topo.high.elev-const.dmax);


%% Derived constants
% Shorelines above base of aquifer
const.sea.Deuteronilus.z = prop.z_aq(const.sea.Deuteronilus.elev); % [m] Deuteronilus shoreline
const.sea.Arabia.z = prop.z_aq(const.sea.Arabia.elev);             % [m] Arabia shoreline
const.sea.Meridiani.z = prop.z_aq(const.sea.Meridiani.elev);       % [m] Meridiani shoreline

%% Scales
scales = 1;