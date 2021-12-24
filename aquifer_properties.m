function [aq_prop,const] = aquifer_properties(const,sim_type)
% author: Marc Hesse
% date: Dec 23, 2021
% Description: Define the aquifer model

%% Check for type of simulation
if ~isfield(sim_type,'heterogeneity_lateral');  sim_type.heterogeneity_lateral  = 'no'; end % no lateral heterogeneity
if ~isfield(sim_type,'heterogeneity_vertical'); sim_type.heterogeneity_vertical = 'no'; end % no vertical heterogeneity

%% Properties (functions)
if strcmp(sim_type.heterogeneity_vertical,'power_law')
    % reference distributions
    aq_prop.ref.clifford1993.phi = @(d,phi_s,d0) phi_s*exp(-d/d0); % [-] porosity
    aq_prop.ref.MI.k0 = 10^(-12.65); % [m^2] ref. permeability at 1km depth Manning and Ingebritsen (1999)
    aq_prop.ref.MI.alpha = 3.2;
    aq_prop.ref.MI.k = @(d) aq_prop.ref.MI.k0*d.^(-aq_prop.ref.MI.alpha);
    aq_prop.ref.MI.K = @(d) aq_prop.ref.MI.k(d)*const.gen.rho*const.gen.gMars/const.gen.mu;
    
    % Porosity as function of depth
    aq_prop.phi0_best = @(phi_s,m) phi_s/(const.aq.dmax^m);  % [-] porosity pre-factor fit to Clifford 1993
    aq_prop.phi       = @(z,phi0,m_exp) phi0*z.^m_exp;
    aq_prop.phi_plot  = @(z,phi_s,m_exp) aq_prop.phi0_best(phi_s,m_exp)*z.^m_exp;
    aq_prop.phibar    = @(h,phi_s,m) (aq_prop.phi0_best(phi_s,m)*h.^m)/(m+1);
    
    % Hydraulic conductivity as function of porosity
    aq_prop.C = @(r_exp) const.aq.K_ref/const.aq.phi_ref^(r_exp); % pre-factor in permeability-porosity relation
    aq_prop.Khat = @(phi,r_exp) aq_prop.C(r_exp)*phi^(r_exp);     % hyd. conductivity as function of porosity
    aq_prop.K0 = @(phi_s,r_exp,n_exp,dmax) aq_prop.C(r_exp)*phi_s^(r_exp)/dmax^n_exp;
    aq_prop.K  = @(z,phi_s,r_exp,n_exp,dmax)  aq_prop.K0(phi_s,r_exp,n_exp,dmax)*z.^n_exp;
    aq_prop.Kbar = @(h,phi_s,r_exp,n_exp,dmax) (aq_prop.K0(phi_s,r_exp,n_exp,dmax)*h.^n_exp)/(n_exp+1); % [m] mean hyd. conductivity
    
    % Hydraulic conductivity as function of depth
    
else
    error('Undefined vertical property variation')
end

% Elevation to height above base of aquifer
aq_prop.z_aq = @(elev) elev-(const.topo.high.elev-const.aq.dmax);

% Convert shorelines to height above aquifer
const.sea.Deuteronilus.z = aq_prop.z_aq(const.sea.Deuteronilus.elev); % [m] elevation of Deuteronilus shoreline
const.sea.Arabia.z       = aq_prop.z_aq(const.sea.Arabia.elev);       % [m] elevation of Arabia shoreline
const.sea.Meridiani.z    = aq_prop.z_aq(const.sea.Meridiani.elev);    % [m] elevation of Meridiani shoreline
