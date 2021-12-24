function [const] = mars_constants(const) % Repo
% date: Nov 26, 2021
% author: Marc Hesse
% Description: Define constants used for the Mars groundwater model.

% General
const.gen.yr2s = 60^2*24*365.25; %   seconds per year
const.gen.dy2s = 60^2*24; %   seconds per day
const.gen.R = 3389508;      % [m]      Mars' mean radius 
const.gen.gEarth = 9.81;    % [m/s^2]  Earth gravity
const.gen.gMars = 3.711;    % [m/s^2]  Mars gravity
const.gen.gMoon = 1.62;     % [m/s^2]  Moon gravity
const.gen.rho = 1e3;        % [kg/m^3] desity of water 
const.gen.mu = 1e-3;        % [Pa s]   water viscosity

% Topograhy
const.topo.high.elev =  1e3;         % [m] mean highland elevation
const.topo.low.elev  = -4e3;         % [m] mean lowland elevation
const.sea.Deuteronilus.elev = -3792; % [m] elevation of Deuteronilus shoreline
const.sea.Arabia.elev = -2090;       % [m] elevation of Arabia shoreline
const.sea.Meridiani.elev = 0;        % [m] elevation of Meridiani shoreline

% Aquifer properties
const.aq.theta_bnd = acos(-1/3);                   % [rad] southern lattitude of dichotomy bnd
const.aq.theta_bnd_deg = rad2deg(const.aq.theta_bnd); % [deg] southern lattitude of dichotomy bnd
const.aq.K_s_min = 1e-6;   % [m/s] min. hyd. conductivity on surface.
const.aq.phi_s_min = 0.3;  % [-] min. porosity on surface
const.aq.phi_s_max = 0.5;  % [-] max. porosity on surface
const.aq.d0 = 2.8375e3;    % [m] porosity decay depth
const.aq.m_best = const.aq.dmax/(2*const.aq.d0*log(2)); % [-] porosity exponent fit to Clifford 1993
const.aq.n_exp = const.aq.m_exp*const.aq.r_exp; % [-] exponent in hyd. cond. power-law with height