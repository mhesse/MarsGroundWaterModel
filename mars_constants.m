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
const.gen.A_surf = 4*pi*const.gen.R^2; % [m2] Mars surface area

% Topograhy
const.topo.high.elev =  1e3;         % [m] mean highland elevation
const.topo.low.elev  = -4e3;         % [m] mean lowland elevation
const.sea.Deuteronilus.elev = -3792; % [m] elevation of Deuteronilus shoreline
const.sea.Arabia.elev = -2090;       % [m] elevation of Arabia shoreline
const.sea.Meridiani.elev = 0;        % [m] elevation of Meridiani shoreline

% Elevation to height above base of aquifer
z_aq = @(elev) elev-(const.topo.high.elev-const.aq.dmax);

% Convert shorelines to height above aquifer
const.sea.Deuteronilus.z = z_aq(const.sea.Deuteronilus.elev); % [m] elevation of Deuteronilus shoreline
const.sea.Arabia.z       = z_aq(const.sea.Arabia.elev);       % [m] elevation of Arabia shoreline
const.sea.Meridiani.z    = z_aq(const.sea.Meridiani.elev);    % [m] elevation of Meridiani shoreline

% Highlands surface area
const.sea.Deuteronilus.highland_area = 1.204814697468799e+14;
const.sea.Arabia.highland_area       = 9.646692944957303e+13;
const.sea.Meridiani.highland_area    = 7.663773386114716e+13;

% Highlands surface area fraction
const.sea.Deuteronilus.highland_area_frac = 0.834520608857827;
const.sea.Arabia.highland_area_frac       = 0.668182757630977;
const.sea.Meridiani.highland_area_frac    = 0.530834894840303;

% Lowlands surface area
const.sea.Deuteronilus.lowland_area = 2.389060263582352e+13;
const.sea.Arabia.lowland_area       = 4.790514293312916e+13;
const.sea.Meridiani.lowland_area    = 6.773433852155270e+13;

% Highlands surface area fraction
const.sea.Deuteronilus.lowland_area_frac = 1-const.sea.Deuteronilus.highland_area_frac;
const.sea.Arabia.lowland_area_frac       = 1-const.sea.Arabia.highland_area_frac;
const.sea.Meridiani.lowland_area_frac    = 1-const.sea.Meridiani.highland_area_frac;

% Southern co-lattitude Dichotomy boundary
const.sea.Deuteronilus.theta_bnd     = acos(1-2*const.sea.Deuteronilus.highland_area_frac);
const.sea.Arabia.theta_bnd           = acos(1-2*const.sea.Arabia.highland_area_frac);
const.sea.Meridiani.theta_bnd        = acos(1-2*const.sea.Meridiani.highland_area_frac);

const.sea.Deuteronilus.theta_bnd_deg = rad2deg(const.sea.Deuteronilus.theta_bnd);
const.sea.Arabia.theta_bnd_deg    = rad2deg(const.sea.Arabia.theta_bnd);
const.sea.Meridiani.theta_bnd_deg       = rad2deg(const.sea.Meridiani.theta_bnd);


% Aquifer properties
% const.aq.theta_bnd = acos(-1/3);                   % [rad] southern lattitude of dichotomy bnd
% const.aq.theta_bnd_deg = rad2deg(const.aq.theta_bnd); % [deg] southern lattitude of dichotomy bnd
const.aq.K_s_min = 1e-6;   % [m/s] min. hyd. conductivity on surface.
const.aq.phi_s_min = 0.3;  % [-] min. porosity on surface
const.aq.phi_s_max = 0.5;  % [-] max. porosity on surface
const.aq.d0 = 2.8375e3;    % [m] porosity decay depth
const.aq.m_best = const.aq.dmax/(2*const.aq.d0*log(2)); % [-] porosity exponent fit to Clifford 1993
const.aq.n_exp = const.aq.m_exp*const.aq.r_exp; % [-] exponent in hyd. cond. power-law with height