function [scales]=comp_mars_scales(scales,aq_prop,const)
% author: Marc Hesse
% date: Dec 23, 2021

%% Compute the dimensionless scales
% Char. aquifer properties 
scales.phi_c = aq_prop.phibar(scales.h_c,const.aq.phi_s,const.aq.m_exp);
scales.K_c = aq_prop.Kbar(scales.h_c,const.aq.phi_s,const.aq.m_exp,const.aq.n_exp,const.aq.dmax);

% char. recharge flux and rate
scales.r_c = const.re.r_mean;
scales.Qr_c = scales.r_c*scales.x_c^2;

% char. head change, flux, and flow rate due to recharge
scales.dh = sqrt(scales.Qr_c/scales.K_c)
scales.q_c = scales.K_c*scales.dh/scales.x_c;
scales.Q_c = scales.x_c*scales.h_c*scales.q_c;
scales.dt = scales.phi_c*scales.x_c^2/(scales.K_c*scales.dh);
scales.dt_yr = scales.dt/const.gen.yr2s;

%% Compute dimensionless gov. paramters
scales.gov.Ss  = scales.dt/scales.t_c*scales.dh/scales.h_c;
scales.gov.Ss2 = scales.phi_c*scales.x_c^2/(scales.K_c*scales.h_c*scales.t_c);
scales.gov.Re  = scales.Qr_c/scales.Q_c*scales.dh/scales.h_c;
scales.gov.Re2 = scales.r_c*scales.x_c^2/(scales.K_c*scales.h_c^2);