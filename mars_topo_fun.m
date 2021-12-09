function [z_int] = mars_topo_fun(x_int,Theta,Phi,mars_topo)
% author: Marc Hesse
% date = Dec 8, 2021
% description: Function that interpolates the discrete Mola topography
theta_int = x_int(1);
phi_int   = x_int(2);
z_int = interp2(Theta,Phi,mars_topo,theta_int,phi_int,'linear');
