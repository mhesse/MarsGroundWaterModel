function [Xc_ext,Yc_ext,F_ext] = extend_field_spherical(Xc,Yc,F,Grid)
% author: Marc Hesse
% date: Mar 30, 2022
% Description: Extends field computed on sperical shell with longitude 
% between 0 and 2*pi to -pi to 3*pi, assuming periodicity. this is helpful
% for plotting and tracing streamlines across the periodicity boundary
% Input:
% Xc, Yc = Ny by Nx coordinate matrices produced by meshgrid
% F = Ny by Nx matrix of field of interest
% Xc_ext, Yc_ext = 2*Ny by Nx extended coordinate matrices
% F_ext = 2*Ny by Nx extended field matrix
Ny = Grid.Ny;
bot = 1:Ny/2; top = Ny/2+1:Ny;

%% Divide the matrices in half
% Coordinate matrices
Xc_bot = Xc(bot,:);
Xc_top = Xc(top,:);

Yc_bot = Yc(bot,:);
Yc_top = Yc(top,:);

% Field matrix
F_bot = F(bot,:);
F_top = F(top,:);

%% Paste together extended coords and field
% The Yc coordinate (longitude) has to be shifted by 2*pi
Xc_ext = [Xc_top;Xc;Xc_bot];
Yc_ext = [Yc_top-2*pi;Yc;Yc_bot+2*pi];
F_ext  = [F_top;F;F_bot];