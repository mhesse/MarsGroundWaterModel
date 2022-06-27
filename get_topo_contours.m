function [topo_contours] = get_topo_contours(topo_contours,theta,phi,mars_topo,dichotomy,hellas,argyre) % repo
% date: Dec 9, 2021
% author: Marc Hesse
% description: Extracts the topographis contours of dichotomy, Hellas and
% Argyre given the elevations specified for each one in topo_contours.hellas = comp_crater_contours(elevations_hellas,theta,phi,mars_topo,hellas);

%% Dichotomy contour
[topo_contours.dichotomy.topo,topo_contours.lowlands,topo_contours.highlands] = comp_dichotomy_contours(topo_contours.dichotomy.z,theta,phi,mars_topo,dichotomy);
topo_contours.lowlands.z = topo_contours.dichotomy.z;
topo_contours.highlands.z = topo_contours.dichotomy.z;

%% Hellas contour
topo_contours.hellas.topo = comp_crater_contours(topo_contours.hellas.z,theta,phi,mars_topo,hellas);
% topo_contours.hellas = comp_crater_contours(topo_contours.hellas.z,theta,phi,mars_topo,hellas);


%% Argyre contour
topo_contours.argyre.topo = comp_crater_contours(topo_contours.argyre.z,theta,phi,mars_topo,argyre);
% topo_contours.argyre = comp_crater_contours(topo_contours.argyre.z,theta,phi,mars_topo,argyre);


