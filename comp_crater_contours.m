function [crater] = comp_crater_contours(elevations,theta,phi,mars_topo,crater)
% file: comp_crater_contours.m
% author: Marc Hesse
% date: Dec 9, 2021
% description: Select the crater contour by checking which contour contains
% the deepest part of the crater
if min(elevations) < crater.z_min
    error('elevations are below crater floor: %3.2f.\n',crater.z_min)
end
for k = 1:length(elevations)
    [cont,n_cont] = get_contours(theta,phi,mars_topo,elevations(k));
    for i=1:n_cont
        if inpolygon(crater.theta_min,crater.phi_min,cont(i).x,cont(i).y)
            crater.topo(k).theta = cont(i).x;
            crater.topo(k).phi   = cont(i).y;
            crater.topo(k).z     = elevations(k);
            i = n_cont;
        end
    end
end