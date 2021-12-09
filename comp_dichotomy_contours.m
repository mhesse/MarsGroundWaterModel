function [dichotomy] = comp_dichotomy_contours(elevations,theta,phi,mars_topo,dichotomy)
% author; Marc Hesse
% date: Dec 8, 2021
% description: Compute the contours of the dichotomy boundary

for k = 1:length(elevations)
    [cont,n_cont] = get_contours(theta,phi,mars_topo,elevations(k));
    if min(cont(1).y) == 0 && max(cont(1).y) == 2*pi;
        if cont(1).y(1) > cont(1).y(end)
%             fprintf('k = %d: 1 > end\n',k)
            dichotomy.topo(k).theta = cont(1).x;
            dichotomy.topo(k).phi   = cont(1).y;
            dichotomy.lowlands(k).theta = [pi;  dichotomy.topo(k).theta;pi;pi];
            dichotomy.lowlands(k).phi   = [2*pi;dichotomy.topo(k).phi  ;0;2*pi];
        else
%             fprintf('k = %d: 1 < end\n',k)
            dichotomy.topo(k).theta = cont(1).x;
            dichotomy.topo(k).phi   = cont(1).y;
            dichotomy.lowlands(k).theta = [pi;dichotomy.topo(k).theta;  pi;pi];
            dichotomy.lowlands(k).phi   = [ 0;dichotomy.topo(k).phi  ;2*pi; 0];
        end
    else
        fprintf('dichotomy elevation too low')
        dichotomy.topo(k).theta = nan;
        dichotomy.topo(k).phi   = nan;
        dichotomy.lowlands(k).theta = nan;
        dichotomy.lowlands(k).phi   = nan;
    end
    dichotomy.topo(k).z     = elevations(k);
end