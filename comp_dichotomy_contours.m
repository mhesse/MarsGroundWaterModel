function [topo,lowlands,highlands] = comp_dichotomy_contours(elevations,theta,phi,mars_topo,dichotomy)
% author; Marc Hesse
% date: Dec 8, 2021
% description: Compute the contours of the dichotomy boundary

for k = 1:length(elevations)
    [cont,n_cont] = get_contours(theta,phi,mars_topo,elevations(k));
    if min(cont(1).y) == 0 && max(cont(1).y) == 2*pi;
        if cont(1).y(1) > cont(1).y(end)
%             fprintf('k = %d: 1 > end\n',k)
            topo(k).theta = cont(1).x;
            topo(k).phi   = cont(1).y;
            
            lowlands(k).theta = [pi;  topo(k).theta;pi;pi];
            lowlands(k).phi   = [2*pi;topo(k).phi  ;0;2*pi];
            
            highlands(k).theta = [0; 0;  topo(k).theta;0];
            highlands(k).phi   = [0;2*pi;topo(k).phi  ;0];
        else
%             fprintf('k = %d: 1 < end\n',k)
            topo(k).theta = cont(1).x;
            topo(k).phi   = cont(1).y;
            
            lowlands(k).theta = [pi;topo(k).theta;  pi;pi];
            lowlands(k).phi   = [ 0;topo(k).phi  ;2*pi; 0];
            
            highlands(k).theta = [0;topo(k).theta;   0;0];
            highlands(k).phi   = [0;topo(k).phi  ;2*pi;0];
        end
    else
        fprintf('dichotomy elevation too low')
        topo(k).theta = nan;
        topo(k).phi   = nan;
        lowlands(k).theta = nan;
        lowlands(k).phi   = nan;
    end
    dichotomy.topo(k).z     = elevations(k);
end