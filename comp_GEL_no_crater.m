% function: comp_GEL_no_crater.m
% author: Marc Hesse
% date: Jun 28
close all, clear, clc
load flow_field_no_crater.mat
load ../MarsTopoProcessing/Mars_topography.mat


%% Interpolate topography to grid
mars_topo_grid = interp2(Theta,Phi,mars_topo,Xc,Yc);
mars_topo_grid_fill = mars_topo_grid; % Initialize the map to be filled in

% Fill Hellas
load ../MarsTopoProcessing/hellas_topo.mat
dof_hellas_fill = Grid.dof(inpolygon(Xc(:),Yc(:),hellas.set.topo(93).theta,hellas.set.topo(93).phi)); 
mars_topo_grid_fill(dof_hellas_fill) = hellas.set.topo(93).z;

% fill Argyre
load ../MarsTopoProcessing/argyre_topo.mat
dof_argyre_fill = Grid.dof(inpolygon(Xc(:),Yc(:),argyre.set.topo(44).theta,argyre.set.topo(44).phi)); 
mars_topo_grid_fill(dof_argyre_fill) = argyre.set.topo(44).z;

save('fill.mat','dof_hellas_fill','dof_argyre_fill','mars_topo_grid_fill')

%% Compute depth of the water table
gw_depth = mars_topo_grid_fill/1e3 - zm;
gw_depth(dof.inactive) = nan;

%% Find seepage areas and set GWT to topo in seepage areas
dof_seepage = Grid.dof(gw_depth<=0);
% dof_seepage_highlands = intersect(dof_seepage,dof.active);
dof_seepage_highlands = dof_seepage;

if length(dof_seepage_highlands) > 0
    % replace the gwt elevation with topography in seepage areas
    zm_seepage = zm;
    zm_seepage(dof_seepage_highlands) = mars_topo_grid_fill(dof_seepage_highlands)/1e3;
end
% Convert back to head (above base of aquifer) 
head = zm_seepage*1e3 - const.aq.z_bot; % [m] thickness of the water layer

%% Compute GW volume
phi = 0.01;
Grid.R_shell = const.gen.R;
Grid = build_grid(Grid);
V_gw_cell = zeros(Grid.N,1); % Initialize
V_gw_cell = phi*head(:).*Grid.A; % [m^3] volume in each grid cell 
V_gw_tot = sum(V_gw_cell);

A_mars = 4*pi*const.gen.R^2; % [m^2] Mars total surface area

fprintf('Total groundwater volume: %3.2e km3.\n',V_gw_tot/1e9)
fprintf('GEL due to groundwater: %3.2f m.\n',V_gw_tot/A_mars)

%% Compute northern ocean volume
ocean_depth = zeros(Grid.N,1); % Initialize
ocean_depth(dof.lowlands) = topo_contours.dichotomy.z - mars_topo_grid(dof.lowlands);
ocean_depth(ocean_depth<0) = 0;

V_oc_cell = zeros(Grid.N,1); % Initialize
V_oc_cell(dof.lowlands) = ocean_depth(dof.lowlands).*Grid.A(dof.lowlands); % [m^3] volume in each grid cell 
V_oc_tot = sum(V_oc_cell);
V_oc_cell(dof.highlands) = nan;

fprintf('Total ocean volume: %3.2e km3.\n',V_oc_tot/1e9)
fprintf('GEL due to ocean: %3.2f m.\n',V_oc_tot/A_mars)

%% Plot all the different pieces
figure
subplot 181
contourf(rad2deg(Xc),rad2deg(Yc),zm,50,'LineColor','none'); hold on
plot(rad2deg(topo_contours.dichotomy.topo.theta),rad2deg(topo_contours.dichotomy.topo.phi),'k-','linewidth',2), hold on
colorbar('location','northoutside')
axis equal
xlim([0 180]), ylim([0 360])
title 'GWT elevation'

subplot 182
contourf(rad2deg(Xc),rad2deg(Yc),mars_topo_grid_fill/1e3,50,'LineColor','none'); hold on
plot(rad2deg(topo_contours.dichotomy.topo.theta),rad2deg(topo_contours.dichotomy.topo.phi),'k-','linewidth',1), hold on
colorbar('location','northoutside')
axis equal
xlim([0 180]), ylim([0 360])
title 'Topography'

subplot 183
contourf(rad2deg(Xc),rad2deg(Yc),gw_depth,50,'LineColor','none'); hold on
% contourf(rad2deg(Xc),rad2deg(Yc),gw_depth,[0 0],'w')
plot(rad2deg(topo_contours.dichotomy.topo.theta),rad2deg(topo_contours.dichotomy.topo.phi),'k-','linewidth',1), hold on
colorbar('location','northoutside')
axis equal
xlim([0 180]), ylim([0 360]), caxis([-5 5])
title 'Depth of GWT'

subplot 184
contour(rad2deg(Xc),rad2deg(Yc),gw_depth,[0 0],'r-','linewidth',2), hold on
plot(rad2deg(topo_contours.dichotomy.topo.theta),rad2deg(topo_contours.dichotomy.topo.phi),'k-','linewidth',2)
plot(rad2deg(Xc(dof_seepage_highlands)),rad2deg(Yc(dof_seepage_highlands)),'.')
axis equal
xlim([0 180]), ylim([0 360])
title 'areas with seepage in highlands'

subplot 185
contourf(rad2deg(Xc),rad2deg(Yc),zm_seepage,50,'LineColor','none'); hold on
plot(rad2deg(topo_contours.dichotomy.topo.theta),rad2deg(topo_contours.dichotomy.topo.phi),'k-','linewidth',2), hold on
colorbar('location','northoutside')
axis equal
xlim([0 180]), ylim([0 360])
caxis([-3 2])
title 'GWT elevation with seepage'

subplot 186
contourf(rad2deg(Xc),rad2deg(Yc),head/1e3,50,'LineColor','none'); hold on
plot(rad2deg(topo_contours.dichotomy.topo.theta),rad2deg(topo_contours.dichotomy.topo.phi),'k-','linewidth',2), hold on
colorbar('location','northoutside')
axis equal
xlim([0 180]), ylim([0 360])
% caxis([-3 2])
title 'GW thickness'

subplot 187
contourf(rad2deg(Xc),rad2deg(Yc),reshape(V_gw_cell,Grid.Ny,Grid.Nx)/1e9,50,'LineColor','none'); hold on
plot(rad2deg(topo_contours.dichotomy.topo.theta),rad2deg(topo_contours.dichotomy.topo.phi),'k-','linewidth',2), hold on
colorbar('location','northoutside')
axis equal
xlim([0 180]), ylim([0 360])
% caxis([-3 2])
title 'GW volume per cell [km^3]'

subplot 188
contourf(rad2deg(Xc),rad2deg(Yc),reshape(V_oc_cell,Grid.Ny,Grid.Nx)/1e9,50,'LineColor','none'); hold on
plot(rad2deg(topo_contours.dichotomy.topo.theta),rad2deg(topo_contours.dichotomy.topo.phi),'k-','linewidth',2), hold on
colorbar('location','northoutside')
axis equal
xlim([0 180]), ylim([0 360])
% caxis([-3 2])
title 'Ocean volume per cell [km^3]'