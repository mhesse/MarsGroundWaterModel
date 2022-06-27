function [dof,dof_f,X_f,Y_f] = get_mars_dofs(Grid,D,Xc,Yc,topo_contours)
% author: Marc Hesse
% date: Dec 8, 2021
% description: Determines the degree of fredom arrays for cells and faces 
% of the various domains on Mars that are needed in the computations.
fprintf('\nDetermining the dofs of each feature: ')
% Hellas
fprintf('Hellas, ')
[dof.hellas,dof_f.hellas,X_f.hellas,Y_f.hellas] = get_feature_dofs(Grid,D,Xc,Yc,topo_contours.hellas.topo);

% Argyre
fprintf('Argyre, ')
[dof.argyre,dof_f.argyre,X_f.argyre,Y_f.argyre] = get_feature_dofs(Grid,D,Xc,Yc,topo_contours.argyre.topo);

% Lowlands
fprintf('Lowlands, ')
[dof.lowlands,dof_f.lowlands,X_f.lowlands,Y_f.lowlands] = get_feature_dofs(Grid,D,Xc,Yc,topo_contours.lowlands);

% Highlands
fprintf('Highlands.\n\n')
[dof.highlands,dof_f.highlands,X_f.highlands,Y_f.highlands] = get_feature_dofs(Grid,D,Xc,Yc,topo_contours.highlands);

% Inactive domain
% Add Hellas and Argyre to the Lowlands
dof.inactive = sort([dof.lowlands;dof.hellas;dof.argyre]);
dof_f.inactive = find_faces(dof.inactive,D,Grid);
[X_f.inactive,Y_f.inactive] = comp_face_coords(dof_f.inactive,Grid);

% Active domain
% Subtract Hellas and Argyre from the Highlands
dof.active = setdiff(dof.highlands,[dof.hellas;dof.argyre]);
dof_f.active = find_faces(dof.active,D,Grid);
[X_f.active,Y_f.active] = comp_face_coords(dof_f.active,Grid);

[dof.active_bnd_in,dof.active_bnd_out] = find_bnd_cells(dof.active,dof.inactive,dof_f.active,D,Grid);
[dof.hellas_bnd.in,dof.hellas_bnd.out] = find_bnd_cells(dof.active,dof.inactive,dof_f.hellas,D,Grid);
[dof.argyre_bnd.in,dof.argyre_bnd.out] = find_bnd_cells(dof.active,dof.inactive,dof_f.argyre,D,Grid);

% Dirichlet boundaries
dof.dir = setdiff(dof.active_bnd_in,[Grid.dof_xmin;Grid.dof_ymin;Grid.dof_ymax]);
dof.dichotomy_bnd.in = setdiff(dof.active_bnd_in,[Grid.dof_xmin;Grid.dof_ymin;Grid.dof_ymax;dof.argyre_bnd.in;dof.hellas_bnd.in]);

save('mars_dofs.mat','dof','dof_f','X_f','Y_f');
topo_contours_old = topo_contours;
save('topo_check.mat','topo_contours_old')
end

function [dof,dof_f,X_f,Y_f] = get_feature_dofs(Grid,D,Xc,Yc,feature)
dof = Grid.dof(inpolygon(Xc(:),Yc(:),feature.theta,feature.phi));
dof_f = find_faces(dof,D,Grid);
[X_f,Y_f] = comp_face_coords(dof_f,Grid);
end
