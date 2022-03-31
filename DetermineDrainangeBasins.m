% file: DetermineDrainageBasins.m
% author: Marc Hesse
% date: Dec 8, 2021, Mar 24, 2022
close all, clear, clc

load ../MarsTopoProcessing/hellas_topo.mat
load ../MarsTopoProcessing/argyre_topo.mat
load ../MarsTopoProcessing/dichotomy_topo.mat
cd ./test2
load flow_field.mat


figure
subplot 121
% Lowlands
plot(topo_contours.dichotomy.topo.theta,topo_contours.dichotomy.topo.phi,'r-','linewidth',2), hold on
plot(target_zones.dichotomy.theta,target_zones.dichotomy.phi,'r:','linewidth',2)
plot(topo_contours.dichotomy.topo.theta,topo_contours.dichotomy.topo.phi-2*pi,'r-','linewidth',2), hold on
plot(target_zones.dichotomy.theta,target_zones.dichotomy.phi-2*pi,'r:','linewidth',2)
plot(topo_contours.dichotomy.topo.theta,topo_contours.dichotomy.topo.phi+2*pi,'r-','linewidth',2), hold on
plot(target_zones.dichotomy.theta,target_zones.dichotomy.phi+2*pi,'r:','linewidth',2)

% Hellas
plot(topo_contours.hellas.topo.theta,topo_contours.hellas.topo.phi,'b-','linewidth',2)
plot(target_zones.hellas.theta,target_zones.hellas.phi,'b:','linewidth',2)
plot(topo_contours.hellas.topo.theta,topo_contours.hellas.topo.phi-2*pi,'b-','linewidth',2)
plot(target_zones.hellas.theta,target_zones.hellas.phi-2*pi,'b:','linewidth',2)

% Argyre
plot(topo_contours.argyre.topo.theta,topo_contours.argyre.topo.phi,'g-','linewidth',2)
plot(target_zones.argyre.theta,target_zones.argyre.phi,'g:','linewidth',2)
plot(topo_contours.argyre.topo.theta,topo_contours.argyre.topo.phi+2*pi,'g-','linewidth',2)
plot(target_zones.argyre.theta,target_zones.argyre.phi+2*pi,'g:','linewidth',2)
axis equal
xlim([0 pi]), ylim([-pi 3*pi])

subplot 122
plot(topo_contours.dichotomy.topo.theta,topo_contours.dichotomy.topo.phi,'linewidth',2), hold on
plot(topo_contours.hellas.topo.theta,topo_contours.hellas.topo.phi,'linewidth',2)
plot(topo_contours.argyre.topo.theta,topo_contours.argyre.topo.phi,'linewidth',2)
axis equal
xlim([0 pi]), ylim([0 2*pi])
drawnow

%% Determine which target zone streamlines end in
drainage_area = 0*Xc(:); drainage_area(dof.inactive) = nan;
for i = 1:Ns
    x_start(i) = S{i}(1,1);
    y_start(i) = S{i}(1,2);
    x_end(i) = S{i}(end,1);
    y_end(i) = S{i}(end,2);
%     subplot 121
%     plot(S{i}(:,1),S{i}(:,2),'c-')
%     if mod(round(i/Ns*100),10) == 0; fprintf('Progress: %d percent\n',round(i/Ns*100)); end
    if inpolygon(x_end(i),y_end(i),target_zones.dichotomy.theta,target_zones.dichotomy.phi) % test Lowlands
        drainage_area(dof.active(i)) = 1;
%         subplot 121
%         plot(S{i}(end,1),S{i}(end,2),'k.')
%         subplot 122
%         plot(x_start(i),y_start(i),'k.')
    elseif inpolygon(x_end(i),y_end(i),target_zones.dichotomy.theta,target_zones.dichotomy.phi-2*pi) % test Lowlands shifted south
        drainage_area(dof.active(i)) = 1;
%         subplot 121
%         plot(S{i}(end,1),S{i}(end,2),'k.')
%         subplot 122
%         plot(x_start(i),y_start(i),'k.')
    elseif inpolygon(x_end(i),y_end(i),target_zones.dichotomy.theta,target_zones.dichotomy.phi+2*pi) % test Lowlands shifted north
        drainage_area(dof.active(i)) = 1;
%         subplot 121
%         plot(S{i}(end,1),S{i}(end,2),'k.')
%         subplot 122
%         plot(x_start(i),y_start(i),'k.')
    elseif inpolygon(x_end(i),y_end(i),target_zones.hellas.theta,target_zones.hellas.phi) % test Hellas
        drainage_area(dof.active(i)) = 2;
%         subplot 121
%         plot(S{i}(end,1),S{i}(end,2),'k.')
%         subplot 122
%         plot(x_start(i),y_start(i),'k.')
    elseif inpolygon(x_end(i),y_end(i),target_zones.hellas.theta,target_zones.hellas.phi-2*pi) % test Hellas shifted
        drainage_area(dof.active(i)) = 2;
%         subplot 121
%         plot(S{i}(end,1),S{i}(end,2),'k.')
%         subplot 122
%         plot(x_start(i),y_start(i),'k.')
    elseif inpolygon(x_end(i),y_end(i),target_zones.argyre.theta,target_zones.argyre.phi) % test Argyre
        drainage_area(dof.active(i)) = 3;
%         subplot 121
%         plot(S{i}(end,1),S{i}(end,2),'k.')
%         subplot 122
%         plot(x_start(i),y_start(i),'k.')
    elseif inpolygon(x_end(i),y_end(i),target_zones.argyre.theta,target_zones.argyre.phi+2*pi) % test Argyre
        drainage_area(dof.active(i)) = 3;
%         subplot 121
%         plot(S{i}(end,1),S{i}(end,2),'k.')
%         subplot 122
%         plot(x_start(i),y_start(i),'k.')
    elseif x_end(i) <= Grid.xmin+Grid.dx % South pole
        drainage_area(dof.active(i)) = -1;
    else
        drainage_area(dof.active(i)) = 0;
    end
     drawnow
end

%% Assign south pole to drainage basin
check_hellas = ismember(2,drainage_area(Grid.dof_xmin));
check_argyre = ismember(3,drainage_area(Grid.dof_xmin));
if check_hellas && ~check_argyre
    fprintf('South pole drains to Hellas\n')
    drainage_area(drainage_area == -1) = 2;
elseif check_argyre && ~check_hellas
    fprintf('South pole drains to Argyre\n')
    drainage_area(drainage_area == -1) = 3;
else
    % Check if the max flux at pole points into Hellas or Argure basin
    dof_hellas = Grid.dof_xmin(drainage_area(Grid.dof_xmin)==2)
    dof_argyre = Grid.dof_xmin(drainage_area(Grid.dof_xmin)==3)
    fprintf('Use qD_pole to determine basin: ')
    [qD_pole_max,ind_max] = max(qD(Grid.dof_f_xmin));
    dof_pole = Grid.dof_xmin(ind_max);
    if ismember(dof_pole,dof_hellas)
        fprintf('Hellas.\n')
        drainage_area(drainage_area == -1) = 2;
    elseif ismember(dof_pole,dof_argyre)
        fprintf('Argyre.\n')
        drainage_area(drainage_area == -1) = 3;
    else
        fprintf('?')
    end
end
save('drainage.mat','Ns','drainage_area')