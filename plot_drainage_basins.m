close all, clear, clc


load ../MarsTopoProcessing/hellas_topo.mat
load ../MarsTopoProcessing/argyre_topo.mat
load ../MarsTopoProcessing/dichotomy_topo.mat
load flow_field_with_crater.mat
load drainage.mat

% drainage_area(drainage_area==4)=2;

%% Geometry
figure
surf(Xc,Yc,reshape(drainage_area,Grid.Ny,Grid.Nx))

return
h_area = plot_mars_map(Xc,Yc,reshape(drainage_area,Grid.Ny,Grid.Nx),Grid,topo_contours);

text(2.8,0,4,'Deuteronilus shoreline: -3.8 km','color','k')
text(2.65,0,4,'Arabia shoreline: -2.1 km','color','k')
text(2.5,0,4,'Meridiani shoreline: 0.0 km','color','k')
text(1.4,-.5,4,'Hellas shorelines:','color','k')
text(1.25,-.5,4,'-3.1 to -5.8 km (Wilson et al. 2010)','color','k')

text(.45,1.3,4,'Argyre shorelines:','color','k')
text(.3,1.3,4,'0 and 1 km (Fairen et al. 2016)','color','k')

% print(gcf, 'Geom.pdf','-dpdf','-r600');
print(gcf, 'Geom.png','-dpng','-r600');


%% Plot drainage areas
figure
h_area = plot_mars_map(Xc,Yc,reshape(drainage_area,Grid.Ny,Grid.Nx),Grid,topo_contours);

text(2.8,0,4,'Deuteronilus: -3.8 km','color','w')
text(.85,-.9,4,'-3.1 km','color','w')
text(.9,1.1,4,'-2.5 km','color','k')
% print(gcf, 'DrainageBasins.pdf','-dpdf','-r600');
print(gcf,'DrainageBasins.png','-dpng','-r600');

disp('done')

%% Head
figure%('position',[10 10 1.25*800 1.25*800 ])
h_head = plot_mars_map(Xc,Yc,reshape(u,Grid.Ny,Grid.Nx)*10-9,Grid,topo_contours);

text(2.5,0,4,'Arabia shoreline: -2.1 km','color','k')
text(.85,-.9,4,'-3.1 km','color','w')
text(.9,1.1,4,'-2.5 km','color','k')
% print(gcf, 'Head.pdf','-dpdf','-r600');
print(gcf,'Head.png','-dpng','-r600');

%% Recharge
figure
% surf(Xc,Yc,reshape(fs,Grid.Ny,Grid.Nx))
h_rech = plot_mars_map(Xc,Yc,reshape(fs,Grid.Ny,Grid.Nx)*10-9,Grid,topo_contours);
% colormap('gray')
text(0.2,3,'Kite, Mischna, Fan, et al., in prep.','color','w')
% print(gcf, 'Recharge.pdf','-dpdf','-r600');
print(gcf,'Recharge.png','-dpng','-r600');
