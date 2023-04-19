% This file plots density, ion velocity and ion acceleration, both over the
% whole region.


%determine how many grid points there are
num_pts = RESX * RESZ;
grid_data = csvread([path folder dataset{d} name '_ion-den.txt']);
grid_pts = grid_data(1:num_pts,:);
density = grid_data(num_pts+1:end,1);

%setting X and Z values
X = reshape(grid_pts(:,1),RESX,RESZ); 
Z = reshape(grid_pts(:,2),RESX,RESZ); 
center_Z = Z(1,:); %Z is changing across the columns

%finding density
den = reshape(density,RESX,RESZ,round(length(density)/num_pts));
final_den = mean(den(:,:,round(size(den,3)/2):end),3); %average over last 1/2
avg_denz = mean(final_den);  %1x128

figure(10)
plot(center_Z*1e3,(avg_denz)/DEN_FAR_PLASMA, ...
       'Marker','.',...
       'MarkerEdgeColor',purple1,...
       'MarkerFaceColor',purple1,...
        'MarkerSize',7, 'LineWidth',2,'Color',purple1)
    
xlabel('z (mm)', 'FontWeight', 'bold', 'FontSize', 20);
ylabel('n_i/n_{io}', 'FontWeight', 'bold', 'FontSize', 20);
set(findobj(gcf,'type','axes'),'FontSize',15);

%% ion velocity and acceleration (position of z) 
bin_edgesz = linspace(-HT_CYL, HT_CYL, RESZ+1);
dz = (bin_edgesz(2)-bin_edgesz(1))/DEBYE; %in units of debye length
bin_centersz = diff(bin_edgesz)/2 + bin_edgesz(1:end-1);

data = csvread([path folder dataset{d} name '_trace.txt']);
ionPosz = data(1:3:end,3);
%ionPos = data(1:3:end,:);
ionVelz = data(2:3:end,3);
% ionAccx = data(3:3:end,1);
% ionAccy = data(3:3:end,2);
ionAccz = data(3:3:end,3);
% ionVel = csvread([path folder name '_ion-vel.txt']);
% ionVelz = ionVel(:,2);
% ionPos = csvread([path folder name '_ion-pos.txt']);
% ionPosz = ionPos(:,3); 

ion_vel_z = zeros(1, RESZ);
ion_acc_z = zeros(1, RESX);

for binz = 1: length(bin_centersz)
    if binz == 1
        qz = find(ionPosz > bin_edgesz(binz) & ionPosz <= bin_edgesz(binz+1));
    else
        qz = find(ionPosz > bin_edgesz(binz) & ionPosz <=bin_edgesz(binz+1));
    end
    ion_vel_z(binz) = (nanmean(ionVelz(qz)));
    ion_acc_z(binz) = (nanmean(ionAccz(qz)));
end
    
figure(11)
plot((center_Z)*1e3,ion_vel_z/SOUND_SPEED,  ...
        'Marker','none','LineStyle','-',...
        'MarkerEdgeColor',purple1,...
         'MarkerFaceColor',purple1,...
         'MarkerSize',7, 'LineWidth',1.5,'Color',purple1);
     
xlabel('z (mm)', 'FontWeight', 'bold', 'FontSize', 20);
ylabel('v_{iz}/C_s', 'FontWeight', 'bold', 'FontSize', 20);
set(findobj(gcf,'type','axes'),'FontSize',15);

figure(12)
plot(center_Z*1e3,ion_acc_z, ...
        'Marker','none','LineStyle','-',...
        'MarkerEdgeColor',purple1,...
        'MarkerFaceColor',purple1,...
        'MarkerSize',7, 'LineWidth',1.5,'Color',purple1)
    
xlabel('z (mm)', 'FontWeight', 'bold', 'FontSize', 20);
ylabel('a_{iz}', 'FontWeight', 'bold', 'FontSize', 20);
set(findobj(gcf,'type','axes'),'FontSize',15);

%% In region of dust
% Density
% if(E_ALPHA > 0)
%     figure
%     plot(center_Z*1e3,(avg_denz), ...
%             'Marker','none','LineStyle','-',...
%             'MarkerEdgeColor',blue,...
%             'MarkerFaceColor',blue,...
%             'MarkerSize',7, 'LineWidth',2,'Color',blue)
%     xline(Z_MAX*1e3)
%     xline(Z_MIN*1e3)
% end

% to set axis
% xlimits = find(center_Z*1e3 > Z_MIN*1e3 & center_Z*1e3<Z_MAX*1e3);
% xlim([center_Z(xlimits(1))*1e3 center_Z(xlimits(end))*1e3])
% hold on
% xline(Z_MIN*1e3)
% xline(Z_MAX*1e3)
% ylimits = find(center_Z*1e3 == Z_MIN*1e3 & center_Z*1e3 == Z_MAX*1e3);
% ylim([center_Z(find(ylimits))*1e3 center_Z(find(ylimits))])

% for shaded region corresponding to area spanned by dust chain
% hold on
% minden = min(avg_denz/DEN_FAR_PLASMA);
% maxden = max(avg_denz/DEN_FAR_PLASMA);
% x2 = [-1.5 1.5 1.5 -1.5];
% y2 = [minden minden maxden maxden];
% patch(x2,y2,blue)
% alpha(0.2)

% ylabel('n_{i}/n_{io}','Fontsize',15)
% xlabel('z (mm)', 'Fontsize', 15)
% title('n_i/n_{io}', 'Fontsize', 15);
% saveas(gcf, [path folder 'zoomed_ni_nio.jpg']);

% Velocity
% figure
% plot(center_Z*1e3,ion_vel_z/SOUND_SPEED, ...
%         'Marker','none','LineStyle','-',...
%         'MarkerEdgeColor',blue,...
%         'MarkerFaceColor',blue,...
%         'MarkerSize',7, 'LineWidth',1.5,'Color',blue)
% to set axis
% limits = find(center_Z*1e3 > -4 & center_Z*1e3<4);
% xlim([center_Z(limits(1))*1e3 center_Z(limits(end))*1e3])
% ylabel('v_{iz}/C_s','Fontsize',15)
% xlabel('z (mm)', 'Fontsize', 15)

% for shaded region corresponding to area spanned by dust chain
% hold on
% x2 = [4.67 7.33 7.33 4.67];
% y2 = [0.0 0.0 1.6 1.6];
% patch(x2,y2,blue)
% title(dataset{d}, 'Interpreter', 'none', 'Units', 'normalized', 'Position', [.5, 1, 1])
%     saveas(gcf, [path folder 'vi.jpg']);

% Acceleration
% figure
% plot(center_Z*1e3,ion_acc_x, ...
%         'Marker','none','LineStyle','-',...
%         'MarkerEdgeColor',red,...
%         'MarkerFaceColor',red,...
%         'MarkerSize',7, 'LineWidth',1.5,'Color',red)
% to set axis
% limits = find(center_Z*1e3 > -4 & center_Z*1e3<4);
% xlim([center_Z(limits(1))*1e3 center_Z(limits(end))*1e3])
% ylabel('a_i','Fontsize',15)
% xlabel('z (mm)', 'Fontsize', 15)
% for shaded region corresponding to area spanned by dust chain
% hold on
% x2 = [4.67 7.33 7.33 4.67];
% y2 = [0.0 0.0 1.6 1.6];
% patch(x2,y2,blue)
% title(dataset{d}, 'Interpreter', 'none', 'Units', 'normalized', 'Position', [.5, 1, 1])