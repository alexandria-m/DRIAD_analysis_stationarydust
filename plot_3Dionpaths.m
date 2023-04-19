% plot ions to show insertion locations
clear all
path = '/Users/alexandriamendoza/Documents/Plasma/changing_density_new/';
folder = 'debug_june1/gaus_quad/';
name = 'linearconstant';

import_without_A;
fileID = fopen([path folder name '_dust-final.txt']);
data = textscan(fileID,'%s%s%s%s%s%s%s%s','Delimiter',[' ','\t'],'MultipleDelimsAsOne',1);
posx = data{1,2};
posy = data{1,3};
posz = data{1,4};
fclose(fileID);
   
positions(:,1) = str2double(posx(2:end));
positions(:,2) = str2double(posy(2:end));
positions(:,3) = str2double(posz(2:end));
command_list = zeros(1,5);
command_list(1) = 1;
command_list(2) = 1;
command_list(3) = 1;

fileID = fopen([path folder name '_trace.txt']);
if(fileID ~= -1)
    if(command_list(2) == 1 && command_list(3) == 1)
        data = textscan(fileID,'%f%f%f%f','Delimiter',',');
        ion_x = data{1,1};
        ion_y = data{1,2};
        ion_z = data{1,3};
        ion_q = data{1,4};
        ion_pos(:,1) = ion_x(1:3:end,1);
        ion_pos(:,2) = ion_y(1:3:end,1);
        ion_pos(:,3) = ion_z(1:3:end,1);
        ion_vel(:,1) = ion_x(2:3:end,1);
        ion_vel(:,2) = ion_y(2:3:end,1);
        ion_vel(:,3) = ion_z(2:3:end,1);
        ion_acc(:,1) = ion_x(3:3:end,1);
        ion_acc(:,2) = ion_y(3:3:end,1);
        ion_acc(:,3) = ion_z(3:3:end,1);
        ion_q = ion_q(1:3:end,1);
    elseif(command_list(2) == 1)
        data = textscan(fileID,'%f%f%f%f','Delimiter',',');
        ion_pos(:,1) = data{1,1};
        ion_pos(:,2) = data{1,2};
        ion_pos(:,3) = data{1,3};
        ion_q = data{1,4};
        ion_vel(:,1) = data{1,5};
        ion_vel(:,2) = data{1,6};
        ion_vel(:,3) = data{1,7};
    elseif(command_list(3) == 1)
        data = textscan(fileID,'%f%f%f%f','Delimiter',',');
        ion_pos(:,1) = data{1,1};
        ion_pos(:,2) = data{1,2};
        ion_pos(:,3) = data{1,3};
        ion_q = data{1,4};
        ion_acc(:,1) = data{1,5};
        ion_acc(:,2) = data{1,6};
        ion_acc(:,3) = data{1,7};
    else
        data = textscan(fileID,'%f%f%f%f','Delimiter',',');
        ion_pos(:,1) = data{1,1};
        ion_pos(:,2) = data{1,2};
        ion_pos(:,3) = data{1,3};
        ion_q = data{1,4};
    end
end
fclose(fileID);
sts = find(diff(ion_pos(:,3)) > DEBYE_I);
starts = ion_pos(sts+1,:);
ends = ion_pos(sts,:);


figure
fig1 = gcf;
set(gcf,'Position',[350  150  450 1050],'color','w')
scatter3(positions(:,1)*1e3,positions(:,2)*1e3,positions(:,3)*1e3,80,'o',...
    'MarkerFaceColor',[0.45 0.45 0.45],'MarkerEdgeColor','k','Linewidth',0.5)
hold on
scatter3(ion_pos(:,1)*1e3,ion_pos(:,2)*1e3,ion_pos(:,3)*1e3,5,'.',...
    'MarkerEdgeColor',[187 204 238]/255)
scatter3(starts(:,1)*1e3,starts(:,2)*1e3,starts(:,3)*1e3,'o','filled',...
    'MarkerFaceColor',[0 153 136]/255)
scatter3(ends(:,1)*1e3,ends(:,2)*1e3,ends(:,3)*1e3,'o','filled',...
    'MarkerFaceColor',[238 51 119]/255)
axis equal
ax = fig1.CurrentAxes;
ax.FontSize = 18;

clear fig1 ax