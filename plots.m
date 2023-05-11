%% This file calls in all other files that plot potential, density, and ion velocity/acceleration.
clear all

%% Set path to data
%The path to your data
path = '/Users/alexandriamendoza/Documents/Plasma/zigzag';
%Folder your data is in
folder = '/debugging';
%Name of each dataset, can read in multiple at a time
dataset = {'/103_change_linear_nodust'};
%Name of run (must be the same for each dataset)
name = '/103test';

%% Different colors for plots
purple1 = [170 60 195]/255;
purple2 = [170 120 245]/255;
blue1 = [140 200 230]/255;
blue2 = [100 150 195]/255;
red1 = [255 100 153]/255;
red2 = [204 20 102]/255;
green1 = [70 200 90]/255;
green2 = [50 250 0]/255;
orange1 = [255 128 0]/255;
orange2 = [255 200 0]/255;


%% Loop over all datasets
for d = 1:length(dataset)
    close all
    %% Read in all parameters
    fid = fopen([path folder dataset{d} name '_params.txt']);
    while ~feof(fid)
        tline = fgetl(fid);
        [~,q] = find(tline == '%');
        tempval = deblank(tline(1:q-1));
        tempvar = tline(q+2:end);
        eval([tempvar ' = ' tempval ';']);
    end
    fclose(fid);
    
    %Read in data from debug file
    import_debug_file;
    %% Set-up grid
    grid_data = csvread([path folder dataset{d} name '_ion-den.txt']);
    % r_pos = sqrt(sum(dustPos(:,1:2).^2,2));
    num_pts = RESX * RESZ;
    grid_pts = grid_data(1:num_pts,:);
    potential = grid_data(num_pts+1:end,2);
    pot = reshape(potential,RESX, RESZ,round(length(potential)/num_pts));

    X = reshape(grid_pts(:,1),RESX, RESZ);
    Z = reshape(grid_pts(:,2),RESX, RESZ);
    %% Plots
    %Dust specific cases
    if(NUM_DUST > 0)
      dustPos = csvread([path folder dataset{d} name '_dust-pos.txt']);
      dustPosz = dustPos(:,3);
      %Plots final dust charge
      dust_charge_avg;
      %Plots dust charge at each time step
      charge_time;
      %Plots Fid (in x, y and z) and Fdd (in x, y and z) for each dust grain
      forces;
    end
    
    % Plot potential map
    plot_potentialscontour;
    % Plot ion density, velocity, and acceleration
    plot_den_vel_acc;
    % Plot density map
    plot_densitycontour;

    %% Save Plots as jpgs
    if(NUM_DUST > 0)
        saveas(figure(1), [path folder dataset{d}  '/charge.jpg']);
        saveas(figure(2), [path folder dataset{d} '/charge_time.jpg']);
        saveas(figure(3), [path folder dataset{d} '/fdd.jpg']);
        saveas(figure(4), [path folder dataset{d} '/fdirect.jpg']);
        saveas(figure(5), [path folder dataset{d} '/forbit.jpg'])
        saveas(figure(9), [path folder dataset{d} '/zoomed_potential.jpg']);
        saveas(figure(14), [path folder dataset{d} '/zoomed_density.jpg']);
    end
    saveas(figure(6), [path folder dataset{d} '/inside_potential.jpg']);
    saveas(figure(7), [path folder dataset{d} '/outside_potential.jpg']);
    saveas(figure(8), [path folder dataset{d} '/insideoutside_potential.jpg']);
    saveas(figure(10), [path folder dataset{d} '/ni_nio.jpg']);
    saveas(figure(11), [path folder dataset{d} '/vi.jpg']);
    saveas(figure(12), [path folder dataset{d} '/acc_iz.jpg']);
    saveas(figure(13), [path folder dataset{d} '/normalized_density.jpg']);

%% Positions and Velocities of ions
data1 = csvread([path folder dataset{d} name '_trace.txt']);
ionPos = data1(1:3:end,:);
vel = data(2:3:end,1:3);
time = 0:ION_TIME_STEP:N_IONDT_PER_DUSTDT*NUM_TIME_STEP*ION_TIME_STEP;
time(1) = [];
% 
figure
insert=diff(ionPos(:,3));
idx = find(diff(ion_pos(:,3))>0);
idx = find(insert>0);
scatter3(ionPos(:,1),ionPos(:,2),ionPos(:,3))
hold on
scatter3(ionPos(idx+1,1), ionPos(idx+1,2), ionPos(idx+1,3), 'MarkerEdgeColor', 'm')
axis equal
% title('ion positions')
% saveas(gcf, [path folder 'ion_pos.jpg']);
% 
figure
scatter(time, ionPos(:,3))
xlabel('t', 'FontSize', 20)
ylabel('x_z', 'FontSize', 20)
saveas(gcf, [path folder 'pos_vel.jpg']);
set(findobj(gcf,'type','axes'),'FontSize',15);
% saveas(gcf, [path folder 'time_ionpos_z.jpg'])
% 
% figure
% scatter(time, vel(:,3))
% xlabel('t')
% ylabel('V_z')
% title([name ' z velocity vs time'])
% saveas(gcf, [path folder 'time_ionvel_z.jpg']);
% 
figure
scatter(ionPos(:,3), vel(:,2))
xlabel('z position', 'FontSize', 20)
ylabel('v_z', 'FontSize', 20)
saveas(gcf, [path folder 'pos_vel.jpg']);
set(findobj(gcf,'type','axes'),'FontSize',15);

end