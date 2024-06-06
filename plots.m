%% This file calls in all other files that plot potential, density, and ion velocity/acceleration.
clear all
%% Set path to data
%The path to your data
path = '/Users/alexandriamendoza/Documents/Plasma/';
%Folder your data is in
folder = 'zigzag/';
%Name of each dataset, can read in multiple at a time
dataset = {'80deg_xz_2/'};
%Name of run (must be the same for each dataset)
name = {'80deg_xz'};

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
    fid = fopen([path folder dataset{d} name{d} '_params.txt']);
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
    grid_data = csvread([path folder dataset{d} name{d} '_ion-den.txt']);
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
      dustPos = csvread([path folder dataset{d} name{d} '_dust-pos.txt']);
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
        export_fig(figure(1), [path folder dataset{d}  '/charge'], '-png','-r300', '-transparent');
        export_fig(figure(2), [path folder dataset{d} '/charge_time'], '-png','-r300', '-transparent');
        export_fig(figure(3), [path folder dataset{d} '/fdd'], '-png','-r300', '-transparent');
        export_fig(figure(4), [path folder dataset{d} '/fdirect'], '-png','-r300', '-transparent');
        export_fig(figure(5), [path folder dataset{d} '/forbit'], '-png','-r300', '-transparent')
        export_fig(figure(9), [path folder dataset{d} '/zoomed_potential'], '-png','-r300', '-transparent');
        export_fig(figure(14), [path folder dataset{d} '/zoomed_density'], '-png','-r300', '-transparent');
       export_fig(figure(15), [path folder dataset{d} '/fid_time'], '-png','-r300', '-transparent');
    end
    export_fig(figure(6), [path folder dataset{d} '/inside_potential'], '-png','-r300', '-transparent');
    export_fig(figure(7), [path folder dataset{d} '/outside_potential'], '-png','-r300', '-transparent');
    export_fig(figure(8), [path folder dataset{d} '/insideoutside_potential'], '-png','-r300', '-transparent');
    export_fig(figure(10), [path folder dataset{d} '/ni_nio'], '-png','-r300', '-transparent');
    export_fig(figure(11), [path folder dataset{d} '/vi'], '-png','-r300', '-transparent');
    export_fig(figure(12), [path folder dataset{d} '/acc_iz'], '-png','-r300', '-transparent');
    export_fig(figure(13), [path folder dataset{d} '/normalized_density'], '-png','-r300', '-transparent');
    export_fig(figure(16), [path folder dataset{d} '/acc_ix'], '-png','-r300', '-transparent');
end