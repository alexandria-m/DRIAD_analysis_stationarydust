clear
format short
%Setting constants
kc=8.99e9; % Coulomb constant in mks units
epsilon0 = 8.85e-12; 
kb = 1.38e-23;

purple = [170 60 195]/255;
red = [200 70 90]/255;
green = [70 200 90]/255;
blue = [100 100 195]/255;
mycolor = blue;
fit = 1; %1 for linear, 2 for polynomial
%% Set path to data
%The path to your data
path = '/Users/alexandriamendoza/Documents/Plasma/';
%Folder your data is in
folder = 'changing_density/debugging/charging/';
%Name of each dataset, can read in multiple at a time
dataset = {'constant_shifted_full/'};
%Name of run (must be the same for each dataset)
name = {'constant'};

%%
for d = 1:length(dataset)
    close all
    import_debug_file
    fprintf(['dataset: ' dataset{d} '\n']);
    fid = fopen([path folder dataset{d} name{d} '_params.txt']);
    while ~feof(fid)
        tline = fgetl(fid);
        %disp(tline);
        [~,q] = find(tline == '%');
        tempval = deblank(tline(1:q-1));
        tempvar = tline(q+2:end);
        eval([tempvar ' = ' tempval ';']);
    end
    fclose(fid);
    E_FIELD = EZ_FIELD;
    EX_FIELD = EX_MAX;
    qe = -1.609e-19;
    dustPos = csvread([path folder dataset{d} name{d} '_dust-pos.txt']);
    dustPos1 = dustPos(:,3); % z-position
    dustPos2 = dustPos(:,1); % x-position
    %% Force due to gravity
    g = 9.8;
    MASS_DUST = repmat(MASS_DUST,[1,NUM_DUST,1]);
    force_grav = -g.*MASS_DUST;

    %% calculation of the dust dust force
    CHARGE_DUST1 = csvread([path folder dataset{d} name{d} '_dust-charge.txt']);
    CHARGE_DUST1(:,end) = [];%get rid of last empty column
    CHARGE_DUST = mean(CHARGE_DUST1(100:end,:)); %average the charge
    POS1 = repmat(reshape(dustPos, [NUM_DUST,1,3]),[1,NUM_DUST,1]);
    POS2 = repmat(reshape(dustPos, [1,NUM_DUST,3]),[NUM_DUST,1,1]);
    disp = POS1-POS2;
    dist = sqrt(sum(disp.^2, 3));
    force_dd = kc.*(CHARGE_DUST'.*CHARGE_DUST)./dist.^3;
    force_dd = repmat(force_dd, [1,1,3]).*disp;
    force_dd = squeeze(nansum(force_dd,2));
    %% Calculate outside ion force
    grid_data = csvread([path folder dataset{d} name{d} '_ion-den.txt']);
    num_pts = RESX * RESZ;
    grid_pts = grid_data(1:num_pts,:);

    X = reshape(grid_pts(:,1),RESX, RESZ);
    Z = reshape(grid_pts(:,2),RESX, RESZ);
    NUMR = RESX/2;
    num_pts = RESX * RESZ;
    dr = RAD_CYL / (NUMR-2);
    dz = 2 * HT_CYL / (RESZ - 1);
    pot_out_data = csvread([path folder dataset{d} name{d} '_outside_potential.txt']);
    pot_out_pts = pot_out_data(1:(num_pts/2),:);
    half_pot_out = pot_out_data((num_pts/2)+1:num_pts,1);
    half_pot_out = reshape(half_pot_out,[NUMR,RESZ]);
    pot_out = [flipud(half_pot_out); half_pot_out];
    pot_out(NUMR-1:NUMR+1,:) = [];

    Xout = reshape(pot_out_pts(:,1),NUMR,RESZ);
    Xout = [-flipud(Xout); Xout];
    Xout(NUMR-1:NUMR+1,:) = [];
    Zout = reshape(pot_out_pts(:,2),NUMR,RESZ);
    Zout = [flipud(Zout); Zout];
    Zout(NUMR-1:NUMR+1,:) = [];

    pot_out1 = interp2(Zout,Xout,pot_out,Z,X);
    pottest = scatteredInterpolant(X,Z, pot_out);
    
    % Calculate the electric field from the outside potential
    dx = X(2,1) - X(1,1);
    dz = Z(1,2) - Z(1,1);
    [Ez_out, Ex_out] = gradient(pot_out1, dz, dx);

    % This loop is for the E_out force in z and x
    % Find the electric field from the outside ions at each dust grain by
    % finding the minimum distance between the dust grain and the grid,
    % then finding what the outside electric field is at that grid value 
    for i = 1:NUM_DUST
        xdiff = X(:,1) - dustPos(i,1);
        [xVal, xIdx] = min(abs(xdiff)); %minimum difference in x
        zdiff = Z(1,:) - dustPos(i,3);
        [zVal, zIdx] = min(abs(zdiff)); %minimum difference in z
        pionout(i) = interp2(Zout, Xout, pot_out,dustPos(i,1),dustPos(i,3));
        Ezionout(i) = interp2(Ez_out, dustPos(i,1), dustPos(i,3));
        %find the outside electric field in z and multiply by charge to get
        %force
        Ezdust_out(i) = Ez_out(xIdx, zIdx);
        QEz(i) = CHARGE_DUST1(i)*Ezdust_out(i);
        
        %find the outside electric field in x and multiply by charge to get
        %force
        Exdust_out(i) = Ex_out(xIdx, zIdx);
        QEx(i) = CHARGE_DUST1(i)*Exdust_out(i);
    end

    %% Calculate ion force from Maxwell stress tensor
    MSTforce = csvread([path folder dataset{d} name{d} '_pts.txt']);
 %   pts = MSTforce(1:100,:); %points on the sphere
%     rh = csvread([path folder datset{d} name{d} '_radius.txt']);
%    MSTforce = MSTforce(101:end,:);
    fd1x = MSTforce(:,1); fd1y = MSTforce(:,2); fd1z = MSTforce(:,3);
    
    for i = 1:NUM_DUST
        %find the average radius of the sphere around each dust grain
%         rH(i) = mean(rh(i:NUM_DUST:end,:);
        %find the force on each dust grain from inside ions
        fx(:,i) = fd1x(i:NUM_DUST:end,:);
        fy(:,i) = fd1y(i:NUM_DUST:end,:);
        fz(:,i) = fd1z(i:NUM_DUST:end,:);
    end
    
    %find the average force on each dust grain from the inside ions
    fdrag2x = mean(fx(length(fx)/2:end,:), 1);
    fdrag2y = mean(fy(length(fy)/2:end,:), 1);
    fdrag2z = mean(fz(length(fz)/2:end,:), 1);
    
    %find the average force on each dust grain from the outside ions
    %first, shift the sphere to around the dust grain
    for i = 1:NUM_DUST
%         sphPtsx = rH(i)*pts(:,1) + dustPos(i,1);
%         sphPtsy = rH(i)*pts(:,2) + dustPos(i,2);
%         sphPtsz = rH(i)*pts(:,3) + dustPos(i,3);
        Ezout = interp2(Ez_out,dustPos(i,1),dustPos(i,3));
    end
    
    qi = 1.602e-19; %charge on an ion
    read_potential_density;
    avg_pot = mean(pot(:,:,10:end),3);

    %Get the collection radius for the ions
    phi_dust = CHARGE_DUST*8.89e9/RAD_DUST;
    vs_sq = 8*kb*TEMP_ION/pi/MASS_SINGLE_ION + (MACH*SOUND_SPEED)^2;
    b_c = RAD_DUST*sqrt(1-2*qi*phi_dust/MASS_SINGLE_ION/vs_sq);

    fdrag1z = calc_force_from_stress_tensor(path, folder, name, dataset, d, X,Z, b_c)
    %% Calculate ion force from DRIAD
    % Read in the data for the direct and orbit contributions to Fid
%     iondustacc = csvread([path folder dataset{d} name{d} '_ion_on_dust_acc.txt']);
%     ax = iondustacc(:,1); ay = iondustacc(:,2); az = iondustacc(:,3);
%     mx = iondustacc(:,4); my = iondustacc(:,5); mz = iondustacc(:,6);
% 
%     for i = 1:NUM_DUST
%         momindx(:,i) = mx(i:NUM_DUST:end,:);
%         momindy(:,i) = my(i:NUM_DUST:end,:);
%         momindz(:,i) = mz(i:NUM_DUST:end,:);
%     end
%     momavgx = mean(momindx(length(momindx)/2:end,:), 1);
%     momavgy = mean(momindy(length(momindy)/2:end,:), 1);
%     momavgz = mean(momindz(length(momindz)/2:end,:), 1);
% 
%     for i = 1:NUM_DUST
%         accindx(:,i) = ax(i:NUM_DUST:end,:);
%         accindy(:,i) = ay(i:NUM_DUST:end,:);
%         accindz(:,i) = az(i:NUM_DUST:end,:);
%     end
%     accavgx = mean(accindx(length(accindx)/2:end,:), 1);
%     accavgy = mean(accindy(length(accindy)/2:end,:), 1);
%     accavgz = mean(accindz(length(accindz)/2:end,:), 1);
% 
%     fdrag2x = accavgx + momavgx;
%     fdrag2y = accavgy + momavgy;
%     fdrag2z = accavgz + momavgz;

    %% calculation of the dust dust force
    kc=8.99e9; % Coulomb constant in mks units  
    POS1 = repmat(reshape(dustPos, [NUM_DUST,1,3]),[1,NUM_DUST,1]);
    POS2 = repmat(reshape(dustPos, [1,NUM_DUST,3]),[NUM_DUST,1,1]);
    disp = POS1-POS2;
    dist = sqrt(sum(disp.^2, 3));
    force_dd = kc.*(CHARGE_DUST'.*CHARGE_DUST)./dist.^3;
    force_dd = repmat(force_dd, [1,1,3]).*disp;
    force_dd = squeeze(nansum(force_dd,2));

    %% Total electric field in z
    if NUM_DUST == 1 
        E_needed1z = ((fdrag1z  + force_grav)./CHARGE_DUST); %MATLAB
        E_needed2z = ((fdrag2z  + force_grav)./CHARGE_DUST); %DRIAD
        E_tot1z = -(E_needed1z);
        E_tot2z = -(E_needed2z);
    else
        E_needed1z = ((force_dd(:,3)' + fdrag1z  + force_grav)./CHARGE_DUST); %MATLAB
        E_needed2z = ((force_dd(:,3)' + fdrag2z  + force_grav)./CHARGE_DUST); %DRIAD
        E_tot1z = -(E_needed1z);
        E_tot2z = -(E_needed2z);
    end

    symbols = 'o';
    figure(1)
    set(gcf, 'Position', [150 550 700 350], 'color', 'w') %one plot
    zmm = dustPos1*1e3; 
    E_Vmm = E_tot1z/1e3;
    plot(dustPos1,E_tot1z,  ...
            'Marker',symbols,'LineStyle','none',...
            'MarkerEdgeColor',mycolor,...
            'MarkerFaceColor',mycolor,...
            'MarkerSize',10, 'LineWidth',1,'Color',mycolor);
        
    if fit == 1
        P = polyfit(dustPos1,E_tot1z,1);
        yfit = polyval(P,dustPos1);
        hold on;
        plot(dustPos1,yfit,'r');
        eqn = string("z-fit (MATLAB): Ez = " + P(1)) + "z + " + string(P(2))
    else
        Ezp = EZ_A*dustPos1.^4 + EZ_B*dustPos1.^3 + EZ_C*dustPos1.^2 + EZ_D*dustPos1 + EZ_E;
        hold on
        plot(dustPos1, Ezp, 'r', 'LineWidth', 2)
    end

    set(findobj(gcf,'type','axes'),'FontSize',20);
    title('MATLAB:', 'Units', 'normalized', 'Position', [.5, 0, 1])
    xlabel('z (m)', 'FontWeight', 'bold', 'FontSize', 20);
    ylabel('E_z (V/m)', 'FontWeight', 'bold', 'FontSize', 20); 

    figure(2)
    set(gcf, 'Position', [150 550 700 350], 'color', 'w') %one plot
    zmm = dustPos1*1e3; 
    E_Vmm = E_tot2z/1e3;
    plot(dustPos1,E_tot2z,  ...
            'Marker',symbols,'LineStyle','none',...
            'MarkerEdgeColor',mycolor,...
            'MarkerFaceColor',mycolor,...
            'MarkerSize',10, 'LineWidth',1,'Color',mycolor);
    if fit == 1
        P = polyfit(dustPos1,E_tot2z,1);
        yfit = polyval(P,dustPos1);
        hold on;
        plot(dustPos1,yfit,'r');
        eqn = string("z-fit (DRIAD): Ez = " + P(1)) + "z + " + string(P(2))
    else
        Ezp = EZ_A*dustPos1.^4 + EZ_B*dustPos1.^3 + EZ_C*dustPos1.^2 + EZ_D*dustPos1 + EZ_E;
        hold on
        plot(dustPos1, Ezp, 'r', 'LineWidth', 2)
    end
    
    set(findobj(gcf,'type','axes'),'FontSize',20);
    title('DRIAD', 'Units', 'normalized', 'Position', [.5, 0, 1])
    xlabel('z (m)', 'FontWeight', 'bold', 'FontSize', 20);
    ylabel('E_z (V/m)', 'FontWeight', 'bold', 'FontSize', 20); 

    %% Force balance in z
    if (EZ_A == 0)
        f_efield = E_FIELD.*CHARGE_DUST;
        total_force1z = force_dd(:,3)' + fdrag1z  + force_grav + f_efield;
        total_force2z = force_dd(:,3)' + fdrag2z  + force_grav + f_efield;
    elseif (EZ_B > 0)
        f_efield = (EZ_FIELD + EZ_A * dustPos1).* CHARGE_DUST';
        total_force1z = force_dd(:,3)' + fdrag1z  + force_grav + f_efield';
        total_force2z = force_dd(:,3)' + fdrag2z  + force_grav + f_efield'; 
    end
    f_t_individual1z = (total_force1z'./force_grav')*100;
    f_t_avg1z = mean(abs((f_t_individual1z)))

    f_t_individual2z = (total_force2z'./force_grav')*100;
    f_t_avg2z = mean(abs((f_t_individual2z)))
    
    all_forcesz = zeros(NUM_DUST,5);
    all_forcesz(:,1) = force_dd(:,3); all_forcesz(:,2) = fdrag1z'; 
    all_forcesz(:,3) = fdrag2z; all_forcesz(:,4)=f_efield;
    all_forcesz(:,5) = force_grav;

end

export_fig(figure(1), [path folder dataset{d} '/dust_ez_mst'], '-png','-r300', '-transparent')
export_fig(figure(2), [path folder dataset{d} '/dust_ez_direct'], '-png','-r300', '-transparent')
export_fig(figure(3), [path folder dataset{d} '/vertical_ez'], '-png','-r300', '-transparent')