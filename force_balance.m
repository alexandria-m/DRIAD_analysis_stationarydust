clear
format short
%Setting constants
kc=8.99e9; % Coulomb constant in mks units
epsilon0 = 8.85e-12; % permitivity of free space in mks unites
kb = 1.38e-23; % Boltzmann constant in mks units
qe = -1.609e-19; %charge of an electron

%specify colors for figures
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
folder = 'zigzag/';
%Name of each dataset, can read in multiple at a time
dataset = {'80deg_xz_3/'};
%Name of run (must be the same for each dataset)
name = {'80deg_xz'};

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
    
    %read in the positions of the dust grains
    dustPos = csvread([path folder dataset{d} name{d} '_dust-pos.txt']);
    dustPos1 = dustPos(:,3); %z-position
    dustPos2 = dustPos(:,1); %x-position
    dustPos3 = dustPos(:,2); %y-position
    
    %% Force due to gravity
    g = 9.8;
    MASS_DUST = repmat(MASS_DUST,[1,NUM_DUST,1]);
    force_grav = -g.*MASS_DUST;

    %% calculation of the dust dust force
    %read in and determine the average charge on the dust particles
    CHARGE_DUST1 = csvread([path folder dataset{d} name{d} '_dust-charge.txt']);
    CHARGE_DUST1(:,end) = [];%get rid of last empty column
    CHARGE_DUST = mean(CHARGE_DUST1(100:end,:)); %average the charge
    
    %calculate distance between each dust grain (r)
    POS1 = repmat(reshape(dustPos, [NUM_DUST,1,3]),[1,NUM_DUST,1]);
    POS2 = repmat(reshape(dustPos, [1,NUM_DUST,3]),[NUM_DUST,1,1]);
    disp = POS1-POS2; %x = (x1-x2)
    dist = sqrt(sum(disp.^2, 3)); %r = (x^2 + y^2 + z^2)^1/2
    
    %force between dust grains is a coulomb force
    force_dd = kc.*(CHARGE_DUST'.*CHARGE_DUST)./dist.^3;
    force_dd = repmat(force_dd, [1,1,3]).*disp;
    force_dd = squeeze(nansum(force_dd,2));

    %% Calculate ion force from Maxwell stress tensor
    MSTforce = csvread([path folder dataset{d} name{d} '_pts.txt']);
    fdx = MSTforce(:,1); fdy = MSTforce(:,2); fdz = MSTforce(:,3);
    
    %read in the force from the ions for each dust grain at each timestep
    for i = 1:NUM_DUST
        fx(:,i) = fdx(i:NUM_DUST:end,:);
        fy(:,i) = fdy(i:NUM_DUST:end,:);
        fz(:,i) = fdz(i:NUM_DUST:end,:);
    end
    
    %calculate the average force from the ions
    fdragx = mean(fx(length(fx)/2:end,:), 1);
    fdragy = mean(fy(length(fy)/2:end,:), 1);
    fdragz = mean(fz(length(fz)/2:end,:), 1);

    %% Total electric field in z
    E_neededz = ((force_dd(:,3)' + fdragz  + force_grav)./CHARGE_DUST); %MST
    E_totz = -(E_neededz);

    %Plot the z electric field at each dust grain for MST method
    symbols = 'o';
    figure(1)
    set(gcf, 'Position', [150 550 700 350], 'color', 'w') %one plot
    plot(dustPos1,E_totz,  ...
            'Marker',symbols,'LineStyle','none',...
            'MarkerEdgeColor',mycolor,...
            'MarkerFaceColor',mycolor,...
            'MarkerSize',10, 'LineWidth',1,'Color',mycolor);
        
    if fit == 1
        P = polyfit(dustPos1,E_totz,1);
        yfit = polyval(P,dustPos1);
        hold on;
        plot(dustPos1,yfit,'r');
        eqn = string("z-fit (MST): Ez = " + P(1)) + "z + " + string(P(2))
    else
        Ezp = EZ_A*dustPos1.^4 + EZ_B*dustPos1.^3 + EZ_C*dustPos1.^2 + EZ_D*dustPos1 + EZ_E;
        hold on
        plot(dustPos1, Ezp, 'r', 'LineWidth', 2)
    end

    set(findobj(gcf,'type','axes'),'FontSize',20);
    xlabel('z (m)', 'FontWeight', 'bold', 'FontSize', 20);
    ylabel('E_z (V/m)', 'FontWeight', 'bold', 'FontSize', 20);

    %% Force balance in z
     if (EZ_A == 0) %constant electric field
        f_efield = E_FIELD.*CHARGE_DUST;
        total_forcez = force_dd(:,3)' + fdragz  + force_grav + f_efield;
        
     elseif (EZ_A > 0) %linear electric field
        f_efield = (EZ_FIELD + EZ_A * dustPos1).* CHARGE_DUST';
        total_forcez = force_dd(:,3)' + fdragz  + force_grav + f_efield';
     else %polynomial electric field
        f_efield = EZ_A*dustPos1.^4 + EZ_B*dustPos1.^3 + EZ_C*dustPos1.^2 + EZ_D*dustPos1 + EZ_E;
        total_forcez = force_dd(:,3)' + fdragz  + force_grav + f_efield';
     end

    %find percentages
    f_t_individualz = (total_forcez'./force_grav')*100;
    f_t_avgz = mean(abs((f_t_individualz)))

    %save all the forces in one place
    all_forcesz = zeros(NUM_DUST,4);
    all_forcesz(:,1) = force_dd(:,3); all_forcesz(:,2) = fdragz';  
    all_forcesz(:,3) = f_efield; all_forcesz(:,4)=force_grav;
    

    %% Total electric field in x
    E_neededx = ((force_dd(:,1)' + fdragx)./CHARGE_DUST);
    E_totx = -(E_neededx);

    %plot the x electric field at each dust for MST method
    symbols = 'o';
    figure(2)
    set(gcf, 'Position', [150 550 700 350], 'color', 'w')
    plot(dustPos2,E_tot1x,  ...
            'Marker',symbols,'LineStyle','none',...
            'MarkerEdgeColor',mycolor,...
            'MarkerFaceColor',mycolor,...
            'MarkerSize',10, 'LineWidth',1,'Color',mycolor);

    set(findobj(gcf,'type','axes'),'FontSize',20);
    xlabel('x (m)', 'FontWeight', 'bold', 'FontSize', 20);
    ylabel('E_x (V/m)', 'FontWeight', 'bold', 'FontSize', 20);
    
    %fit field with linear equation with y-intercept at 0
    dlm = fitlm(dustPos2,E_tot1x,'Intercept',false);
    hold on;
    lsline
    eqn = string("x-fit (MST): Ex = " + dlm.Coefficients.Estimate) + "x"
    
    %% Force balance in x
    if (EX_B > 0) %if we are applying external Ex
        f_efield = (EX_FIELD + EX_B * dustPos2).* CHARGE_DUST';
        total_forcex = force_dd(:,1)' + fdragx + f_efield';
    else
        total_forcex = force_dd(:,1)' + fdragx;
    end
    
    %find percentages
    f_t_individualx = (total_forcex'./force_grav')*100;
    f_t_avgx = mean(abs((f_t_individualx)))

    all_forcesx = zeros(NUM_DUST,3);
    all_forcesx(:,1) = force_dd(:,1); all_forcesx(:,2) = fdragx';
    
    if(EX_B > 0) %only if we are applying external Ex
        all_forcesx(:,4) = f_efield';
    end

end

export_fig(figure(1), [path folder dataset{d} '/dust_ez'], '-png','-r300', '-transparent')
export_fig(figure(2), [path folder dataset{d} '/dust_ex'], '-png','-r300', '-transparent')