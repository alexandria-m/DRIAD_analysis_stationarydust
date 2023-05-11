% This scripts plots the outside, inside, outside+inside+dust potential.
% The outside and inside are plotted over only the entire simulation region, but
% the outside+inside+dust is plotted over both the entire simulation as
% well as a region set by 1.5*debye above and below the dust.
%% Outside Potential
NUMR = RESX/2;
dr = RAD_CYL / (NUMR-2);
dz = 2 * HT_CYL / (RESZ - 1);
pot_out_data = csvread([path folder dataset{d} name '_outside_potential.txt']);
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

pot_out = interp2(Zout,Xout,pot_out,Z,X);
max_out = max(max(pot_out));
min_out = min(min(pot_out));

figure(6)
colormap cool
levout = linspace(-min_out, -max_out, RESZ);
contourf(Z*1e3,X*1e3,-pot_out, levout, 'Linestyle', 'none')
axis equal
colorbar
%caxis([min_out max_out])

xlabel('x (mm)', 'FontWeight', 'bold', 'FontSize', 20);
ylabel('z (mm)', 'FontWeight', 'bold', 'FontSize', 20);
set(findobj(gcf,'type','axes'),'FontSize',15);
title('outside potential', 'FontWeight', 'bold', 'FontSize', 20)

%title(dataset{d}, 'Interpreter', 'none', 'Units', 'normalized', 'Position', [.5, 0, 1])
%% Inside Potential
avg_pot = mean(pot(:,:,5:end),3);
max_in = max(max(avg_pot));
min_in = min(min(avg_pot));

figure(7)
colormap cool
levin = linspace(min_in, max_in, RESZ);
contourf(Z*1e3,X*1e3,avg_pot, levin, 'Linestyle', 'none')
axis equal
colorbar
%caxis([min_in max_in])

xlabel('x (mm)', 'FontWeight', 'bold', 'FontSize', 20);
ylabel('z (mm)', 'FontWeight', 'bold', 'FontSize', 20);
set(findobj(gcf,'type','axes'),'FontSize',15);
title('inside potential', 'FontWeight', 'bold', 'FontSize', 20)

%% Inside + Outside + Dust Potential
epsilon0 = 8.85e-12; 
maxV = max(max(avg_pot - pot_out));
minV = min(min(avg_pot - pot_out));
    
potlev = linspace(minV, maxV, 128);
V_ions = (avg_pot-pot_out);
 
figure(8)
load('pinkblue.mat')
colormap(pinkblue)
    contourf(X*1e3, Z*1e3, V_ions, potlev, 'Linestyle', 'none')
    colorbar

    % Potential due to dust
    if(NUM_DUST >0)
    % Potential due to dust
    POS1=repmat(reshape(grid_pts,[NUM_GRID_PTS,1,2]),[1,NUM_DUST,1]);
    POS2=repmat(reshape(dustPos(:,[1 3]),[1,NUM_DUST,2]),[NUM_GRID_PTS,1,1]);
    %find the distance between every grid point and dust position
    %dist_pts_inv is num_pts x num_dust in size
    dist_pts=sqrt(sum((POS1-POS2).^2,3));
    V_temp = (1/(4*pi*epsilon0))./dist_pts;
    qq = repmat(CHARGE_DUST1, num_pts, 1);
    V_dust = qq.*V_temp;
    V_dust = sum(V_dust,2);
    V_dust = reshape(V_dust,RESX,RESZ);
    V_adj = V_ions(6,:)+V_dust(6,:);
    %levels for potential contour
    levmin = min(min(avg_pot - pot_out+V_dust));
    levmax = max(max(avg_pot - pot_out+V_dust));
    potlev = linspace(levmin, levmax, 1000);
    %plot
    contourf(X'*1e3, Z'*1e3, (V_ions+V_dust-V_adj)', potlev, 'Linestyle', 'none')
    hold on
    plot((dustPos(:,1))*1e3,(dustPos(:,3))*1e3,'k.','MarkerSize',12)
    colorbar
    caxis([minV+mean(mean(V_dust)) maxV])
    end
    
hold on
contour(X*1e3,Z*1e3,V_ions, 'LineWidth', 1, 'ShowText', 'on', 'color',[.6 .6 .6],...
        'LabelSpacing',250)
axis equal

CB = colorbar;
CB.Title.String='total potential (V)';
CB.Title.FontSize = 15;
CB.Title.FontWeight = 'bold';

xlabel('x (mm)', 'FontWeight', 'bold', 'FontSize', 20);
ylabel('z (mm)', 'FontWeight', 'bold', 'FontSize', 20);
set(findobj(gcf,'type','axes'),'FontSize',15);
%title('total potential', 'FontWeight', 'bold', 'FontSize', 20);

%% In region of dust (1.5 debye lengths above and below)
if(NUM_DUST >0)
    figure(9)
    colormap(pinkblue)
    contourf(X'*1e3, Z'*1e3, (V_ions+V_dust-V_adj)', potlev,'Linestyle', 'none')
    hold on
    plot((dustPos(:,1))*1e3,(dustPos(:,3))*1e3,'.','MarkerSize',12, 'Color', 'b')
    caxis([minV+mean(mean(V_dust)) maxV])
    limitup = find(Z*1e3 > (dustPos(1,3)+1.5*DEBYE)*1e3);
    limitdown = find(Z*1e3 < (dustPos(end,3)-1.5*DEBYE)*1e3);
    ylim([Z(limitdown(end))*1e3 Z(limitup(1))*1e3])
    axis equal
    CB = colorbar;
    CB.Title.String = 'total potential (V)';
    CB.Title.FontSize = 15;
    CB.Title.FontWeight = 'bold';
    
    xlabel('x (mm)', 'FontWeight', 'bold', 'FontSize', 20);
    ylabel('z (mm)', 'FontWeight', 'bold', 'FontSize', 20);
    set(findobj(gcf,'type','axes'),'FontSize',15);
end