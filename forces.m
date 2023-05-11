% This script calculates Fid and Fdd for each dust grain.
% %% Set-up grid
% grid_data = csvread([path folder dataset{d} name '_ion-den.txt']);
% % r_pos = sqrt(sum(dustPos(:,1:2).^2,2));
% num_pts = RESX * RESZ;
% grid_pts = grid_data(1:num_pts,:);
% potential = grid_data(num_pts+1:end,2);
% pot = reshape(potential,RESX, RESZ,round(length(potential)/num_pts));
% 
% X = reshape(grid_pts(:,1),RESX, RESZ);
% Z = reshape(grid_pts(:,2),RESX, RESZ);
%% Determine potential from outside ions
half_X = RESX/2;
pot_out_data = csvread([path folder dataset{d} name '_outside_potential.txt']);
pot_out_pts = pot_out_data(1:(NUM_GRID_PTS/2),:);
half_pot_out = pot_out_data((NUM_GRID_PTS/2)+1:NUM_GRID_PTS,1);
half_pot_out = reshape(half_pot_out,[half_X,RESZ]);
pot_out = zeros(RESX,RESZ);
pot_out(half_X+1:end,:) = half_pot_out;

for i = 1:half_X
    pot_out(i,:) = half_pot_out(end-i+1,:);
end
%% Calculate the electric field from the potentials
%outisde ions
dx = X(2,1) - X(1,1);
dz = Z(1,2) - Z(1,1);
[Ez_out, Ex_out] = gradient(pot_out, dz, dx);

%inside ions
avg_pot = mean(pot(:,:,10:end),3);
[Ez,Ex] = gradient(avg_pot,dz,dx);


% dustPos = csvread([path folder dataset{d} name '_dust-pos.txt']);
% dustPosz = dustPos(:,3);

%% calculate ion force directly from DRIAD output
% This loop is for the E_out term in the drag force
% for i = 1:NUM_DUST
%     xdiff = dustPos(i,1)-X(:,1);
%     xdiffmin = abs(min(xdiff));
%     zdiff = dustPos(i,3)-Z(1,:);
%     zdiffmin = abs(min(zdiff));
%     
%     Edust_out(i) = E_out(find(abs(xdiff) == xdiffmin), find(abs(zdiff) == zdiffmin));
%    Edust_out(i) = Ez_out(64, 128);
%    QE(i) = CHARGE_DUST1(i)*Edust_out(i);
% end

% Read in the data for the direct and orbit contributions to Fid
iondustacc = csvread([path folder dataset{d} name '_ion_on_dust_acc.txt']);
ax = iondustacc(:,1); ay = iondustacc(:,2); az = iondustacc(:,3);
mx = iondustacc(:,4); my = iondustacc(:,5); mz = iondustacc(:,6);

for i = 1:NUM_DUST
    momindx(:,i) = mx(i:NUM_DUST:end,:);
    momindy(:,i) = my(i:NUM_DUST:end,:);
    momindz(:,i) = mz(i:NUM_DUST:end,:);
end
momavgx = mean(momindx(length(momindx)/2:end,:), 1);
momavgy = mean(momindy(length(momindy)/2:end,:), 1);
momavgz = mean(momindz(length(momindz)/2:end,:), 1);

for i = 1:NUM_DUST
    accindx(:,i) = ax(i:NUM_DUST:end,:);
    accindy(:,i) = ay(i:NUM_DUST:end,:);
    accindz(:,i) = az(i:NUM_DUST:end,:);
end
accavgx = mean(accindx(length(accindx)/2:end,:), 1);
accavgy = mean(accindy(length(accindy)/2:end,:), 1);
accavgz = mean(accindz(length(accindz)/2:end,:), 1);

%fdrag = accavg + momavg + QE;

%% calculation of the dust dust force
kc=8.99e9; % Coulomb constant in mks units  
POS1 = repmat(reshape(dustPos, [NUM_DUST,1,3]),[1,NUM_DUST,1]);
POS2 = repmat(reshape(dustPos, [1,NUM_DUST,3]),[NUM_DUST,1,1]);
disp = POS1-POS2;
dist = sqrt(sum(disp.^2, 3));
force_dd = kc.*(CHARGE_DUST1'.*CHARGE_DUST1)./dist.^3;
force_dd = repmat(force_dd, [1,1,3]).*disp;
force_dd = squeeze(nansum(force_dd,2));

%% Plot the forces
if(NUM_DUST > 1) %Fdd will be zero for one dust grain
    figure(3)
    hold on
    plot(dustPos1*1e3, force_dd(:,1), '*', 'MarkerSize', 15, 'Color', 'm') 
    plot(dustPos1*1e3, force_dd(:,2), '+', 'MarkerSize', 15, 'Color', 'k') 
    plot(dustPos1*1e3, force_dd(:,3), '.', 'MarkerSize', 15, 'Color', 'c')

    legend({'x', 'y', 'z'}, 'Location', 'bestoutside', 'FontSize', 15)
    ylabel('F_{dd}', 'FontWeight', 'bold', 'FontSize', 20);
    xlabel('z (mm)', 'FontWeight', 'bold', 'FontSize', 20);
    set(findobj(gcf,'type','axes'),'FontSize',15);
end

figure(4)
hold on
plot(dustPos1*1e3, accavgx, '*', 'MarkerSize', 15, 'Color', 'm')
plot(dustPos1*1e3, accavgy, '+', 'MarkerSize', 15, 'Color', 'k')
plot(dustPos1*1e3, accavgz, '.', 'MarkerSize', 15, 'Color', 'c')

legend({'x', 'y', 'z'}, 'Location', 'bestoutside', 'FontSize', 15)
ylabel('accDustIon', 'FontWeight', 'bold', 'FontSize', 20);
xlabel('z (mm)', 'FontWeight', 'bold', 'FontSize', 20);
set(findobj(gcf,'type','axes'),'FontSize',15);

figure(5)
hold on
plot(dustPos1*1e3, momavgx, '*', 'MarkerSize', 15, 'Color', 'm')
plot(dustPos1*1e3, momavgy, '+', 'MarkerSize', 15, 'Color', 'k')
plot(dustPos1*1e3, momavgz, '.', 'MarkerSize', 15, 'Color', 'c')

legend({'x', 'y', 'z'}, 'Location', 'bestoutside', 'FontSize', 15)
ylabel('momDustIon', 'FontWeight', 'bold', 'FontSize', 20);
xlabel('z (mm)', 'FontWeight', 'bold', 'FontSize', 20);
set(findobj(gcf,'type','axes'),'FontSize',15);
