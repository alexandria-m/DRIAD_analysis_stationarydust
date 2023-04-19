% This scripts creates a contour plot for the density, over the whole
% simulation cylinder as well as a region 1.5*debye above and below the dust.
%% Calculate density limits
maxden= max(max(final_den/DEN_FAR_PLASMA));
minden = min(min(final_den/DEN_FAR_PLASMA));
denlev = linspace(minden, maxden, RESZ);

%% Plot
figure(13)
colormap hot
contourf(X'*1e3, Z'*1e3,final_den'/DEN_FAR_PLASMA, denlev, 'Linestyle', 'none')
hold on
if(NUM_DUST >0)
    hold on
    plot((dustPos(:,1))*1e3,(dustPos(:,3))*1e3,'.','MarkerSize',12)
end
%caxis([minden maxden])
axis equal
CB=colorbar;
CB.Title.String='n_i/n_{io}';
CB.Title.FontSize = 15;
CB.Title.FontWeight = 'bold';
% set(CB.XLabel,{'String','Rotation','Position'},{'XLabel',0,[0.5 -0.01]})

xlabel('x (mm)', 'FontWeight', 'bold', 'FontSize', 20);
ylabel('z (mm)', 'FontWeight', 'bold', 'FontSize', 20);
set(findobj(gcf,'type','axes'),'FontSize',15);

%% In region of dust
figure(14)
fig1 = gcf;
colormap hot
denlev = linspace(0, maxden, 256);
contourf(X'*1e3, Z'*1e3,final_den'/DEN_FAR_PLASMA, denlev, 'Linestyle', 'none')
if(NUM_DUST >0)
    hold on
    plot((dustPos(:,1))*1e3,(dustPos(:,3))*1e3,'.','MarkerSize',12)
end
CB=colorbar;
CB.Title.String='n_i/n_{io}';
CB.Title.FontSize = 15;
CB.Title.FontWeight = 'bold';
set(CB.XLabel,{'String','Rotation','Position'},{'XLabel',0,[0.5 -0.01]})
%caxis([minden maxden])

limitup = find(Z*1e3 > (dustPos(1,3)+1.5*DEBYE)*1e3);
limitdown = find(Z*1e3 < (dustPos(end,3)-1.5*DEBYE)*1e3);
ylim([Z(limitdown(end))*1e3 Z(limitup(1))*1e3])
axis equal

xlabel('x (mm)', 'FontWeight', 'bold', 'FontSize', 20);
ylabel('z (mm)', 'FontWeight', 'bold', 'FontSize', 20);
set(findobj(gcf,'type','axes'),'FontSize',15);