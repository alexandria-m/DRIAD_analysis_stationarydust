symbols = '^';
mycolor = red2;

%% Calculation of the charge
qe = -1.602e-19;

% read the dust charge
CHARGE_DUST = csvread([path folder dataset{d} name '_dust-charge.txt']);
CHARGE_DUST(:,end) = [];
CHARGE_DUST1 = mean(CHARGE_DUST(100:end,:));

dustPos1 = dustPos(:,3);
%dustPosbox = (0.0060+ dustPos1); %shift the data so electrode is at zero

%% plot in coulombs
% figure
% elec = -1.602e-19; %charge of electron in coulombs
%  h2 = plot(dustPos1,CHARGE_DUST1, ...
%         'Marker',symbols,'LineStyle','-',...
%         'MarkerEdgeColor',mycolor,...
%         'MarkerFaceColor',mycolor,...
%         'MarkerSize',7, 'LineWidth',1,'Color',mycolor);
% 
% xlabel('z (mm)','Fontsize',15)
% ylabel('Charge (C)','Fontsize',15)

%% plot in electrons
figure(1)
coulomb = -6.24e18; %coulomb to electrons
h2 = plot(dustPos1,CHARGE_DUST1*coulomb, ...
         'Marker',symbols,...
         'MarkerEdgeColor',mycolor,...
         'MarkerFaceColor',mycolor,...
         'MarkerSize',7, 'LineWidth',1.5,'Color',red1);
xlabel('z (mm)','Fontsize',20)
ylabel('Charge (e)','Fontsize',20)
set(findobj(gcf,'type','axes'),'FontSize',15);
