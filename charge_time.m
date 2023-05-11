% This script plots the charge at each timestep. If the dust moves, timesteps = 
% dust_dt. If the dust doesn't move, timesteps = ION_TIME_STEP.

figure(2)
%% Read in dust charge
if(NUM_DUST > 1)
    time = CHARGE_DUST(1:end,1);
    CHARGE_DUST(end,:) = [];
end
qe = -1.602e-19;

%% plot dust charge
colors = [purple1; red1; blue1; orange1; red2; green1; blue2; orange2; purple2; green2;];
for i = 1:NUM_DUST
    plot((CHARGE_DUST(:,i)/qe)/1000, '.', 'Color', colors(i,:))
    hold on
end

%ylim([(min(CHARGE_DUST(:,1)/qe) - 300) max(CHARGE_DUST(:,1)/qe) + 300]/1000)
xlabel('Timesteps [10^{-9} s]', 'FontSize', 20) 
ylabel('Charge \times 10^3 [e^{-1}]', 'FontSize', 20)
set(findobj(gcf,'type','axes'),'FontSize',15);

figure
x = 1:10;
y = rand(10, 10);