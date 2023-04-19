

%% read in specific data from debug file
% import_without_A;
grid_data = csvread([path folder dataset{d} name '_ion-den.txt']);

%%
%determine how many grid points there are
num_pts = NUM_GRID_PTS;
grid_pts = grid_data(1:num_pts,:);
density =grid_data(num_pts+1:end,1);
potential = grid_data(num_pts+1:end,2);

X = reshape(grid_pts(:,1),RESX,RESZ);
Z = reshape(grid_pts(:,2),RESX,RESZ);
%%
den = reshape(density,RESX,RESZ,round(length(density)/num_pts));
pot = reshape(potential,RESX,RESZ,round(length(potential)/num_pts));
