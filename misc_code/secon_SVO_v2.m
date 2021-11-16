function scores = secon_SVO_v2( secondary_data )

%% own analysis script for SVO secondary items

%% get list of all variables
% this is hard-coded in the m-file EP_gambles_v2
% so simply run it to get the relevant variables
EP_gambles_v2;

%% get endpoints for self
% get the slopes and intercepts for the 2ndary items
% compute these from the outer items
% slope = (O_left - O_right) / (S_left - S_right)
% intercept = O_left - slope*S_left

SL_endpoints_self = zeros( 9, 2 );
SL_slope_intercep = zeros( 9, 2 );

for i_i = 1:9 % these are the 2ndary items
    SL_ind = 6+i_i; % here "6+" because, we want 2ndary items
    SL_endpoints_self(i_i,1:2)  = SL(1,[1,9],SL_ind);
    SL_slope_intercep(i_i,1)    = (SL(2,1,SL_ind) - SL(2,9,SL_ind)) / ...
                                  (SL(1,1,SL_ind) - SL(1,9,SL_ind));
    SL_slope_intercep(i_i,2)    = SL(2,1,SL_ind) - SL_slope_intercep(i_i,1) * SL(1,1,SL_ind);
end

%% get ranges with 101 entries & minimum difference on these
item_range_self     = zeros( 9, 101 );
item_range_othe     = zeros( 9, 101 );
ineq_index          = zeros( 9,   1 );
for i_i = 1:9
    item_range_self(i_i,:) = linspace( SL_endpoints_self(i_i,1), SL_endpoints_self(i_i,2), 101);
    item_range_othe(i_i,:) = item_range_self(i_i,:) * SL_slope_intercep(i_i,1) + SL_slope_intercep(i_i,2);
    % item(6+count,item_range_self(count,:));
    [~, ineq_index(i_i,1)] = min( abs( item_range_self(i_i,:) - item_range_othe(i_i,:) ) );
end
ineq_index_t = ineq_index';

%% get data 
% example_data = [74 94;94 96;81 81;90 100;85 85;100 90;63 88;100 90;90 100];
% secondary_data = example_data;

% init
secondary_data_option = nan(1,9);
% for 1:9 options 
for i_i = 1:9    
    [~, secondary_data_option(1,i_i)] = min( abs( secondary_data(i_i,1) - item_range_self(i_i,:) ) );
end
max_diff = 101-1;

%% joint
ideal_joint_maximizer=[
    101
    NaN
    1
    101
    NaN
    101
    NaN
    1
    1]';

%% Calculate differences
scores(1,1) = nanmean(( abs( secondary_data_option - ineq_index_t )./ max_diff ), 2 ); % Ineq_diff
scores(1,2) = nanmean(( abs( secondary_data_option - ideal_joint_maximizer)./ max_diff), 2 ); % joint_diff
scores(1,3) = scores(1,1) + scores(1,2);
scores(1,4) = scores(1,1) / scores(1,3);
scores(1,5) = scores(1,2) / scores(1,3);