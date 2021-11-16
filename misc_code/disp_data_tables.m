%==========================================================================
% Display tables  
%
% T1: Number of participants for each svo category (primary scores)
% T2: Number of participants for each prosocial category (secondary scores)
% T3: Group mean SVO (primary & secondary scores) 
% T4: Group mean parameter estimates for Inequality & Joint-gain model
% 
% T1 & T2 will be displayed for all samples/subgroups
% T3 & T4 will be displayed for patient & control groups (control sample &
% online subgroup)
%==========================================================================

% sample name (+1 because dataset indices start with 0)
sample = sampleNames{which_dataset+1};

%% primary SVO categories

rowNamesT1 = {'Altruistic', 'Prosocial', 'Individualistic', 'Competitive', 'Total N'};
varNameT1 = {'SVO_categories_primary'};

% svo primary categories
alt = sum(SVO_out_acker(:,2) ==1);
pro = sum(SVO_out_acker(:,2) ==2);
ind = sum(SVO_out_acker(:,2) ==3);
com = sum(SVO_out_acker(:,2) ==4);
% total N
totalNpri = alt + pro + ind + com;

T1 = table([alt; pro; ind; com; totalNpri],...
    'VariableNames', varNameT1,...
    'RowNames', rowNamesT1);

disp(sample)
disp(T1)
disp(' ')

%% secondary SVO categories

rowNamesT2 = {'Prosocial - inequality averse', 'Prosocial - joint gain maximizing', 'Total N'};
varNameT2 = {'SVO_categories_secondary'};

% svo secondary categories
ia = sum(SVO_out_acker(:,2) == 2 & (SVO_out_acker(:,8) ==1));
jg = sum(SVO_out_acker(:,2) == 2 & (SVO_out_acker(:,8) ==2));
% total N
totalNsec = sum(incl_pro_sn);

T2 = table([ia; jg; totalNsec],...
    'VariableNames', varNameT2,...
    'RowNames', rowNamesT2);

disp(sample)
disp(T2)
disp(' ')

%% Group mean SVO scores

if which_dataset == 2 || which_dataset == 5
    
    % SVO score/primary scores
    priMean = round(mean(SVO_out_acker(:,1)),2); 
    priSD = round(std(SVO_out_acker(:,1)),2);   
    priN = length(SVO_out_acker(:,1));           
    
    % Prosocial motivation score/secondary scores
    % incl_pro_sn = index for prosocial participants
    secMean = round(mean(SVO_out_acker(incl_pro_sn,9)),2);
    secSD = round(std(SVO_out_acker(incl_pro_sn,9)),2);
    secN = length(SVO_out_acker(incl_pro_sn,9));
    
    % Inequality distance score
    ineMean = round(mean(scores_2nd(incl_pro_sn,1)),2);
    ineSD = round(std(scores_2nd(incl_pro_sn,1)),2);
    ineN = length(scores_2nd(incl_pro_sn,1));
    
    % Joint-gain distance score
    joiMean = round(mean(scores_2nd(incl_pro_sn,2)),2);
    joiSD = round(std(scores_2nd(incl_pro_sn,2)),2);
    joiN = length(scores_2nd(incl_pro_sn,2));
    
elseif which_dataset == 3 % if patient (single score, no SD) 
    
    % SVO score/primary scores
    priMean = round(SVO_out_acker(:,1),2);
    priSD = NaN;
    priN = length(SVO_out_acker(:,1));
    
    % Prosocial motivation score/secondary scores
    secMean = round(SVO_out_acker(:,9),2);
    secSD = NaN;
    secN = length(SVO_out_acker(:,9));
    
    % Inequality distance score
    ineMean = round(scores_2nd(:,1),2);
    ineSD = NaN;
    ineN = length(scores_2nd(:,1));
    
    % Joint-gain distance score
    joiMean = round(scores_2nd(:,2),2);
    joiSD = NaN;
    joiN = length(scores_2nd(:,2));
    
end

% make table if dataset == patient, control, or online subgroup
if which_dataset == 2 || which_dataset == 3 || which_dataset == 5
    
    rowNamesT3 = {'Mean', 'SD', 'N'};
    varNameT3 = {'SVO_score', 'Prosocial_motivation_score',...
        'Inequality_distance_score', 'Joint_gain_distance_score'};

    dataT3 = [priMean, secMean, ineMean, joiMean;...
        priSD, secSD, ineSD, joiSD;...
        priN, secN, ineN, joiN];
    
    T3 = array2table(dataT3,...
        'VariableNames', varNameT3,...
        'RowNames', rowNamesT3);

    disp(sample) 
    disp(T3)
    disp(' ')
end


%% JPE model parameter
% Model: Inequality & Joint-gain
% Parameter: Inequality, Joint-gain
% Samples: patient, control, online subgroup

modelT4 = 6; % model idx
parIn = 2; % inequality parameter idx 
parJg = 3; % joint-gain parameter idx

if which_dataset == 2 || which_dataset == 5
    
    % inequality parameter estimates    
    parInM = round(mean(HF_B{modelT4}(:,parIn)),3);
    parInSD = round(std(HF_B{modelT4}(:,parIn)),3);
    parInN = length(HF_B{modelT4}(:,parIn));
    
    % joint-gain parameter estimates    
    parJgM = round(mean(HF_B{modelT4}(:,parJg)),3);
    parJgSD = round(std(HF_B{modelT4}(:,parJg)),3);
    parJgN = length(HF_B{modelT4}(:,parJg));
    
elseif which_dataset == 3
    
    % inequality parameter estimates    
    parInM = round(HF_B{modelT4}(:,parIn),3);
    parInSD = NaN;
    parInN = length(HF_B{modelT4}(:,parIn));
    
    % joint-gain parameter estimates    
    parJgM = round(HF_B{modelT4}(:,parJg),3);
    parJgSD = NaN;
    parJgN = length(HF_B{modelT4}(:,parJg));
    
end

% make table if dataset == patient, control, or online subgroup
if which_dataset == 2 || which_dataset == 3 || which_dataset == 5
    
    tableNameT4 = ['Model parameter estimates: ' mod_vec{modelT4} ' model'];
    rowNamesT4 = {'Mean', 'SD', 'N'};
    varNameT4 = {'Inequality', 'Joint_gain'};
    
    dataT4 = [parInM parJgM;...
        parInSD parJgSD;...
        parInN parJgN];
    
    T4 = array2table(dataT4,...
        'VariableNames', varNameT4,...
        'RowNames', rowNamesT4);
    
    disp(sample)
    disp(tableNameT4)
    disp(T4)
    disp('')
    
end