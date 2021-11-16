%==========================================================================
% Main analyses for the 'social motives in a patient with amygdala lesions'
% paper (Doppelhofer, Hurlemann, Bach, & Korn, 2021, Neuropsychologia)
%
% To run this script change the directory, add SVO_slider.m to your path, 
% set index for sample/subgroup
%
% OUTPUT variables:
% SVO_out_acker (SVO scores, for column description see SVO_Slider.m)
% scores_2nd (own calculations of secondary score)
% LBF (Log-group Bayes factors of JPE models) 
% HF_B (regression coefficients of linear regression models) 
%
% Figures: 
% Figure 1: scatterplots & boxplots for social value orientation scores
% Figure 2: barplot for model comparison (Log-group Bayes factors) 
% 
% lisamdoppelhofer@gmail.com, christoph.w.korn@gmail.com
%==========================================================================

clear 
close all

% set directory 
dir = pwd; 
addpath(genpath(dir));

%% ========================================================================
% Sample indices 
% 0 = online sample
% 1 = student sample
% 2 = UW-control 
% 3 = UW patient
% 4 = UW + control (UW = 1. row)
% 5 = age- & gender-matched subgroup of online sample 
%==========================================================================

% set which sample (subgroup) to analyse 
which_dataset = 2;

if which_dataset == 0
    L = load('SVO_DATA_online.mat');
elseif which_dataset > 0 && which_dataset < 4
    L = load(['SVO_DATA_lab_' num2str(which_dataset) '.mat']);
elseif which_dataset == 4
    L1 = load('SVO_DATA_lab_2.mat');
    L2 = load('SVO_DATA_lab_3.mat');
    L.SVO_DA_3 = vertcat(L2.SVO_DA_3, L1.SVO_DA_3);
elseif which_dataset == 5
    L = load('SVO_DATA_online.mat');
    % age between 36 & 49 and gender = female 
    subgroupIdx = L.SVO_DA_3(:,2) > 36 & L.SVO_DA_3(:,2) < 49 & L.SVO_DA_3(:,3) == 1 ;
    L.SVO_DA_3 = L.SVO_DA_3(subgroupIdx,:);
end

sampleNames = {'Online sample',...
    'Student sample',...
    'Control sample',...
    'UW patient',...
    'UW patient & Control sample',...
    'Online subgroup'};

%% Get SVO & HF items by running this file

EP_gambles_v2

%% Initialize 

n_sub = size(L.SVO_DA_3, 1);
SVO_out_item = cell(1, n_sub);
SVO_out_item_for_ana = cell(1, n_sub);
SVO_out_acker = nan(n_sub, 14);
p_mean = nan(n_sub,2);
SVO_out = nan(n_sub,1);
scores_2nd = nan(n_sub,5);

%% ========================================================================
% Anylase SVO data 
%==========================================================================

SVO_data = L.SVO_DA_3(:,7:21);

for i_sub = 1:n_sub

    % get chosen payoffs for self and other (in SL) according to chosen
    % options (1:9, in SVO_data); one cell for each subject (15 items)
    for i_ite = 1:15 
        SVO_out_item{i_sub}(i_ite,1:2) = SL(1:2, SVO_data(i_sub,i_ite), i_ite);
    end
    
    %% calculate SVO scores via SVO_Slider.m (by K. Ackermann) 
    % using full payoff variant
    % http://ryanomurphy.com/styled-2/downloads/index.html
    
    % reshape input (from 15x2 double to 1x30 double)
    SVO_out_item_for_ana{i_sub} = reshape(SVO_out_item{i_sub}', [1,30]);
    % get output (1 row per subject, 14 columns)
    SVO_out_acker(i_sub,:) = SVO_Slider(SVO_out_item_for_ana{i_sub}); 

    %%  calculate SVO primary and secondary scores on our own
     
    % primary items 
    % svo score = inverse tangent of the ratio of mean allocation for 
    % other - 50 and mean allocation for self - 50
    p_mean(i_sub,1:2) = mean(SVO_out_item{i_sub}(1:6,:));
    SVO_out(i_sub,1) = atand((p_mean(i_sub,2)-50)/(p_mean(i_sub,1)-50));
    
    % secondary items
    % output columns: 1=abs inequality score, 2=abs joint gain score,
    % 3=sum of 1&2, 4&5 = normalized scores
    scores_2nd(i_sub,:) = secon_SVO_v2(SVO_out_item{i_sub}(7:15,:)); 

end

% compare own and SVO_Slider.m calculations
check_own_vs_acker = sum(SVO_out(:,1) - SVO_out_acker(:,1)); 

% compare secondary scores: max difference (rounding errors) 
check_own_vs_acker_sec = max(scores_2nd(:,4) - SVO_out_acker(:,9)); 

% check transitivity (1=transitive, 0=intransitive)
SVO_sum_trans = sum(SVO_out_acker(:,3)); 

%% exclude intransitive subjects

% from SVO scores 
SVO_out_acker_ori = SVO_out_acker;
excl_intra =  SVO_out_acker(:,3) == 1;
SVO_out_acker = SVO_out_acker(excl_intra,:);

% from all secondary scores 
scores_2nd_ori = scores_2nd; 
scores_2nd = scores_2nd(excl_intra,:); 

% from original data 
L.SVO_DA_3_ori = L.SVO_DA_3;
L.SVO_DA_3 = L.SVO_DA_3(excl_intra,:);

% subject number after exclusion 
n_sub = size(L.SVO_DA_3, 1);  

%% include only prosocials 
% inclusion idx for prosocial subjects according to 
% (1) primary items SVO_out_acker(:,2) 

incl_pro = SVO_out_acker(:,2) == 2; 

% (1) and (2) secondary items SVO_out_acker(:,11:12)
%   preference for both inequality aversion and joint gain maximization over
%   both individualism and altruism in the secondary items
%   see http://ryanomurphy.com/styled-2/downloads/files/SVO_Slider_Tutorial.pdf

incl_pro_sn = SVO_out_acker(:,2) == 2 & (SVO_out_acker(:,11) + SVO_out_acker(:,12) == 3);

%% PLOT SVO primary score and secondary scores (boxplots) 

figTitles = {'SVO score','Prosocial motivation score',...
    'Inequality distance score', 'Joint-gain distance score'};
col = [0 0.4470 0.7410]; 

if which_dataset ~= 3
    fig1 = figure;
    for iFig = 1:4
        
        if iFig == 1
            data = SVO_out_acker(:,1);
        elseif iFig == 2
            data = SVO_out_acker(incl_pro_sn,9);
        elseif iFig == 3
            data = scores_2nd(incl_pro_sn,1);
        elseif iFig == 4
            data = scores_2nd(incl_pro_sn,2);
        end
        
        subplot(2,2,iFig)
        scatter(ones(length(data),1) - 0.15, data ,...
            'MarkerFaceColor', col, 'MarkerFaceAlpha', 0.3, 'MarkerEdgeColor', col,...
            'jitter', 'on', 'jitterAmount', 0.025)
       
        hold on
        % add UW patient's data if dataset == 4
        if which_dataset == 4
            scatter(1-0.15, data(1,:), 'MarkerFaceColor', 'r', 'MarkerEdgeColor', 'r')
        end

        boxplot(data, 'Symbol','', 'jitter', 0.25, 'Color', col, 'widths', 0.1);
        title(figTitles{iFig})
        xticklabels(sampleNames{which_dataset+1})
        ylabel('Score')
    end
end

%% ========================================================================
% Analyse JPE task (= HF/Haruno&Frith) task 
%==========================================================================

% get items
HF_data = L.SVO_DA_3(:,22:57);

% online data only completed every 3rd item
if which_dataset == 0 || which_dataset == 5  
    HF_trials = 3:3:36; 
else
    HF_trials = 1:36;   
end

%%  HF models: get necessary variables from items

% absolute inequality 
HF(:,3) = abs(HF(:,1) - HF(:,2));

% joint gain 
HF(:,4) = sum(HF(:,1:2), 2); 

% for Fehr-Schmidt-Model = 
% alpha*self + beta*disadvantageous inequal. + gamma*advantageous inequal.
% disadvantageous inequality = max(other - self, 0)
HF(:,5) = HF(:,2) - HF(:,1);
% set negative values to zero (corresponds to taking max)
HF(HF(:,5)<0,5) = 0; 
% advantageous inequality = max(self - other, 0)
HF(:,6) = HF(:,1) - HF(:,2);
HF(HF(:,6)<0,6) = 0; 

% for ERC model (Equity, Reciprocity, Competition) = 
% alpha * self + beta * 0.5 * (self/(other+self) - 0.5 )^2
HF(:,7) = 0.5 * ( HF(:,1) ./ HF(:,4) - 0.5 ).^2;

% intercept 
HF(:,8) = 1; 

%% HF models: regression models 

% model names 
mod_vec = {'Individualism',...
    'Prosocial', ...
    'Altruism', ...
    'Joint-gain',...
    'Inequality',...
    'Inequality & Joint-gain',...
    'Fehr-Schmidt', ...
    'ERC'};

% model parameter names 
mod_xla = { '1=self', ...
    '1=self; 2=other', ...
    '1=other', ...
    '1=(self+other)',...
    '1=abs(self-other)', ...
    '1=abs(self-other); 2=(self+other)',...
    '1=self; 2=disadvantageous inequality; 3=advantageous inequality',...
    '1=self; 2=quadratic normalized inequality'};

X_mo{1} = HF(:,[8,1]);        % Individualism
X_mo{2} = HF(:,[8,1,2]);      % Altruism
X_mo{3} = HF(:,[8,2]);        % Simple altruism 
X_mo{4} = HF(:,[8,4]);        % joint gain model
X_mo{5} = HF(:,[8,3]);        % Inequality model 
X_mo{6} = HF(:,[8,3,4]);      % IA + JG
X_mo{7} = HF(:,[8,1,5,6]);    % Fehr-Schmidt
X_mo{8} = HF(:,[8,1,7]);      % ERC

% N of models
n_mod = size(X_mo,2);
% N of model parameter for each model 
k_mod = nan(n_mod,1);
for i_mod = 1:n_mod
    k_mod(i_mod,1) = size(X_mo{i_mod},2);
end

%% Regressions

% initialize 
n_obs = nan(n_sub, 1);
HF_y_var = nan(n_sub, 1); 
HF_B = cell(1, n_mod);
HF_R = cell(1, n_mod);
HF_STATS = cell(1, n_mod);
HF_SSR = nan(n_sub, n_mod);
HF_BIC = nan(n_sub, n_mod); 

for i_sub = 1:n_sub
    
    % response variable 
    HF_y = HF_data(i_sub, HF_trials)';
    
    % number of observations per subject 
    n_obs(i_sub,1) = size(HF_y, 1);
    % variance of observations per subject
    HF_y_var(i_sub,1) = var(HF_y);
    
    for i_mod = 1:n_mod
        
        % regress 
        [HF_B{i_mod}(i_sub,:),~,HF_R{i_mod}(i_sub,:),~,HF_STATS{i_mod}(i_sub,:)] = regress(HF_y, X_mo{i_mod}(HF_trials,:));
        
        % residual sum of squares (RSS) aka sum of squared residuals (SSR)
        HF_SSR(i_sub,i_mod) = sum(HF_R{i_mod}(i_sub,:) .^2);
        % Bayesian Information Criterion (gaussian special case) 
        HF_BIC(i_sub,i_mod) = n_obs(i_sub,1) * log(HF_SSR(i_sub,i_mod)/n_obs(i_sub,1)) + k_mod(i_mod,1) * log(n_obs(i_sub,1)); 
        
    end

end

%% exclude subject with zero variance in JPE(=HF) task 

HF_exclude = HF_y_var ~= 0; 

%% Plot Log Group Bays Factors (barplot)

if which_dataset == 3
    HF_BIC_sum = HF_BIC; % only one value if patient's data 
else
    HF_BIC_sum = sum(HF_BIC(HF_exclude,:)); % all other samples 
end 

% log-group bayes factors: LBF = 0.5 * (BIC_1 - BIC_2)
% substract 1st model 
LBF = 0.5 * (HF_BIC_sum - (HF_BIC_sum(1,1)));

% plot LBFs 
figure
barh(fliplr(LBF))
xlabel('Log-group Bayes factors')
set(gca,'yticklabel',fliplr(mod_vec),'FontSize',10.5)
title(sampleNames{which_dataset+1})

%% Display group mean svo scores and parameter estimates 

disp_data_tables; 
