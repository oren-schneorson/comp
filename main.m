clear; clc
%{

main.m, part of a package for comp.

1. main.m (script)
2. preamble.m (script)
3. load_data.m (script)
4. merge2T.m (script)
5. get_indices.m (function)

This script loads daily/monthly U.S. stock market data. 

Quality checks: 
%}

%{
***************************
***** SETUP
***************************
%}


% set root
root = mfilename('fullpath');
if contains(root, 'LiveEditor') || numel(root) == 0
    % in case you're running via LiveEditor, default to root
    % change this line if running on a different system
    root = '/home/oren/Documents/MATLAB/comp';
    cd(root)
    
else
    root = fileparts(root);
end

% flags: currently working with only with daily
load_fpath_equity = false; % load the data from .mat file
monthly = false; % work with monthly frequency
dates_CBS = false;
monthly_avg = false; % end of month values or monthly average
glob = '';



if ~monthly && ~dates_CBS
    freq = 'D';
else
    freq = 'M';
end


% for Fama-Macbeth regressions, use returns_var 
if strcmp(freq, 'D')
    returns_var = {'xretd'};
elseif strcmp(freq, 'M')
    returns_var = {'xretm'};
else
    error('freq must be D or M.')
end


%{
Note: 3 Jan 2023

** daily + ~dates_CBS: V
    plain vanilla daily returns
** monthly + ~dates_CBS: monthly_avg: V
    computes monthly returns from first day of the month to last day
    of the month
** monthly + dates_CBS:
    same as with ~dates_CBS only using CBS dates instead of the last day of
    month
** daily + dates_CBS:
    really this is monthly, only using the return of days right after
    CBS inflation announcement.



%}




% Sample period
if strcmp(freq, 'D')
    period = 7;
elseif strcmp(freq, 'M')
    period = 5;
else
    error('freq must be D or M.')
end

% flag to override period in downstream scripts: false by default
override_period = false;

if period == 1
    min_date = datetime(2000,1,1);
    max_date = datetime(2007,12,31);
elseif period == 2
    min_date = datetime(2011,1,1);
    max_date = datetime(2019,12,31);
elseif period == 3
    min_date = datetime(2000,1,1);
    max_date = datetime(2019,12,31);
elseif period == 4
    min_date = datetime(2000,1,1);
    max_date = datetime(2022,12,31);
elseif period == 5 % monthly
    min_date = datetime(1951,12,31);
    max_date = datetime(2022,9,16);
elseif period == 6 % % daily
    min_date = datetime(1983,12,30);
    max_date = datetime(2022,9,16);
elseif period == 7 % % daily & TIPS
    min_date = datetime(2003,1,2);
    max_date = datetime(2022,9,16);
end

% preamble
fpath_m = fullfile(root, 'preamble.m');
run(fpath_m)

params = containers.Map();
vars = {...
    'min_date',...
    'max_date',...
    'holidays',...
    'returns_var',...
    'freq',...
    'monthly',...
    'dates_CBS'};

for v = vars
    v = char(v);
    eval(sprintf("params('%s') = %s;", v, v))
end



if exist(fpath_equity_full, 'file') && load_fpath_equity
    load(fpath_equity_full)
    % adj: for adusting weight of returns by public ownership
    %[~, ret_idx_ew] = get_indices(T, 'adj');

else

% load data
fpath_m = fullfile(root, 'load_data.m');
run(fpath_m)


% merge data to iids table
%fpath_m = fullfile(root, 'merge2T.m');
%run(fpath_m)

%{
% drop 2% of most illiquid firms
fname = 'liquidity_busdays.xlsx';
fpath = fullfile(lib_data, fname);
liq = readtable(fpath);
%idx = liq{:, 'busdays_median'} <= prctile(liq{:, 'busdays_median'}, 98);
%idx = ismember(T.iid, liq.iid(idx));
%T = T(idx, :);
%}

clc
% adj: for adusting weight by public ownership
%{
[~,...
    ret_idx_ew,...
    ret_idx,...
    xret_idx_ew,...
    xret_idx,...
    mcap_idx] = get_indices(T, 'reg', params);

save(fpath_equity_full, '-v7.3');
%}

end


%
clc

factors_tab = containers.Map();
factors_tab('epsilon_101') = 'full_arima';
factors_tab('epsilon_011') = 'full_arima';
factors_tab('innov') = 'innov';
factors_tab('Pi') = 'innov';
factors_tab('Pi12') = 'innov';
factors_tab('EPi12') = 'innov';
factors_tab('dEPi12') = 'innov';
factors_tab('e1') = 'innov';
factors_tab('e12') = 'innov';
factors_tab('ret_m') = 'mkt';
factors_tab('xret_m') = 'mkt';
factors_tab('retr_m') = 'mkt';

vars = {'ret_idx', 'ret_idx_ew', 'xret_idx', 'xret_idx_ew'};
for v = vars
    v = char(v);
    factors_tab(v) = v;
end

vars = {'HML_Epi_ret_idx', 'HML_Epi_ret_idx_ew'};
for v = vars
    v = char(v);
    factors_tab(v) = 'HML_Epi';
end

vars = {'HML_ret_idx', 'HML_ret_idx_ew'};
for v = vars
    v = char(v);
    factors_tab(v) = 'HML';
end

vars = {'SMB_ret_idx', 'SMB_ret_idx_ew'};
for v = vars
    v = char(v);
    factors_tab(v) = 'SMB';
end

vars = ex_PI_TIPS.Properties.VariableNames;
for v = vars
    v = char(v);
    factors_tab(v) = 'ex_PI_TIPS';
end


    
vars = ex_PI_Cleveland.Properties.VariableNames;
for v = vars
    v = char(v);
    factors_tab(v) = 'ex_PI_Cleveland';
end









return

%%

clc
%{

daily: dT5YIE
daily: dT10YIE
daily: dT5YIFR


%}
%ex_PI_D.Properties.VariableNames'

if strcmp(freq, 'D')
    returns_var = {'xretd'};
elseif strcmp(freq, 'M')
    returns_var = {'xretm'};
else
    error('freq must be D or M.')
end
if strcmp(freq, 'D')
    factors = {'dT5YIE'};
elseif strcmp(freq, 'M')
    factors = {'epsilon_101'};
    factors = {'UITB1'};
else
    error('freq must be D or M.')
end




params('factors') = factors;
params('returns_var') = returns_var;


% Fama-MacBeth regressions: stage 1 and 2
%mw = ''; % estimate betas in sample
mw = '_mw'; % % estimate betas OOS

load_fpath_stage1 = false; % if false, recreate (and overwrite) stage 1
load_fpath_stage2 = false; % if false, recreate (and overwrite) stage 2


% file name for excel file of stage%d
added_str_ = added_str('', factors,...
    min_date,...
    max_date,...
    freq,...
    monthly_avg,...
    mw,...
    '',...
    dates_CBS);
fname_stage = [added_str_, '.xlsx'];

if numel(mw) > 0
    fname_stage = strrep(fname_stage, '.xlsx', '.mat');
end

fpath_stage1 = fullfile(root, 'stage1', fname_stage);
fpath_stage2 = fullfile(root, 'stage2', fname_stage);


% ********************************************
% time-series regression (for each iid i)
% ********************************************
if load_fpath_stage1 && exist(fpath_stage1, 'file') == 2

    if numel(mw) > 0
        load(fpath_stage1)
    else
        stage1 = readtable(fpath_stage1);
    end
    
else
    % Fama-Macbeth, stage 1
    fpath_m = fullfile(root, ['stage1_regs', mw, '.m']);
    run(fpath_m)
end
%%
% join beta_Epi
idx = contains(stage1.Properties.VariableNames, 'beta');
stage1.Properties.VariableNames(idx) = cellfun(@(c)...
    strrep(c, 'beta_1', 'beta_Epi'),...
    stage1.Properties.VariableNames(idx),...
    'UniformOutput', false);
%stage1 = renamevars(stage1, 'id', 'iid');

idx_T = ~contains(T.Properties.VariableNames, 'beta_Epi');
T = outerjoin(T(:, idx_T), stage1(:, {'gvkey', 'beta_Epi'}),...
    'Type', 'left', 'Keys', 'gvkey', 'MergeKeys', true);


fprintf('Finished running stage1_regs%s. factors: %s\n',...
    mw, char(strjoin(factors, ', ')))
return
% ******************************************
% cross-section regression (for each date t)
% ******************************************
if load_fpath_stage2 && exist(fpath_stage2, 'file') == 2
    if numel(mw) > 0
        load(fpath_stage2)
    else
        stage2 = readtable(fpath_stage2);
    end

else
    % Fama-Macbeth, stage 2
    fpath_m = fullfile(root, ['stage2_regs', mw, '.m']);
    run(fpath_m)
end

%%
%{
monthly, end of month values, there are outlier dates
Issue was ret was not set to be < 20
A further check on the iid shows the usual prcbd=1 and then shoots up
%}
clc
idx = stage2.lambda_1 > prctile(stage2.lambda_1, 99);

% e.g. 31/08/2015 has an extreme risk premium, check stage 2
stage2(find(idx, 1), :)
return


%dates = stage2.date(find(idx, 1));
T(ismember(T.datadate, dates), :)


%%

% Chart 1: for each factor f,
% upper panel - excess return against beta
% lower panel - histogram of betas
clf
%subplot(1,2,1)
fpath_m = fullfile(root, 'chart1.m');
run(fpath_m)

%%

% Chart 1b: for each factor f,
% upper panel - excess return against beta
% lower panel - histogram of betas
clc
clf
%subplot(1,2,1)
for i_f = 1:numel(factors)
    
chart1a(stage1, i_f)
set(gcf, 'Color', 'w')
set(gca, 'FontSize', 16)

fname = strrep(fname_stage, '.xlsx', '.png');
fname = strrep(fname, '.mat', '.png');
fname = [sprintf('chart1a_%s', factors_nm{i_f}), fname];


fpath = fullfile(root, 'gfx', fname);
export_fig(fpath)

end

return
clf
%subplot(1,2,2)
chart1b(stage1, 1)

set(gcf, 'Color', 'w')
set(gca, 'FontSize', 16)

fname = strrep(fname_stage, '.xlsx', '.png');
fname = strrep(fname, '.mat', '.png');
fname = ['chart1b_', fname];


fpath = fullfile(root, 'gfx', fname);
export_fig(fpath)

%fpath_m = fullfile(root, 'chart1b.m');
%run(fpath_m)

%%


% Chart 2: for` each factor f,
% upper panel - price of risk factor over time + 12-months conditional mean
% lower panel - comparing factors' 12-months conditional mean price
%subplot(1,2,2)
clf
fpath_m = fullfile(root, 'chart2.m');
run(fpath_m)

%% Chart 3: plot of headline CPI

fpath_m = fullfile(root, 'chart3.m');
run(fpath_m)


%% Chart 4: market cap & size distribution

fpath_m = fullfile(root, 'chart4.m');
run(fpath_m)


%% chart number? plots individual firm daily and total return
clc
InputVariables = {'adj_prccd'};
GroupingVariables = {'iid'};

func = @(x) (x(end)/x(1));
aux = varfun(@(x) func(x), T,...
    'InputVariables', InputVariables,...
    'GroupingVariables', GroupingVariables);
aux.Properties.VariableNames(end) = {'ret'};

func = @(x) days252bus(x(1), x(end));
InputVariables = {'datadate'};
aux_ = varfun(@(x) func(x), T,...
    'InputVariables', InputVariables,...
    'GroupingVariables', GroupingVariables);
aux_.Properties.VariableNames(end) = {'days252bus'};
aux = [aux, aux_(:,end)];
clear aux_

aux.ret = aux.ret.^(1./aux.days252bus);
aux.ret = aux.ret.^252;
aux.ret = (aux.ret-1)*100;

for iid = iids'
clf
clc
fprintf('%.3f%%', find(iid==iids)/niids*100)


conm = securitiesinfo.conm(securitiesinfo.iid==iid);
idx = T.iid == iid;

subplot(1,2,1)
plot(T.datadate(idx), T.adj_prccd(idx))
title('Adjusted closing price')

subplot(1,2,2)
plot(T.datadate(idx), T.ret(idx))
title(sprintf('Returns, daily volatility: %.3f', nanstd(T.ret(idx))))

tit = [conm,...
    sprintf('Annual increase: %.3f%%', aux.ret(aux.iid==iid)),...
    sprintf('Mean daily return: %.3f%%', mean(T.ret(idx), 'omitnan')),...
    sprintf('Mean daily return (annualized): %.3f%%', ((mean(T.ret(idx), 'omitnan')/100+1)^252-1)*100)]
sgtitle(tit)


fname = [num2str(iid), '.png'];
fpath = fullfile(root, 'gfx', fname);
return
export_fig(fpath)


end





%% ARMA(1,1) correlation with dex_PI_M
clc
%head(ex_PI_M)
%head(full_arima)
tab = full_arima(:, [1, end]);
tab.datadate = tab.datadate+calmonths(0); % shifting full_arima by one month
tab.datadate = eomdate(tab.datadate);
tab = innerjoin(ex_PI_Cleveland(:, [1, (size(ex_PI_Cleveland, 2)-1)/2+2:end]), tab);
a=corr(tab{:, 2:end}, 'rows', 'pairwise');
sum(a,2)/size(a,1)
[a(end,:)', a(:,end)]


head(tab)
%help corr
% correlation of about 25-37% between epsilon and Israeli TIPS (monthly)
