function [stage1_] = stage1_reg(T_, params)
%{
This function runs Fama-Macbeth

%}

if strcmp(params('freq'), 'D')
    var_ret = 'retd';
    var_xret = 'xretd';
    var_retr = 'retrd';
elseif strcmp(params('freq'), 'M')
    var_ret = 'retm';
    var_xret = 'xretm';
    var_retr = 'retrm';
else
    error('freq must be D or M.')
end


T_ = sortrows(T_, 'datadate', 'ascend'); % sorted by date


ngrps = 1; % single entity regression
factors = params('factors');
holidays = params('holidays');

% Fama MacBeth: stage 1
lead = 0; % lead X on Y, 0 by default.


% variable names for table of results of stage1
VariableNames = ...
    {...
    'min_date', 'max_date', 'gvkey', 'mcap', 'nobs', 'Eret',...
    'alpha_int1', 'alpha', 'alpha_int2', 't_alpha',...
    'R_sq', 'F_stat', 'pval', 'sigma_e'};
% add beta_%d for each factor
VariableNames_ = arrayfun(@(f) {...
    sprintf('beta_%d_int1', f),...
    sprintf('beta_%d', f),...
    sprintf('beta_%d_int2', f),...
    sprintf('t_beta_%d', f)},...
    1:numel(factors), 'UniformOutput', false);
VariableNames = [VariableNames, VariableNames_{:}];


aux = NaN(ngrps, 1);
aux_ = repmat({aux}, 1, numel(VariableNames));
aux_{1} = NaT(ngrps, 1);
aux_{2} = NaT(ngrps, 1);

% allocate memory for table of stage1 results
stage1_ = table(aux_{:},...
    'VariableNames', VariableNames);


% sample 
idx_sample = true(size(T_, 1), 1);
idx_sample = idx_sample & T_.datadate <= params('max_date');
idx_sample = idx_sample & T_.datadate >= params('min_date');
idx_sample = idx_sample & T_.(var_ret) ~= 0;

if strcmp(params('freq'), 'D')
    idx_sample = idx_sample & abs(T_{:, var_ret}) < 20;
elseif strcmp(params('freq'), 'M')
    idx_sample = idx_sample & abs(T_{:, var_ret}) < 240;
else
    error('freq must be D or M.')
end

%idx_sample = idx_sample & idx_sample;


%idx = idx & ~(T.vol == 0 & T.ret == 0);
T_ = T_(idx_sample, :); % portfolio returns + date
stage1_.mcap = mean(T_.mcap, 'omitnan');
clear idx_sample

% define regression variables
Y = T_{:, params('returns_var')}; % dependent variable
X = T_{:, params('factors')}; % independent variable
X = [NaN(lead, size(X,2)); X(1:end-lead,:)]; % lead/lag
X = [ones(size(X, 1), 1), X]; % add constant

idx = all(~isnan([X, Y]), 2); % valid rows
idx = idx & all(isreal([X, Y]), 2); % valid rows
idx = idx & all(~isinf([X, Y]), 2); % valid rows
%idx = idx & T_.ret ~= 0; % should not be an issue after illiquid
%securities are dealt with

if sum(idx) == 0
    stage1_.min_date = NaT;
    stage1_.max_date = NaT;
    return
else
    min_date_i = min(T_.datadate(idx)); % min_date of valid rows
    max_date_i = max(T_.datadate(idx)); % max_date of valid rows
    
end


stage1_.min_date = min_date_i;
stage1_.max_date = max_date_i;

% overall mean (excess) return

%xr = T_{:, 'xret'}; % excess returns
%stage1.Eret(i) = mean(xr(idx), 'omitnan');

if strcmp(params('freq'), 'D')
    stage1_.Eret = mean(T_.(var_ret), 'omitnan');
elseif strcmp(params('freq'), 'M')
    stage1_.Eret = mean(T_.retm, 'omitnan');
end


% stage 1 regression
N = sum(idx); % # of obs
%K = size(X, 2); % # of regressors
stage1_.nobs = N;

% skip firms with less than 100 valid obs (daily) or 12 (monthly)
if strcmp(params('freq'), 'D') && N < 100
    return
elseif strcmp(params('freq'), 'M') && N < 12
    return
end


% firms without any trades
if all(Y == 0)
    return
end


[B,BINT,~,~,STATS] = regress(Y(idx), X(idx,:));
% debug
%all(~any(isnan([Y(idx), X(idx,:)]), 2))
%all(~any(isinf([Y(idx), X(idx,:)]), 2))
%[Y(idx), X(idx,:)]

%[B,BINT,R,~,STATS] = regress(Y(idx),X(idx,:));
%[EstCoeffCov, se, coeff] = hac(X(idx,:),Y(idx), 'intercept', false);
[~, se, ~] = hac(X(idx,:), Y(idx), 'intercept', false, 'display', 'off');

% save regression results
stage1_.R_sq = STATS(1);
stage1_.F_stat = STATS(2);
stage1_.pval = STATS(3);
stage1_.sigma_e = STATS(4);

stage1_.alpha = B(1);
stage1_.alpha_int1 = BINT(1, 1);
stage1_.alpha_int2 = BINT(1, 2);

t = B./se; % t-stats (Newey-West robust)


stage1_.t_alpha = B(1)/se(1); % Newey-West t-stat (alpha)
for f = 1:numel(factors)
    eval(sprintf('stage1_.beta_%d = B(%d);', f, f+1))
    eval(sprintf('stage1_.beta_%d_int1 = BINT(%d, 1);', f, f+1))
    eval(sprintf('stage1_.beta_%d_int2 = BINT(%d, 2);', f, f+1))
    eval(sprintf('stage1_.t_beta_%d = t(%d);', f, f+1))  % Newey-West t-stat (beta)
end


end