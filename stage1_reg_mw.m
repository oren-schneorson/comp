function [stage1_] = stage1_reg_mw(T_, params)
%{
This function runs Fama-Macbeth first stage regression for
a single firm. 

All relevant data is contained i T_ (table) and params (dict).
params contains such relevant information as the
factors to be used in the regression and the returns variable.

The output collects regression results in a table, stage1_.

_mw variant of the function does this for each date in T_, using
expanding window with exponentially decreasing weights, with
half life equalling 5 years.

%}
%%

i = params('gvkey_i');
ngvkeys = params('ngvkeys');

if strcmp(params('freq'), 'D')
    freq_ = 252;
elseif strcmp(params('freq'), 'M')
    freq_ = 12;
else
    error('freq must be D or M.')
end


hh = freq_*5; % num busdays a year times 5 years
h = log(2)/hh;

weight_func = @(t_, tau_) exp(-abs(t_-tau_).*h)./sum(...
    exp(-abs(t_-[1:t_-1]).*h));


T_ = sortrows(T_, 'datadate', 'ascend'); % sorted by date

ngrps = 1; % single entity regression
factors = params('factors');
holidays = params('holidays');

% Fama MacBeth: stage 1
lead = 0; % lead X on Y, 0 by default.

% variable names for table of results of stage1
VariableNames = ...
    {...
    'min_date', 'max_date', 'nobs', 'Exret',...
    'alpha', 't_alpha',...
    'R_sq', 'F_stat', 'pval', 'sigma_e'};

% add beta_%d for each factor
VariableNames = [VariableNames,...
    arrayfun(@(f) sprintf('beta_%d', f),...
    1:numel(factors), 'UniformOutput', false)];
VariableNames = [VariableNames,...
    arrayfun(@(f) sprintf('t_beta_%d', f),...
    1:numel(factors), 'UniformOutput', false)];
VariableNames = [VariableNames,...
    arrayfun(@(f) sprintf('sigma_beta_%d', f),...
    1:numel(factors), 'UniformOutput', false)];



% sample 
idx_sample = true(size(T_, 1), 1);
idx_sample = idx_sample & T_.datadate <= params('max_date');
idx_sample = idx_sample & T_.datadate >= params('min_date');

if strcmp(params('freq'), 'D')
    idx_sample = idx_sample & abs(T_{:, params('returns_var')}) < 20;
elseif strcmp(params('freq'), 'M')
    idx_sample = idx_sample & abs(T_{:, params('returns_var')}) < 240;
end

idx_sample = idx_sample & idx_sample;


%idx = idx & ~(T.vol == 0 & T.ret == 0);
T_ = T_(idx_sample, :); % portfolio returns + date

clear idx_sample

% define regression variables
Y = T_{:, params('returns_var')}; % dependent variable
X = T_{:, params('factors')}; % independent variable
X = [NaN(lead, size(X,2)); X(1:end-lead,:)]; % lead/lag
X = [ones(size(X, 1), 1), X]; % add constant

aux = NaN(size(X,1), 1);
aux_ = repmat({aux}, 1, numel(VariableNames));
aux_{1} = NaT(size(X,1), 1);
aux_{2} = NaT(size(X,1), 1);

% allocate memory for table of stage1 results
stage1_ = table(aux_{:},...
    'VariableNames', VariableNames);

    for t = 1:size(X,1)
        
        clc
        fprintf('Computing conditional factor model, i: %.2f%%, t: %.2f%%',...
            i/ngvkeys*100, t/size(X,1)*100)
        
        idx = all(~isnan([X, Y]), 2); % valid rows
        idx = idx & T_.datadate <= T_.datadate(t);
        %idx = idx & (T_.ret ~= 0 | T_.VOL ~= 0);
        
    if sum(idx) == 0
        stage1_.min_date(t) = NaT;
        stage1_.max_date(t) = NaT;
        continue
    else
        min_date_it = min(T_.datadate(idx)); % min_date of valid rows
        max_date_it = max(T_.datadate(idx)); % max_date of valid rows
    end
    
    
    if ismember(max_date_it, stage1_.max_date)
        % valid rows causes the last rows to repeat regression
        continue
    end

    stage1_.min_date(t) = min_date_it;
    stage1_.max_date(t) = max_date_it;
    
    % overall mean excess return
    if strcmp(params('freq'), 'D')
        stage1_.Exret(t) = ...
            mean(T_.xretd(idx), 'omitnan');
    elseif strcmp(params('freq'), 'M')
        stage1_.Exret(t) = ...
            mean(T_.xretm(idx), 'omitnan');
    end
    
    %stage1_.mcap(t) = mean(T_.mcap(idx), 'omitnan');

    % stage 1 regression
    N = sum(idx); % # of obs
    K = size(X, 2); % # of regressors
    stage1_.nobs(t) = N;

    % skip dates with less than 100 valid obs (daily) or 13 (monthly)
    if strcmp(params('freq'), 'D') && N < 100
        continue
    elseif strcmp(params('freq'), 'M') && N < 24
        continue
    end


    % firms without any trades
    if all(Y == 0)
        continue
    end
    
    W = weight_func(t, 1:1:t)'; % weights decline over time
    %W = weight_func(t, t-2:-1:0)';
    %sum(W) ~= 1
    W = W(idx);
    W = W./sum(W, 'omitnan');
    %W_ = ones(sum(idx), 1); % equal weights
    %[W, W_]
    %return
    D = X(idx, :)'*diag(W)*X(idx, :); % inverse information matrix
    %[B, STDB, MSE, ~] = lscov(X(idx,:),Y(idx), W);    
    [B, ~, ~, ~] = lscov(X(idx,:),Y(idx), W);    
    % I think it doesn't make sense to use heterosk. adjusted with wls
    %[EstCoeffCov,se,coeff] = hac(X(idx,:),Y(idx), 'intercept', false);

    Y_hat = X*B; % expected Y
    R = Y(idx) - Y_hat(idx); % residuals
    SSE = sum(R.^2); % sum of sq resid
    RSS = sum((Y_hat(idx)-mean(Y(idx))).^2);
    TSS = sum((Y(idx)-mean(Y(idx))).^2); 
    
    sigma_hat = SSE/(N-K);
    Sigma_hat = inv(D)*sigma_hat;
    Sigma_hat = sqrt(diag(Sigma_hat));
    
    %tval = tinv((1-0.05/2), N-K);
    %BINT = [B-tval*sqrt(diag(Sigma_hat)), B+tval*sqrt(diag(Sigma_hat))];

    STATS(1) = 1-SSE/TSS; % R-squared
    STATS(2) = (RSS/(K-1))/sigma_hat; % F-stat
    STATS(3) = 1-fcdf(STATS(2), K-1, N-K); % p-val of F test
    STATS(4) = sigma_hat; % Cov of B
    
    % save regression results
    stage1_.R_sq(t) = STATS(1);
    stage1_.F_stat(t) = STATS(2);
    stage1_.pval(t) = STATS(3);
    stage1_.sigma_e(t) = STATS(4);

    stage1_.alpha(t) = B(1);
    %stage1_.alpha_int1(t) = BINT(1, 1);
    %stage1_.alpha_int2(t) = BINT(1, 2);
    stage1_.t_alpha(t) = B(1)/Sigma_hat(1); % regular OLS t-stat
    
    for f = 1:numel(factors)
        eval(sprintf('stage1_.beta_%d(t) = B(%d);', f, f+1))
        eval(sprintf('stage1_.sigma_beta_%d(t) = Sigma_hat(%d);', f, f+1))
        eval(sprintf('stage1_.t_beta_%d(t) = B(%d)/Sigma_hat(%d);', f, f+1, f+1))
    end 

    end

   stage1_ = stage1_(~isnan(stage1_.R_sq), :);
%%
end
