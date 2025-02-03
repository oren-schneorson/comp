
%{
This script shows how the code of get_factor is used.
Two main scripts are based on this:
1. compute_factors.m -- constructs and saves factor models
(equal and value weighted portfolios on sorts)
2. Fama-French replication, the cross section of stock returns.

%}

% TODO: b2m and sz are not exactly the same sample.
% all companies have mcap, not all companies have b2m
% TODO: monthly, get_factor, dates_CBS...



%%
clc
InputVariables = {'beta_Epi'};
model = 'inflation';
nptiles = 10;
bp_Epi = {[1:nptiles-1]*100/nptiles};

PycharmProjects_dir = fullfile('/home', username, 'PycharmProjects');



[ret_idx_ew_Epi, ret_idx_Epi, xret_idx_ew_Epi, xret_idx_Epi,...
    GroupCount_Epi,...
    prctiles_Epi, ~,...
    mcaps_Epi, T_InputVariables_Epi, HML_Epi] = ...
    get_factor(T, InputVariables, bp_Epi, params);


% Boons works with HML not LMH
%HML_Epi{:, 2:end} = -HML_Epi{:, 2:end}; % high inflation, low return

HML_Epi.Properties.VariableNames(2:end) = strcat('HML_Epi_',...
    HML_Epi.Properties.VariableNames(2:end));


clc
added_str_ = added_str(model, factors, min_date, max_date, freq, monthly_avg, mw, '', dates_CBS);
fname = [added_str_, '.mat'];
fname = strrep(fname, model, [model, '_ptiles']);
fpath = fullfile(root,'generated_data', fname);


save(fpath,...
    'ret_idx_ew_Epi',...
    'ret_idx_Epi',...
    'xret_idx_ew_Epi',...
    'xret_idx_Epi',...
    'GroupCount_Epi',...
    'prctiles_Epi',...
    'mcaps_Epi',...
    'bp_Epi',...
    'T_InputVariables_Epi', 'HML_Epi')


%%
ma_days = 365;
clc
clf
subplot(2,2,1)
hold on


idx_sample = true(size(ret_idx_ew_Epi.datadate));
%idx_sample = idx_sample & ret_idx_ew_Epi.datadate > datetime(1999,12,31);
%idx_sample = idx_sample & ret_idx_Epi.datadate <= datetime(2006,12,31);
idx_sample = idx_sample & ret_idx_ew_Epi.datadate > datetime(2009,12,31);
idx_sample = idx_sample & ret_idx_Epi.datadate <= datetime(2019,12,31);
plot(mean(ret_idx_ew_Epi{idx_sample, 2:end}, 'omitnan'))
plot(mean(ret_idx_Epi{idx_sample, 2:end}, 'omitnan'))
legend({'equally-weighted', 'value-weighted', })
title('IRP, mean return')



%{

The result I get is the following:
ARMA(1,1)


Posit IRP for the first decade of 2000's, thence about zero.

A negative IRP: stocks that covaried positively with inflation earned
higher returns than stocks that covaried negatively with inflation.

This is the opposite of what the literature finds in the U.S., e.g.
Boons et al. (2020) and Ang et al. (2012). Specifically, Boons
found that when the relationship between the whole stock
market and 

This is the result whether I use:
(frequency) daily data or monthly data
(measure)   inflation innovtions from capital markets, actual inflation
            or inflation innovations coming from ARMA(1,1)

It is important to note that this is for the In-Sample exercise.
The out-of-sample exercise?

Note: Boons has a

%}
size(ret_idx_ew_Epi.datadate)
head(ret_idx_ew_Epi)
size(movmean(ret_idx_ew_Epi{:, 2:end}, ma_days))
subplot(2,2,2)
plot(ret_idx_ew_Epi.datadate, movmean(ret_idx_ew_Epi{:, 2:end}, ma_days))
%
leg = ret_idx_ew_Epi.Properties.VariableNames(2:end);
leg = cellfun(@(c) strrep(c, '_', '\_'), leg, 'UniformOutput', false);
legend(leg)
title('return index')

subplot(2,2,3)
plot(GroupCount_Epi.date, GroupCount_Epi{:, 2:end})
title('GroupCount')
return
plot(ret_idx_ew_Epi.date, movmean(...
    ret_idx_ew_Epi{:, 2}-ret_idx_ew_Epi{:, end},...
    ma_days))
%}
%plot(HML_Epi.date+ma_days/2, movmean(HML_Epi.HML_Epi_ret_idx, ma_days))



%%


clc
InputVariables = {'mcap', 'beta_Epi'};
model = 'sz_inflation';
nptiles = 5;
%bp = {[1:nptiles-1]*100/nptiles};
bp_sz_Epi = {[20, 50], [1:nptiles-1]*100/nptiles};
[...
    ret_idx_ew_sz_Epi, ret_idx_sz_Epi,...
    xret_idx_ew_sz_Epi, xret_idx_sz_Epi,...
    GroupCount_sz_Epi,...
    prctiles_sz_Epi, ~,...
    mcaps_sz_Epi, T_InputVariables_sz_Epi, ~] =...
    get_factor_cond(T, InputVariables, bp_sz_Epi, params);



added_str_ = added_str(model, [InputVariables(1), factors], min_date, max_date, freq, monthly_avg, '', '', dates_CBS);
fname = [added_str_, '.mat'];
fname = strrep(fname, model, [model, '_ptiles']);
fpath = fullfile(root,'generated_data', fname);
save(fpath,...
    'ret_idx_ew_sz_Epi',...
    'ret_idx_sz_Epi',...
    'xret_idx_ew_sz_Epi',...
    'xret_idx_sz_Epi',...
    'GroupCount_sz_Epi',...
    'prctiles_sz_Epi',...
    'mcaps_sz_Epi',...
    'bp_sz_Epi',...
    'T_InputVariables_sz_Epi')


%{
2023-06-28

Note: data loss.
I'm losing many observations, but this is necessary for the procedure.
First, some companies do not have beta_Epi at all. Second, The first month
of observations is deleted for all companies. All in all it's about 1.41%
of the data, a figure I can live with.
    
The firms that do not have any beta_Epi are omitted in the beginning of the
procedure, so I'm not reall seeing this in T_InputVariables...
    
On the whole 0.48% can be attributed to first month of firm. On average,
this means firms have about 208 months of data.
    
Tough I have much more, 162 (168) on average (median)...
I checked that I'm not losing data at the end of firm life.
    
I've found an extra reason to lose observations.
The long tail of any group is omitted.
       
%}
    
%{
% Quality check: c should be empty....
a=varfun(@min, T_InputVariables_sz_Epi(T_InputVariables_sz_Epi.grp==0,:), 'InputVariables', 'date', 'GroupingVariables', 'iid');
b=varfun(@max, T_InputVariables_sz_Epi(T_InputVariables_sz_Epi.grp==0,:), 'InputVariables', 'date', 'GroupingVariables', 'iid');
c=join(a,b);
d=between(c.min_date,c.max_date);
c.d=calmonths(d);
c(c.d>0,:)
    
%}
%% average mcap of inflation portfolios
clc

factor_nm = 'sz_Epi';
returns_var = 'ret';
w = '_ew';

eval(sprintf("nptiles = cellfun(@numel, bp_%s)+1;;", factor_nm));
eval(sprintf("a = innerjoin(mcaps_%s, GroupCount_%s);", factor_nm, factor_nm));



b = a(:, 1:prod(nptiles)+1);

b{:, 2:prod(nptiles)+1} =...
    b{:, 2:prod(nptiles)+1}./a{:, prod(nptiles)+2:end};

clf
nsubplots = ceil(sqrt(nptiles(1)));
for i = 1:nptiles(1)
subplot(nsubplots,nsubplots,i)
plot(b.date, b{:,1+i:nptiles(1):end})
end

subplot(nsubplots,nsubplots,nsubplots^2)
eval(sprintf("a = GroupCount_%s;", factor_nm));
plot(a.date, a{:,2:end})


%% inflation portfolios, IRP in mean return
clc
factor_nm = 'sz_Epi';
returns_var = 'ret';
w = '_ew';

eval(sprintf("nptiles = cellfun(@numel, bp_%s)+1;", factor_nm));
eval(sprintf("a = %s_idx%s_%s;", returns_var, w, factor_nm));
b = a(:, 'date');

for i = 1:nptiles(2)
    
    b = [b, array2table(mean(a{:,...
        2+(i-1)*nptiles(1):1+i*nptiles(1)},...
        2, 'omitnan'))];
    b = renamevars(b, 'Var1', sprintf("%s_idx%s_%d;", returns_var, w, i));
    
end

clf
plot(b.date, b{:, 2:end})

idx_sample = true(size(b, 1), 1);
idx_sample = idx_sample & b.date >= datetime(2010,1,1);
idx_sample = idx_sample & b.date <= datetime(2019,12,31);

subplot(1,2,1)
c = mean(b{idx_sample, 2:end}, 'omitnan')
plot(c)

subplot(1,2,2)
ma_days = 252*3;
%ma_days = 12;
d = movmean(b{:, 2:end}, ma_days, 'omitnan');

plot(b.date(idx_sample), d(idx_sample, :))


    
%%
clc
bp = flip(bp_sz_Epi);
bp = bp_sz_Epi;
cellfun(@numel, bp)+1
[ia, ib] = ind2sub(cellfun(@numel, bp)+1, [1:size(mcaps_sz_Epi, 2)-1]');

sub2ind(cellfun(@numel, bp)+1, ia ,ib)
return

[num2str(ia), num2str(ib)]

idx = find(ib==2)+1;
clf
plot(mcaps_sz_Epi.date, mcaps_sz_Epi{:, idx})




%%
clc
InputVariables = {'mcap', 'b2m'};
model = 'bm_sz';
nptiles = 5;
%bp = {[1:nptiles-1]*100/nptiles};
bp_bm_sz = {[20, 50], [1:nptiles-1]*100/nptiles};
[ret_idx_ew_bm_sz, ret_idx_bm_sz, xret_idx_ew_bm_sz, xret_idx_bm_sz,...
    GroupCount_bm_sz,...
    prctiles_bm_sz, ~,...
    mcaps_bm_sz, T_InputVariables_bm_sz, ~] = get_factor(T, InputVariables, bp_bm_sz, params);


added_str_ = added_str(model, InputVariables, min_date, max_date, freq, monthly_avg, '', '', dates_CBS);
fname = 'bm_sz_ptiles_';
fname = [fname, added_str_, '.mat'];
fpath = fullfile(root,'generated_data', fname);

save(fpath,...
    'ret_idx_ew_bm_sz',...
    'ret_idx_bm_sz',...
    'xret_idx_ew_bm_sz',...
    'xret_idx_bm_sz',...
    'GroupCount_bm_sz',...
    'prctiles_bm_sz',...
    'mcaps_bm_sz',...
    'bp_bm_sz',...
    'T_InputVariables_bm_sz')


InputVariables = {'b2m'};
model = 'bm';
nptiles = 10;
bp_bm = {[1:nptiles-1]*100/nptiles};

[ret_idx_ew_bm, ret_idx_bm, xret_idx_ew_bm, xret_idx_bm,...
    GroupCount_bm,...
    prctiles_bm, ~,...
    mcaps_bm, T_InputVariables_bm, HML] = get_factor(T, InputVariables, bp_bm, params);

HML.Properties.VariableNames(2:end) = strcat('HML_',...
    HML.Properties.VariableNames(2:end));


added_str_ = added_str(model, InputVariables, min_date, max_date, freq, monthly_avg, '', '', dates_CBS);
fname = 'bm_ptiles_';
fname = [fname, added_str_, '.mat'];
fpath = fullfile(root,'generated_data', fname);

save(fpath,...
    'ret_idx_ew_bm',...
    'ret_idx_bm',...
    'xret_idx_ew_bm',...
    'xret_idx_bm',...
    'GroupCount_bm',...
    'prctiles_bm',...
    'mcaps_bm',...
    'bp_bm',...
    'T_InputVariables_bm', 'HML')

InputVariables = {'mcap'};
model = 'size';

bp_sz = {[1:nptiles-1]*100/nptiles};

[ret_idx_ew_sz, ret_idx_sz, xret_idx_ew_sz, xret_idx_sz,...
    GroupCount_sz,...
    prctiles_sz, ~,...
    mcaps_sz, T_InputVariables_sz, SMB] = get_factor(T, InputVariables, bp_sz, params);

SMB{:, 2:end} = -SMB{:, 2:end}; % factor comes out as HML, revert to SMB

SMB.Properties.VariableNames(2:end) = strcat('SMB_',...
    SMB.Properties.VariableNames(2:end));

added_str_ = added_str(model, InputVariables, min_date, max_date, freq, monthly_avg, '', '', dates_CBS);
fname = 'sz_ptiles_';
fname = [fname, added_str_, '.mat'];
fpath = fullfile(root,'generated_data', fname);

save(fpath,...
    'ret_idx_ew_sz',...
    'ret_idx_sz',...
    'xret_idx_ew_sz',...
    'xret_idx_sz',...
    'GroupCount_sz',...
    'prctiles_sz',...
    'mcaps_sz',...
    'bp_sz',...
    'T_InputVariables_sz', 'SMB')


%%
clf
clc
idx = stage1.nobs>1000;
histogram(stage1.beta_Epi(idx))
sortrows(stage1(idx,:), 'nobs')

%%
clc

stage1(stage1.nobs>=252,:)
sum(stage1.nobs>=252)
return

head(sortrows(stage1, 'nobs'))
min(T.date)
min(stage1.min_date)


%%
clc
cmap = flip(bone);


ind = cellfun(@(c) str2double(strrep(c, 'ret_idx_ew_', '')),...
    ret_idx_ew_bm_sz.Properties.VariableNames(2:end))';
[ia, ib] = ind2sub(nptiles*ones(1, numel(factors_)), ind);

mean_ret_bm_sz = varfun(@(x) mean(x, 'omitnan'), ret_idx_ew_bm_sz,...
    'InputVariables', ret_idx_ew_bm_sz.Properties.VariableNames(2:end));
mean_ret_bm_sz = mean_ret_bm_sz{1,:};
mean_GroupCount = varfun(@(x) mean(x, 'omitnan'), GroupCount_bm_sz,...
    'InputVariables', GroupCount_bm_sz.Properties.VariableNames(2:end));
mean_GroupCount = mean_GroupCount{1,:};

mean_ret_bm = varfun(@(x) mean(x, 'omitnan'), ret_idx_ew_bm,...
    'InputVariables', ret_idx_ew_bm.Properties.VariableNames(2:end));
mean_ret_sz = varfun(@(x) mean(x, 'omitnan'), ret_idx_ew_sz,...
    'InputVariables', ret_idx_ew_sz.Properties.VariableNames(2:end));
mean_ret_al = mean(ret_idx_ew{:,end}, 'omitnan');

clf

subplot(1,2,1)

aux = reshape(mean_ret_bm_sz, nptiles*ones(1, numel(factors_)));
[~, h] = contourf(1:nptiles, 1:nptiles, aux);
colormap(cmap)
colorbar
h.LineStyle = 'none';

title('Avg daily return')
xlabel(sprintf('Size quintiles'))
ylabel(sprintf('Book to market quintiles'))

VariableNames = split(num2str(1:nptiles), '  ');
VariableNames = cellfun(@(c) ['BM-', c], VariableNames, 'UniformOutput', false);
RowNames = split(num2str(1:nptiles), '  ');
RowNames = cellfun(@(c) ['ME-', c], RowNames, 'UniformOutput', false);

VariableNames{1} = 'Low-BM';
VariableNames{end} = 'High-BM';
RowNames{1} = 'Small-ME';
RowNames{end} = 'Big-ME';



%{
% axis 1 flipped
aux = flip(aux);
RowNames = flip(RowNames);

%}

aux = [mean_ret_bm{:,:}; aux];
aux = [[mean_ret_al; mean_ret_sz{:,:}'], aux];

VariableNames = [{'All'}; VariableNames];
RowNames = [{'All'}; RowNames];

aux = mat2cell(aux, size(aux, 1), ones(1, size(aux, 2)));
tab = table(aux{:},...
    'VariableNames', VariableNames,...
    'RowNames', RowNames);
fpath = fullfile(PycharmProjects_dir, 'Stocks_inflation', 'data_tab4.csv');
writetable(tab, fpath, 'WriteRowNames', true)

subplot(1,2,2)
aux = reshape(mean_GroupCount, nptiles*ones(1, numel(factors_)));
[~, h] = contourf(1:nptiles, 1:nptiles, aux);

colormap(cmap)
colorbar
h.LineStyle = 'none';

title('Average group count')
xlabel(sprintf('Size quintiles'))
ylabel(sprintf('Book to market quintiles'))



%%
clc
factors_ = {'beta_1'};
nptiles = 10;
[ret_idx_ew_beta, GroupCount, prctiles, ret_idx_beta] =...
    get_factor(T, factors_, nptiles, params);

mean_GroupCount = varfun(@(x) mean(x, 'omitnan'), GroupCount,...
    'InputVariables', GroupCount.Properties.VariableNames(2:end));
mean_GroupCount = mean_GroupCount{1,:};



%
VariableNames = arrayfun(@(i) ['decile', num2str(i)], 1:nptiles, 'UniformOutput', false);
RowNames = {'alpha_post', 'beta_post', 't_stat_alpha', 't_stat_beta',...
    'N', 'R_sq', 'mean_GroupCount', 'mean_ret'};
nrows = numel(RowNames);

aux = repmat({NaN(nrows,1)}, 1, nptiles);
tab = table(aux{:}, 'VariableNames', VariableNames, 'RowNames', RowNames);

%%

%aux = innerjoin(ret_idx_ew_beta, ex_PI_D(:, [{'date'}, factors]));
aux = innerjoin(ret_idx_ew_beta, mkt(:, {'date', 'ret_m'}));
aux = innerjoin(aux, R_f);
aux.xret_m = aux.ret_m-aux.(R_f.Properties.VariableNames{end});


for f = 1:nptiles

X = aux{:,factors};
X = [ones(size(X)), X]; % include a constant to measure alpha
Y = aux{:,f+1};

idx = ~any(isnan([Y, X]), 2);

N = sum(idx);
K = size(X, 2);

D = X(idx, :)'*X(idx, :);

[B,BINT,R,~,STATS] = regress(Y(idx),X(idx,:));
[EstCoeffCov,se,coeff] = hac(X(idx,:),Y(idx), 'intercept', false);


sigma_hat = R'*R/(N-K);
Sigma_hat = inv(D)*sigma_hat;

[BINT(:,1), B, BINT(:,2)]
%t = B./sqrt(diag(Sigma_hat)); % regular OLS
t = B./se; % regular OLS
tab{1, f} = B(1); % alpha
tab{2, f} = B(2); % beta_ex_PI
tab{3, f} = t(1); % t-stat alpha
tab{4, f} = t(2); % t-stat beta
tab{5, f} = N; % # of obs
tab{6, f} = STATS(1); % R squared
tab{7, f} = mean_GroupCount(f); % mean group count
tab{8, f} = mean(Y(idx)); % mean return, portfolio


end




mean_ret_beta = varfun(@(x) mean(x, 'omitnan'), ret_idx_ew_beta,...
    'InputVariables', ret_idx_ew_beta.Properties.VariableNames(2:end));
mean_ret_beta = mean_ret_beta{:,:}';
cov_ret_beta = cov(ret_idx_ew_beta{2:end,2:end});
t = mean_ret_beta./sqrt(diag(cov_ret_beta));

%
mean_ret = mean_ret_beta';
mean_ret = mat2cell(mean_ret, size(mean_ret,1), ones(1, size(mean_ret, 2)));
mean_ret = table(mean_ret{:});

mean_ret.Properties.VariableNames = tab.Properties.VariableNames;

%{
tab = [tab; mean_ret];
tab.Properties.RowNames{end} = 'mean_ret';
%}

sharpe = (mean_ret_beta'-mean(aux{:,end}))./sqrt(diag(cov_ret_beta))';
sharpe = mat2cell(sharpe, size(sharpe,1), ones(1, size(sharpe, 2)));
sharpe = table(sharpe{:});
sharpe.Properties.VariableNames = tab.Properties.VariableNames;
tab = [tab; sharpe];

tab.Properties.RowNames{end} = 'Sharpe-ratio';


X = aux{:,end};
Y = aux{:,end-1}-aux{:,2};
X = [ones(size(X)), X];

idx = ~any(isnan([Y, X]), 2);

N = sum(idx);
K = size(X, 2);

D = X(idx, :)'*X(idx, :);


[B,BINT,R,~,STATS] = regress(Y(idx),X(idx,:));
[EstCoeffCov,se,coeff] = hac(X(idx,:),Y(idx), 'intercept', false);

sigma_hat = R'*R/(N-K);
Sigma_hat = inv(D)*sigma_hat;

[BINT(:,1), B, BINT(:,2)]

% t = B./sqrt(diag(Sigma_hat)); % regular OLS
t = B./se; % NW se
sharpe = (mean_ret_beta(end)-mean_ret_beta(1)-mean(aux{:,end}))./...
    sum(sum((cov_ret_beta([1, end], [1, end]))));


tab.HML = [B(1); B(2); t(1); t(2); N; STATS(1); mean_GroupCount(1)+mean_GroupCount(end); sharpe; mean_ret{1, end}-mean_ret{1, 1}];


fpath = fullfile(PycharmProjects_dir, 'Stocks_inflation', 'data_tab5.csv');
writetable(tab, fpath, 'WriteRowNames', true)



%% for interacting two factors, plot contourf of returns table
clc
cmap = flip(bone);


ind = cellfun(@(c) str2num(strrep(c, 'ret_idx_ew_', '')),...
    ret_idx_ew_bm_sz.Properties.VariableNames(2:end))';
[ia, ib] = ind2sub(nptiles*ones(1, numel(factors_)), ind);

mean_ret_beta = varfun(@(x) mean(x, 'omitnan'), ret_idx_ew_bm_sz,...
    'InputVariables', ret_idx_ew_bm_sz.Properties.VariableNames(2:end));
mean_GroupCount = varfun(@(x) mean(x, 'omitnan'), GroupCount,...
    'InputVariables', GroupCount.Properties.VariableNames(2:end));
mean_GroupCount = mean_GroupCount{1,:}
return

mean_ret_al = mean(ret_idx_ew{:,end}, 'omitnan');

clf

subplot(1,2,1)

aux = reshape(mean_ret_bm_sz, nptiles*ones(1, numel(factors_)));
[~, h] = contourf(1:nptiles, 1:nptiles, aux);
colormap(cmap)
colorbar
h.LineStyle = 'none';

title('Avg daily return')
xlabel(sprintf('Size quintiles'))
ylabel(sprintf('Book to market quintiles'))

VariableNames = split(num2str(1:nptiles), '  ');
VariableNames = cellfun(@(c) ['BM-', c], VariableNames, 'UniformOutput', false);
RowNames = split(num2str(1:nptiles), '  ');
RowNames = cellfun(@(c) ['ME-', c], RowNames, 'UniformOutput', false);

VariableNames{1} = 'Low-BM';
VariableNames{end} = 'High-BM';
RowNames{1} = 'Small-ME';
RowNames{end} = 'Big-ME';



%{
% axis 1 flipped
aux = flip(aux);
RowNames = flip(RowNames);

%}

aux = [mean_ret_bm{:,:}; aux];
aux = [[mean_ret_al; mean_ret_sz{:,:}'], aux];

VariableNames = [{'All'}; VariableNames];
RowNames = [{'All'}; RowNames];

aux = mat2cell(aux, size(aux, 1), ones(1, size(aux, 2)));
tab = table(aux{:},...
    'VariableNames', VariableNames,...
    'RowNames', RowNames);

fpath = fullfile(PycharmProjects_dir, 'Stocks_inflation', 'data_tab4.csv');
writetable(tab, fpath, 'WriteRowNames', true)

subplot(1,2,2)
aux = reshape(mean_GroupCount, nptiles*ones(1, numel(factors_)));
[~, h] = contourf(1:nptiles, 1:nptiles, aux);

colormap(cmap)
colorbar
h.LineStyle = 'none';

title('Average group count')
xlabel(sprintf('Size quintiles'))
ylabel(sprintf('Book to market quintiles'))

