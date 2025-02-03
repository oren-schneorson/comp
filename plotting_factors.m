clear; clc
%{
Plotting factors I generated.
Plotting the IRP

TODO:
# compute quantity of risk of HML (new script)
# compute quantity of risk of MKT (new script)
# get correlation of HML and MKT

# run a two factor model with MKT and HML-inflation

%}

% set root
username = getenv("USER");
root = mfilename('fullpath');
if contains(root, 'LiveEditor') || numel(root) == 0
    % in case you're running via LiveEditor, default to root
    % change this line if running on a different system
    root = fullfile('/home', username, 'Documents/MATLAB/comp');
    cd(root)
    
else
    root = fileparts(root);
end

cd(root)


lib_matlab = fullfile('/home', username, 'Documents', 'MATLAB');
lib_data = fullfile('/media', username, 'D', 'data');
lib_pycharm = fullfile('/home', username, 'PycharmProjects');

addpath(fullfile(matlab_dir, 'my_functions'))
addpath(fullfile(matlab_dir, 'altmany-export_fig-410f0ad'))



%dims = {'mcap', 'stage1'}
dims = {'stage1'}
dims = {'mcap', 'stage1'};
%dims = {'mcap', };
freq = 'M';
months_lag = 1;

% dims2columns = {'mcap': 'mcap', 'stage1': 'beta_1'}
%dims2grp = ['grp_%s' % dim for dim in dims]
dims2grp = cellfun(@(dim) ['grp_', dim], dims, 'UniformOutput',false);
dims2columns = containers.Map({'mcap', 'stage1'}, {'mcap', 'beta_1_vasicek'});

joined_dims = char(join(dims, '_'));

lib_comp = fullfile(lib_data, 'comp');
lib_comp_pycharm = fullfile(lib_pycharm, 'comp');
lib_secd = fullfile(lib_comp, 'secd');
lib_secm = fullfile(lib_comp, 'secm');
lib_secd_proc = fullfile(lib_comp, 'secd_proc');
lib_secm_proc = fullfile(lib_comp, 'secm_proc');
lib_fundq = fullfile(lib_comp, 'fundq');
lib_matlab_comp = fullfile(lib_matlab, 'comp');
lib_stage1 = fullfile(lib_matlab_comp, 'stage1', 'mw', freq);
lib_grps = fullfile(lib_matlab_comp, 'stage1', 'mw', 'grps', joined_dims, freq);


VariableNames_0 = {'S', 'M', 'B'};
VariableNames_1 = strcat('grp', arrayfun(@num2str, 0:9, 'UniformOutput', false));

if all(cellfun(@(dim) strcmp(dim, 'mcap'), dims))
    VariableNames = VariableNames_0;
elseif all(cellfun(@(dim) strcmp(dim, 'stage1'), dims))
    VariableNames = VariableNames_1;

elseif all(ismember({'mcap', 'stage1'}, dims))

    
    VariableNames = {};
    for VariableName_0 = VariableNames_0
        VariableName = strcat([char(VariableName_0), '_'], VariableNames_1);
        VariableNames = [VariableNames, VariableName];
    end
    
else
    VariableNames = {};
end

VariableNames = [{'datadate'}, VariableNames];

if numel(dims) > 1
    DataLines = 4;
    %, header=[0, 1], index_col=0)
else
    DataLines = 2;
end

for index_typ = {'', 'x'}
    for weight_typ = {'', 'w'}
        fname = sprintf('ret_index_%s_all_grps.csv', joined_dims);
        fpath = fullfile(lib_comp_pycharm, 'index', fname);
        
        opts = detectImportOptions(fpath);
        opts.DataLines = [DataLines inf];
        opts.VariableNames = VariableNames;
        
        eval(sprintf('%sret_%sindex = readtable(fpath, opts);',...
            char(index_typ), char(weight_typ)))

    end
end

min_date = min(ret_index.datadate);
max_date = max(ret_index.datadate);

head(xret_index)

[MKT, ~] = load_data_files('MKT', fullfile(lib_data, 'fama_french_factors'), 'FRED');
MKT.Properties.VariableNames{1} = 'datadate';
[HML, ~] = load_data_files('HML', fullfile(lib_data, 'fama_french_factors'), 'FRED');
HML.Properties.VariableNames{1} = 'datadate';
[SMB, ~] = load_data_files('SMB', fullfile(lib_data, 'fama_french_factors'), 'FRED');
SMB.Properties.VariableNames{1} = 'datadate';




clc
if numel(dims) > 1
    for index_typ = {'', 'x'}
        for weight_typ = {'', 'w'}
            clear U S
            eval(sprintf('U = %sret_%sindex;',...
                char(index_typ), char(weight_typ)))

S = stack(U, 2:size(U, 2));
S.Properties.VariableNames(2:3) = {'keys', 'val'};

keys = cellfun(@(c) strsplit(c, '_')', cellstr(S.keys), 'UniformOutput', false);
keys = [keys{:}]';

S = [S, cell2table(keys, 'VariableNames', dims)];
S = S(:, [{'datadate'}, dims,{'val'}]);
S = varfun(@mean, S, 'InputVariables', 'val',...
    'GroupingVariables', {'datadate', dims{2}});
U = unstack(S, {'mean_val'}, dims(2));
U = removevars(U, 'GroupCount');

            eval(sprintf('%sret_%sindex = U;',...
                char(index_typ), char(weight_typ)))
    
        end
    end



end

clear U S



%% get portfolio beta, gen table6A-B, full sample

clc

% load stage1 regressor:
fpath_full_arima = fullfile(lib_comp, 'full_arima_101.csv');
full_arima = readtable(fpath_full_arima);
full_arima.Properties.VariableNames = [{'datadate'}, full_arima.Properties.VariableNames(2:end)];
%full_arima.loc[:, 'datadate'] = pd.to_datetime(full_arima.loc[:, 'datadate'])
% full_arima.loc[:, 'datadate'] = full_arima.loc[:, 'datadate'].apply(lambda t: t + relativedelta(months=1))
% full_arima.loc[:, 'datadate'] = full_arima.loc[:, 'datadate'] - pd.Timedelta(days=1)

freq = 'M';
model = 'inflation';
factor = 'epsilon_101';
y0 = 1967;
m0 = 7;
y1 = 2022;
m1 = 12;
mw = '_mw';
w = 'ew';

dates_CBS = '';
monthly_avg = '';

min_date = datetime(y0, m0, 1);
min_date = min_date + calmonths(1)-1;
max_date = datetime(y1, m1, 1);
max_date = max_date + calmonths(1)-1;

added_str_ = added_str(model, factor, min_date, max_date, freq, monthly_avg, mw, w, dates_CBS)
%added_str_ = added_str(model, factor, y0, m0, y1, m1, mw, w, freq, dates_CBS, monthly_avg)


idx_sample = ~isnat(xret_index.datadate);

idx_sample = idx_sample & xret_index.datadate >= min_date;
idx_sample = idx_sample & xret_index.datadate <= max_date;

data = xret_index(idx_sample, :);
grps = xret_index.Properties.VariableNames(2:end);

data = data(:, [{'datadate'}, flip(grps)]);
%data.HML_pi = data{:, end} - data{:, 2};
data.HML_pi = data{:, 2} - data{:, end};

data = innerjoin(data, full_arima(:, [{'datadate', factor}]), 'Keys', 'datadate');
data = innerjoin(data, MKT, 'Keys', 'datadate');
%data = innerjoin(data, Rf(:, [{'datadate', 'Rf'}]), 'Keys', 'datadate');


% Y = data.loc[:, 'xret'].values
B_grps = [];
t_stat_grps = [];
Exrets = [];
t_stat_Exrets = [];
Sharpe_ratios = [];
alpha_CAPMs = [];
t_alpha_CAPMs = [];

idx_regression = ~isnan(data{:, factor});
for grp = [grps, {'HML_pi'}]
    
    % regression of portfolio excess returns on inflation innovation
    X = data{idx_regression, factor};
    X = [ones(size(X, 1), 1), X];

    Y = data{idx_regression, grp};
    
    % regressions for EI
    idx = all(~isnan([X, Y]), 2);
    
    Y = Y(idx, :);
    X = X(idx, :);
    
    N = size(X, 1);
    K = size(X, 2);

    [B,BINT,R,~,STATS] = regress(Y, X);
    R_sq_adj = 1 - (1-STATS(1))*(N-1)/(N-K-1);

    sigma_hat = sqrt(R'*R/(N-K));
    D = X'*X; % information matrix
    
    F = 1/(sigma_hat^2 * K)*...
        B'*D*B;

    V_hat = sigma_hat^2*inv(D);
    STDB = sqrt(diag(V_hat));
    tstat_B = B./STDB;

    DW = sum(diff(R).^2)/sum(R.^2); % Durbin-Watson stat
    B_grps = [B_grps, B];
    t_stat_grps = [t_stat_grps, tstat_B];




    % regression of portfolio excess returns on MKT
    X = data{idx_regression, 'MKT'};
    X = [ones(size(X, 1), 1), X];

    Y = data{idx_regression, grp};
    
    % regressions for EI
    idx = all(~isnan([X, Y]), 2);
    
    Y = Y(idx, :);
    X = X(idx, :);
    
    N = size(X, 1);
    K = size(X, 2);

    [B,BINT,R,~,STATS] = regress(Y, X);
    R_sq_adj = 1 - (1-STATS(1))*(N-1)/(N-K-1);

    sigma_hat = sqrt(R'*R/(N-K));
    D = X'*X; % information matrix
    
    F = 1/(sigma_hat^2 * K)*...
        B'*D*B;

    V_hat = sigma_hat^2*inv(D);
    STDB = sqrt(diag(V_hat));
    tstat_B = B./STDB;

    DW = sum(diff(R).^2)/sum(R.^2); % Durbin-Watson stat

    alpha_CAPMs = [alpha_CAPMs, B(1, :)*12];
    t_alpha_CAPMs = [t_alpha_CAPMs, tstat_B(1, :)];



    % regression of portfolio excess returns on a constant
    X = data{idx_regression, {}};
    X = [ones(size(X, 1), 1), X];

    Y = data{idx_regression, grp};
    
    % regressions for EI
    idx = all(~isnan([X, Y]), 2);
    
    Y = Y(idx, :);
    X = X(idx, :);
    
    N = size(X, 1);
    K = size(X, 2);

    [B,BINT,R,~,STATS] = regress(Y, X);
    R_sq_adj = 1 - (1-STATS(1))*(N-1)/(N-K-1);

    sigma_hat = sqrt(R'*R/(N-K));
    D = X'*X; % information matrix
    
    F = 1/(sigma_hat^2 * K)*...
        B'*D*B;

    V_hat = sigma_hat^2*inv(D);
    STDB = sqrt(diag(V_hat));
    tstat_B = B./STDB;

    DW = sum(diff(R).^2)/sum(R.^2); % Durbin-Watson stat

    Exrets = [Exrets, B(1, :)*12];
    t_stat_Exrets = [t_stat_Exrets, tstat_B(1, :)];
    Sharpe_ratios = [Sharpe_ratios, sqrt(12)*B(1, :)/std(Y(idx, :))];



    
end


tab_6A = [B_grps(2, :); t_stat_grps(2, :)];
tab_6A = array2table(tab_6A, 'RowNames', {'beta_post', 't-stat_beta'},...
    'VariableNames', [grps, {'HML_pi'}]);

fname = ['tab_6A__', added_str_, '.xlsx'];
fpath_save = fullfile(lib_pycharm, 'Stocks_inflation', 'tabs', fname);

writetable(tab_6A, fpath_save, 'WriteRowNames', true)  % full sample


tab_6B = [Exrets; t_stat_Exrets; Sharpe_ratios; alpha_CAPMs; t_alpha_CAPMs];
tab_6B = array2table(tab_6B, 'RowNames',...
    {'Exret', 't_Exret', 'Sharpe_ratio', 'alpha_CAPM', 't_alpha_CAPM'},...
    'VariableNames', [grps, {'HML_pi'}]);

fname = ['tab_6BC__', added_str_, '.xlsx'];
fpath_save = fullfile(lib_pycharm, 'Stocks_inflation', 'tabs', fname);

writetable(tab_6B, fpath_save, 'WriteRowNames', true)  % full sample

tab_6A
tab_6B


%%
%{
***************************************************

                    PLOTS

***************************************************
%}
%% plot average return over period, all groups [for printing]
clc
clf
hold on
%min_date_sample = datetime(1967, 7, 31)
%max_date_sample = datetime(1999, 12, 31)

%min_date_sample = datetime(2000, 1, 1);
%max_date_sample = datetime(2022, 12, 31);

min_date_samples = [datetime(1967, 7, 31), datetime(2000, 1, 1)];
max_date_samples = [datetime(1999, 12, 31), datetime(2022, 12, 31)];

for i_fig = 1:2

min_date_sample = min_date_samples(i_fig);
max_date_sample = max_date_samples(i_fig);

idx_sample = true(size(xret_index.datadate));
idx_sample = idx_sample & xret_index.datadate >= min_date_sample;
idx_sample = idx_sample & xret_index.datadate <= max_date_sample;

plot(1:size(xret_index, 2)-1, mean(xret_index{idx_sample, 2:end}, 1))

end


min_date_samples.Format = 'MMM yy';
max_date_samples.Format = 'MMM yy';
leg0 = char(min_date_samples);
leg1 = char(max_date_samples);
%{
leg = [leg0, [' - '; ' - '], leg1, [' (\lambda='; ' (\lambda='],...
    char(arrayfun(@(c) num2str(c, '%.2f'),...
    price_of_risks, 'UniformOutput', false)),...
    [')'; ')']];
%}
leg = [leg0, [' - '; ' - '], leg1];


leg = cellstr(leg);
legend(leg)

xlabel('Inflation exposure decile')
ylabel('Average excess return (% p.m.)')

ax = gca;
ax.XTick = 1:size(xret_index, 2)-1;
%ax.XTick = ax.XTick + [.25, 0, -.1];
%ax.XTickLabel = {'Small', '', 'Big'};
ax.XTickLabel = xret_index.Properties.VariableNames(2:end);

%title('Average return', 'equally weighted index')
xlim([1, size(xret_index, 2)-1])

set(gcf, 'Position', [169 225 931 595]) % [width=0.88\textwidth]
set(gca, 'FontSize', 14)
set(gca, 'FontName', 'Lucida')
ax.XLabel.FontSize = 12;
ax.YLabel.FontSize = 12;

set(gcf, 'Color', 'w')

fpath_fig = fullfile(root, 'gfx', 'IRP_xax-grp.png');

export_fig(gcf, fpath_fig)


%% plot average return over period, all groups; xaxis is beta [for printing]
clc
clf
hold on

min_date_samples = [datetime(1967, 7, 31), datetime(2001, 1, 1)];
max_date_samples = [datetime(2000, 12, 31), datetime(2022, 12, 31)];

price_of_risks = [];
quant_of_risks = [];

for i_fig = 1:2

min_date_sample = min_date_samples(i_fig);
max_date_sample = max_date_samples(i_fig);


added_str_ = added_str(model, factor, min_date_sample, max_date_sample, freq, monthly_avg, mw, w, dates_CBS)
%added_str_ = added_str(model, factor, y0, m0, y1, m1, mw, w, freq, dates_CBS, monthly_avg)


idx_sample = ~isnat(xret_index.datadate);
idx_sample = idx_sample & xret_index.datadate >= min_date_sample;
idx_sample = idx_sample & xret_index.datadate <= max_date_sample;

data = xret_index(idx_sample, :);
grps = xret_index.Properties.VariableNames(2:end);

data = data(:, [{'datadate'}, flip(grps)]);
%data.HML_pi = data{:, end} - data{:, 2};
data.HML_pi = data{:, 2} - data{:, end};
data = innerjoin(data, MKT, 'Keys', 'datadate');
data = innerjoin(data, full_arima(:, [{'datadate', factor}]), 'Keys', 'datadate');
%data = innerjoin(data, Rf(:, [{'datadate', 'Rf'}]), 'Keys', 'datadate');




% Y = data.loc[:, 'xret'].values
B_grps = [];
t_stat_grps = [];
Exrets = [];
t_stat_Exrets = [];
Sharpe_ratios = [];
alpha_CAPMs = [];
t_alpha_CAPMs = [];

idx_regression = ~isnan(data{:, factor});


for grp = [grps, {'HML_pi'}]
    
    % regression of portfolio excess returns on a inflation innovations
    X = data{idx_regression, factor};
    X = [ones(size(X, 1), 1), X];

    Y = data{idx_regression, grp};
    
    % regressions for EI
    idx = all(~isnan([X, Y]), 2);
    
    Y = Y(idx, :);
    X = X(idx, :);
    
    N = size(X, 1);
    K = size(X, 2);

    [B,BINT,R,~,STATS] = regress(Y, X);
    R_sq_adj = 1 - (1-STATS(1))*(N-1)/(N-K-1);

    sigma_hat = sqrt(R'*R/(N-K));
    D = X'*X; % information matrix
    
    F = 1/(sigma_hat^2 * K)*...
        B'*D*B;

    V_hat = sigma_hat^2*inv(D);
    STDB = sqrt(diag(V_hat));
    tstat_B = B./STDB;

    DW = sum(diff(R).^2)/sum(R.^2); % Durbin-Watson stat
    B_grps = [B_grps, B];
    t_stat_grps = [t_stat_grps, tstat_B];


    % regression of portfolio excess returns on MKT
    X = data{idx_regression, 'MKT'};
    X = [ones(size(X, 1), 1), X];

    Y = data{idx_regression, grp};
    
    % regressions for EI
    idx = all(~isnan([X, Y]), 2);
    
    Y = Y(idx, :);
    X = X(idx, :);
    
    N = size(X, 1);
    K = size(X, 2);

    [B,BINT,R,~,STATS] = regress(Y, X);
    R_sq_adj = 1 - (1-STATS(1))*(N-1)/(N-K-1);

    sigma_hat = sqrt(R'*R/(N-K));
    D = X'*X; % information matrix
    
    F = 1/(sigma_hat^2 * K)*...
        B'*D*B;

    V_hat = sigma_hat^2*inv(D);
    STDB = sqrt(diag(V_hat));
    tstat_B = B./STDB;

    DW = sum(diff(R).^2)/sum(R.^2); % Durbin-Watson stat

    alpha_CAPMs = [alpha_CAPMs, B(1, :)*12];
    t_alpha_CAPMs = [t_alpha_CAPMs, tstat_B(1, :)];



    % regression of portfolio excess returns on a constant
    X = data{idx_regression, {}};
    X = [ones(size(X, 1), 1), X];

    Y = data{idx_regression, grp};
    
    % regressions for EI
    idx = all(~isnan([X, Y]), 2);
    
    Y = Y(idx, :);
    X = X(idx, :);
    
    N = size(X, 1);
    K = size(X, 2);

    [B,BINT,R,~,STATS] = regress(Y, X);
    R_sq_adj = 1 - (1-STATS(1))*(N-1)/(N-K-1);

    sigma_hat = sqrt(R'*R/(N-K));
    D = X'*X; % information matrix
    
    F = 1/(sigma_hat^2 * K)*...
        B'*D*B;

    V_hat = sigma_hat^2*inv(D);
    STDB = sqrt(diag(V_hat));
    tstat_B = B./STDB;

    DW = sum(diff(R).^2)/sum(R.^2); % Durbin-Watson stat

    Exrets = [Exrets, B(1, :)*12];
    t_stat_Exrets = [t_stat_Exrets, tstat_B(1, :)];
    Sharpe_ratios = [Sharpe_ratios, sqrt(12)* B(1, :)/std(Y(idx, :))];

    
end




idx_sample = true(size(xret_index.datadate));
idx_sample = idx_sample & xret_index.datadate >= min_date_sample;
idx_sample = idx_sample & xret_index.datadate <= max_date_sample;
IRP = mean(xret_index{idx_sample, 2:end}, 1) * 12;

quant_of_risk = B_grps(2, 1:end-1);

price_of_risk =...
    -( IRP(end)-IRP(1) )./...
    ( quant_of_risk(end)-quant_of_risk(1) );
price_of_risks = [price_of_risks; price_of_risk];
quant_of_risks = [quant_of_risks; quant_of_risk]

tab_6C = [Exrets; t_stat_Exrets; Sharpe_ratios; alpha_CAPMs; t_alpha_CAPMs];

tab_6C = array2table(tab_6C, 'RowNames',...
    {'Exret', 't_Exret', 'Sharpe_ratio', 'alpha_CAPM', 't_alpha_CAPM'},...
    'VariableNames', [grps, {'HML_pi'}]);

fname = ['tab_6BC__', added_str_, '.xlsx'];
fpath_save = fullfile(lib_pycharm, 'Stocks_inflation', 'tabs', fname);


writetable(tab_6C, fpath_save, 'WriteRowNames', true)  % subsample

[~, idx_sort] = sort(B_grps(2,1:end-1));

plot(B_grps(2,idx_sort), IRP(idx_sort))


end


min_date_samples.Format = 'MMM yy';
max_date_samples.Format = 'MMM yy';
leg0 = char(min_date_samples);
leg1 = char(max_date_samples);
leg = [leg0, [' - '; ' - '], leg1, [' (\lambda='; ' (\lambda='],...
    char(arrayfun(@num2str, price_of_risks, 'UniformOutput', false)),...
    [')'; ')']];
leg = cellstr(leg);
legend(leg)

ax = gca;
title('Average return', 'equally weighted index')
xlabel('Inflation beta')
%ylabel(['% p.', lower(freq), '.'])
%ylabel('%')
ylabel('% p.a.')

set(gcf, 'Position', [169 225 931 595]) % [width=0.88\textwidth]
set(gca, 'FontSize', 14)
set(gca, 'FontName', 'Lucida')

set(gcf, 'Color', 'w')


fpath_fig = fullfile(root, 'gfx', 'IRP_xax-beta.png');

export_fig(gcf, fpath_fig)

return
% Figure: quantity of risks, by period
clf
plot(1:size(xret_index, 2)-1, quant_of_risks(1,:))
ax = gca;
ax.YAxis.Color = ax.Children(1).Color;

yyaxis right
plot(1:size(xret_index, 2)-1, quant_of_risks(2,:))


yyaxis left
ax.YAxis(1).Color = ax.Children(1).Color;



ax = gca;
ax.XTick = 1:size(xret_index, 2)-1;
ax.XTickLabel = xret_index.Properties.VariableNames(2:end);

%title('Average return', 'equally weighted index')
xlim([1, size(xret_index, 2)-1])
xlabel('Portfolio, by inflation deciles')



min_date_samples.Format = 'MMM yy';
max_date_samples.Format = 'MMM yy';
leg0 = char(min_date_samples);
leg1 = char(max_date_samples);
%{
leg = [leg0, [' - '; ' - '], leg1, [' (\lambda='; ' (\lambda='],...
    char(arrayfun(@(c) num2str(c, '%.2f'),...
    price_of_risks, 'UniformOutput', false)),...
    [')'; ')']];
%}
leg = [leg0, [' - '; ' - '], leg1];


leg = cellstr(leg);
legend(leg, 'Location', 'northwest')

set(gcf, 'Position', [169 225 931 595]) % [width=0.88\textwidth]
set(gca, 'FontSize', 14)
set(gca, 'FontName', 'Lucida')

set(gcf, 'Color', 'w')

fpath_fig = fullfile(root, 'gfx', 'qrisk_xax-grp.png');

export_fig(gcf, fpath_fig)


















%% plot average return over period, all groups; xaxis is beta
clc
%min_date_sample = datetime(1967, 7, 31)
%max_date_sample = datetime(1999, 12, 31)

min_date_sample = datetime(2000, 1, 1)
max_date_sample = datetime(2022, 12, 31)

idx_sample = true(size(ret_index.datadate));
idx_sample = idx_sample & ret_index.datadate >= min_date_sample;
idx_sample = idx_sample & ret_index.datadate <= max_date_sample;

subplot(2,2,1)
plot(B_grps(2,1:end-1),12*mean(ret_index{idx_sample, 2:end}, 1))
title('Average return', 'equally weighted index')
xlabel('Inflation beta')
ylabel('% p.a.')


subplot(2,2,2)
plot(B_grps(2,1:end-1), 12*mean(xret_index{idx_sample, 2:end}, 1))
title('Average excess return', 'equally weighted index')
xlabel('Inflation beta')
ylabel('% p.a.')



subplot(2,2,3)
plot(B_grps(2,1:end-1), 12*mean(ret_windex{idx_sample, 2:end}, 1))
title('Average return', 'value weighted index')
xlabel('Inflation beta')
ylabel('% p.a.')



subplot(2,2,4)
plot(B_grps(2,1:end-1), 12*mean(xret_windex{idx_sample, 2:end}, 1))
title('Average excess return', 'value weighted index')
ylabel('% p.a.')

sgtitle({...
    'Portfolio performance, by inflation exposure',...
    sprintf('%s -- %s', char(min_date_sample), char(max_date_sample))})
xlabel('Inflation beta')

return








%% plot average return over period, all groups
clc
%min_date_sample = datetime(1967, 7, 31)
%max_date_sample = datetime(1999, 12, 31)

min_date_sample = datetime(2000, 1, 1)
max_date_sample = datetime(2022, 12, 31)

idx_sample = true(size(ret_index.datadate));
idx_sample = idx_sample & ret_index.datadate >= min_date_sample;
idx_sample = idx_sample & ret_index.datadate <= max_date_sample;

subplot(2,2,1)
plot(1:size(ret_index, 2)-1, mean(ret_index{idx_sample, 2:end}, 1))
ax = gca;
ax.XTick = 1:size(ret_index, 2)-1;
%ax.XTick = ax.XTick + [.25, 0, -.1];
%ax.XTickLabel = {'Small', '', 'Big'};
ax.XTickLabel = ret_index.Properties.VariableNames(2:end);
title('Average return', 'equally weighted index')
xlim([1, size(ret_index, 2)-1])



subplot(2,2,2)
plot(1:size(xret_index, 2)-1, mean(xret_index{idx_sample, 2:end}, 1))
ax = gca;
ax.XTick = 1:size(xret_index, 2)-1;
%ax.XTick = ax.XTick + [.25, 0, -.1];
ax.XTickLabel = ret_index.Properties.VariableNames(2:end);
title('Average excess return', 'equally weighted index')



subplot(2,2,3)
plot(1:size(ret_windex, 2)-1, mean(ret_windex{idx_sample, 2:end}, 1))
ax = gca;
ax.XTick = 1:size(ret_windex, 2)-1;
%ax.XTick = ax.XTick + [.25, 0, -.1];
ax.XTickLabel = ret_index.Properties.VariableNames(2:end);
title('Average return', 'value weighted index')



subplot(2,2,4)
plot(1:size(xret_windex, 2)-1, mean(xret_windex{idx_sample, 2:end}, 1))
ax = gca;
ax.XTick = 1:size(xret_windex, 2)-1;
%ax.XTick = ax.XTick + [.25, 0, -.1];
ax.XTickLabel = ret_index.Properties.VariableNames(2:end);
title('Average excess return', 'value weighted index')

sgtitle({...
    'Portfolio performance, by inflation exposure',...
    sprintf('%s -- %s', char(min_date_sample), char(max_date_sample))})

return









%%
mm=6; % movemean
clc

%period = 'early';
period = 'late';
print_each_subplot = true;

if strcmp(period, 'early')
    min_date_sample = datetime(1967, 7, 31);
    max_date_sample = datetime(1999, 12, 31);
    fpath_fig = fullfile(root, 'gfx', 'IHML and other factors - early period.png');

else
    min_date_sample = datetime(2000, 1, 1);
    max_date_sample = datetime(2022, 12, 31);
    fpath_fig = fullfile(root, 'gfx', 'IHML and other factors - late period.png');

end



params = containers.Map()
params('weight') = false;
params('constant') = true;

clf
for mm = 1:1
L = 0; % lag market


data = xret_index;
grps = xret_index.Properties.VariableNames(2:end);

data = data(:, [{'datadate'}, flip(grps)]);
%data.HML_pi = data{:, end} - data{:, 2};
data.HML_pi = data{:, 2} - data{:, end};

data = innerjoin(data, MKT);
data = innerjoin(data, HML);
data = innerjoin(data, SMB);

idx_sample = true(size(data.datadate));
idx_sample = idx_sample & data.datadate >= min_date_sample;
idx_sample = idx_sample & data.datadate <= max_date_sample;


if print_each_subplot
    clf
    fpath_fig = strrep(fpath_fig, 'period.png', 'period_1.png');
else
    subplot(3,3,1)
end

factors = {'HML_pi', 'HML','SMB', };
output=predict(data(idx_sample, :), 'MKT', factors, params);
plot(output.date, output{:, {'beta_1', 'beta_2', 'beta_3'}})

title(sprintf('%d factor model', numel(factors)))
legend(str4fig(factors))
ylim([-1 1.5])
ylabel('factor loading')

ax = gca;
Colors = reshape([ax.Children.Color]', 3, 3)';
Colors = Colors(end:-1:1, :);


yyaxis right
plot(output.date, output{:, {'R_sq'}},'k--')
ax = gca;
ax.YAxis(2).Color = 'k';
ax.Legend.String(end) = {'R\_sq (right axis)'};
ylim([0, 1])



if print_each_subplot
    fig = gcf;
    for child = fig.Children'
        set(child, 'FontSize', 12)
        set(child, 'FontName', 'Lucida')
    end
    
    set(gcf, 'Color', 'w')
    set(gcf, 'Position', [169 225 931 595]) % [width=0.88\textwidth]
    export_fig(gcf, fpath_fig)

    clf
    fpath_fig = strrep(fpath_fig, 'period_1.png', 'period_2.png');
else
    subplot(3,3,4)
end

factors = {'HML_pi', 'HML', };
output=predict(data(idx_sample, :), 'MKT', factors, params);
plot(output.date, output{:, {'beta_1', 'beta_2'}})

title(sprintf('%d factor model', numel(factors)))
legend(str4fig(factors))
ylim([-1 1.5])
ylabel('factor loading')

ax = gca;
ax.Children(1).Color = Colors(2, :);
ax.Children(2).Color = Colors(1, :);


yyaxis right
plot(output.date, output{:, {'R_sq'}},'k--')
ax = gca;
ax.YAxis(2).Color = 'k';
ax.Legend.String(end) = {'R\_sq (right axis)'};
ylim([0, 1])


if print_each_subplot
    fig = gcf;
    for child = fig.Children'
        set(child, 'FontSize', 12)
        set(child, 'FontName', 'Lucida')
    end
    
    set(gcf, 'Color', 'w')
    set(gcf, 'Position', [169 225 931 595]) % [width=0.88\textwidth]
    export_fig(gcf, fpath_fig)

    clf
    fpath_fig = strrep(fpath_fig, 'period_2.png', 'period_3.png');
else
    subplot(3,3,5)
end

factors = {'HML_pi', 'SMB', };
output=predict(data(idx_sample, :), 'MKT', factors, params);
plot(output.date, output{:, {'beta_1', 'beta_2'}})

ax = gca;
ax.Children(1).Color = Colors(3, :);
ax.Children(2).Color = Colors(1, :);

title(sprintf('%d factor model', numel(factors)))
legend(str4fig(factors))
ylim([-1 1.5])
ylabel('factor loading')


yyaxis right
plot(output.date, output{:, {'R_sq'}}, 'k--')
ax = gca;
ax.YAxis(2).Color = 'k';
ax.Legend.String(end) = {'R\_sq (right axis)'};
ylim([0, 1])


if print_each_subplot
    fig = gcf;
    for child = fig.Children'
        set(child, 'FontSize', 12)
        set(child, 'FontName', 'Lucida')
    end
    
    set(gcf, 'Color', 'w')
    set(gcf, 'Position', [169 225 931 595]) % [width=0.88\textwidth]
    export_fig(gcf, fpath_fig)

    clf
    fpath_fig = strrep(fpath_fig, 'period_3.png', 'period_4.png');
else
    subplot(3,3,6)
end

factors = {'HML', 'SMB', };
output=predict(data(idx_sample, :), 'MKT', factors, params);
plot(output.date, output{:, {'beta_1', 'beta_2'}})

ax = gca;
ax.Children(1).Color = Colors(3, :);
ax.Children(2).Color = Colors(2, :);

title(sprintf('%d factor model', numel(factors)))
legend(str4fig(factors))
ylim([-1 1.5])
ylabel('factor loading')


yyaxis right
plot(output.date, output{:, {'R_sq'}}, 'k--')
ax = gca;
ax.YAxis(2).Color = 'k';
ax.Legend.String(end) = {'R\_sq (right axis)'};
ylim([0, 1])


if print_each_subplot
    fig = gcf;
    for child = fig.Children'
        set(child, 'FontSize', 12)
        set(child, 'FontName', 'Lucida')
    end
    
    set(gcf, 'Color', 'w')
    set(gcf, 'Position', [169 225 931 595]) % [width=0.88\textwidth]
    export_fig(gcf, fpath_fig)

    clf
    fpath_fig = strrep(fpath_fig, 'period_4.png', 'period_5.png');
else
    subplot(3,3,7)
end

factors = {'HML_pi',};
output=predict(data(idx_sample, :), 'MKT', factors, params);
plot(output.date, output{:, {'beta_1'}})

ax = gca;
ax.Children(1).Color = Colors(1, :);


title(sprintf('%d factor model', numel(factors)))
legend(str4fig(factors))
ylim([-1 1.5])
ylabel('factor loading')

yyaxis right
plot(output.date, output{:, {'R_sq'}},'k--')
ax = gca;
ax.YAxis(2).Color = 'k';
ax.Legend.String(end) = {'R\_sq (right axis)'};
ylim([0, 1])



if print_each_subplot
    fig = gcf;
    for child = fig.Children'
        set(child, 'FontSize', 12)
        set(child, 'FontName', 'Lucida')
    end
    
    set(gcf, 'Color', 'w')
    set(gcf, 'Position', [169 225 931 595]) % [width=0.88\textwidth]
    export_fig(gcf, fpath_fig)

    clf
    fpath_fig = strrep(fpath_fig, 'period_5.png', 'period_6.png');
else
    subplot(3,3,8)
end

factors = {'HML',};
output=predict(data(idx_sample, :), 'MKT', factors, params);
plot(output.date, output{:, {'beta_1'}})

ax = gca;
ax.Children(1).Color = Colors(2, :);

title(sprintf('%d factor model', numel(factors)))
legend(str4fig(factors))
ylim([-1 1.5])
ylabel('factor loading')

yyaxis right
plot(output.date, output{:, {'R_sq'}},'k--')
ax = gca;
ax.YAxis(2).Color = 'k';
ax.Legend.String(end) = {'R\_sq (right axis)'};
ylim([0, 1])


if print_each_subplot
    fig = gcf;
    for child = fig.Children'
        set(child, 'FontSize', 12)
        set(child, 'FontName', 'Lucida')
    end
    
    set(gcf, 'Color', 'w')
    set(gcf, 'Position', [169 225 931 595]) % [width=0.88\textwidth]
    export_fig(gcf, fpath_fig)

    clf
    fpath_fig = strrep(fpath_fig, 'period_6.png', 'period_7.png');
else
    subplot(3,3,9)
end

factors = {'SMB',};
output=predict(data(idx_sample, :), 'MKT', factors, params);
plot(output.date, output{:, {'beta_1'}})

ax = gca;
ax.Children(1).Color = Colors(3, :);

title(sprintf('%d factor model', numel(factors)))
legend(str4fig(factors))
ylim([-1 1.5])
ylabel('factor loading')

yyaxis right
plot(output.date, output{:, {'R_sq'}},'k--')
ax = gca;
ax.YAxis(2).Color = 'k';
ax.Legend.String(end) = {'R\_sq (right axis)'};
ylim([0, 1])


min_date_sample.Format = 'MMM yyyy';
max_date_sample.Format = 'MMM yyyy';
tit0 = char(min_date_sample);
tit1 = char(max_date_sample);


fig = gcf;
for child = fig.Children'
    set(child, 'FontSize', 12)
    set(child, 'FontName', 'Lucida')
end

set(gcf, 'Color', 'w')
set(gcf, 'Position', [169 225 931 595]) % [width=0.88\textwidth]
export_fig(gcf, fpath_fig)
%sgtitle({'Factor models', [tit0, ' - ', tit1]})




return
clf
hold on
MKT_movmean = 12*movmean(data.MKT, mm);
HML_movmean = 12*movmean(data.HML_pi, mm);

[mm,...
    corr(data{idx_sample, {'MKT'}}, data{idx_sample, {'HML_pi'}}),...
    corr(MKT_movmean(idx_sample), HML_movmean(idx_sample))...
    ]
end


plot(data.datadate, MKT_movmean, 'k')
plot(data.datadate, HML_movmean)
legend({'MKT', 'HML'}, 'Location', 'northwest')


%{
Summary of results:
HML_inf tightly connected to MKT?
HML and SMB are not?!

no no no, MKT comes from a 3 factor model.
{MKT, HML,SMB} are independent by construction.
Need to run my own model, get SMB and HML 
(compare methodology via FF 3 factor model)

rerun test.

%}










