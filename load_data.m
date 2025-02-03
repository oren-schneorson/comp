%{
Series of scripts to load stock market data and compute relevant measures.

data being loaded:
CPI, consumer price index (Israel)
ex_PI (expected / break-even inflation expectation, BEI)
mkt (market data)
Risk free rate: DTB3
Exchange rates: fx_daily
Last quality check: 2 Jul 2023

%}

clc
% CPI (to get real price series)
fpath_m = fullfile(root, 'load_CPI.m');
run(fpath_m)


% CPI (to get real price series)
fpath_m = fullfile(root, 'load_currency.m');
run(fpath_m)

% expected inflation (needed for applications with expectations data)
fpath_m = fullfile(root, 'load_ex_PI.m');
run(fpath_m)
% Risk free rate
fpath_m = fullfile(root, 'load_Rf.m');
run(fpath_m)


% market prices and returns
fpath_m = fullfile(root, 'load_mkt.m');
run(fpath_m)
%%
% quarterly reports (file from Artium)
%fpath_m = fullfile(root, 'load_quarterly_reports.m');
%run(fpath_m)
%% TODO: generate full ARIMA to the U.S.
%fpath = './ARMA(1,1)_US.csv';
%innov = readtable(fpath);

fpath = './full_arima_101.csv';
full_arima = readtable(fpath);
full_arima = renamevars(full_arima, 'date', 'datadate');
%full_arima.datadate = full_arima.datadate-caldays(1);
% lag inflation one month
full_arima.datadate+caldays(1)+calmonths(1)-caldays(1);
%head(full_arima)

%%
T = table;
%{
if exist(fpath_equity, 'file') && load_fpath_equity
    T = readtable(fpath_equity);
else

    % load gvkeys
    fpath_m = fullfile(root, 'load_gvkeys.m');
    run(fpath_m)
    
    writetable(T, fpath_equity)
end

%}