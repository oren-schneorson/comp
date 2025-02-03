clc


series = 'adj_prccd';
CPI_var = 'CPI_SA';
fname = 'SP500_YAHOO.D.xlsx';

fpath = fullfile(lib_data, 'stock_markets', fname);

opts = detectImportOptions(fpath);
opts.VariableNames(1) = {'datadate'};
opts.VariableTypes(1) = {'datetime'};
opts.VariableTypes(2:end) = repmat({'double'}, 1, numel(opts.VariableTypes)-1);


opts = setvaropts(opts,'datadate', ...
                       'DatetimeFormat',date_format);
opts = setvaropts(opts,{'datadate'}, ...
                       'InputFormat',date_format);

                   
mkt = readtable(fpath, opts);

Keys = {'datadate'};
mkt = outerjoin(mkt, CPI_D(:, [Keys, {'CPI_SA'}]), 'Keys', Keys, 'MergeKeys', true, 'Type', 'left');


R_f_ = R_f;
R_f_{:,end} = R_f_{:,end}/100+1;
R_f_{:,end} = R_f_{:,end}.^(1/252);
R_f_{:,end} = (R_f_{:,end} - 1)*100;


mkt = outerjoin(mkt, R_f_,...
    'Keys', 'datadate', 'MergeKeys', true,...
    'Type', 'left');

mkt.ret_m = [NaN; mkt{2:end, series}./mkt{1:end-1, series}];
mkt.retr_m = ...
    [NaN; mkt{2:end, series}./...
          mkt{1:end-1, series}]./(...
    [NaN; mkt{2:end, CPI_var}./...
          mkt{1:end-1, CPI_var}]);

% returns in percent
mkt.ret_m = (mkt.ret_m-1)*100;
mkt.retr_m = (mkt.retr_m-1)*100;

% excess return
mkt.xret_m = mkt.ret_m-mkt{:, R_f.Properties.VariableNames(end)};
mkt = renamevars(mkt, 'mcap', 'mcap_m');
mkt = renamevars(mkt, 'adj_prccd', 'adj_prccd_m');

vars_mkt = {'datadate', 'ret_m', 'xret_m', 'retr_m', 'adj_prccd_m', 'mcap_m'};
mkt = mkt(:, vars_mkt);



