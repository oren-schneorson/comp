
if ~exist('var_type', 'var')
    run('var_types.m') % Variables types mapping
end

fx_dir = fullfile(lib_data, 'fx');

%%
fpath = fullfile(fx_dir, 'fx_daily.csv');
opts=detectImportOptions(fpath);
for v=1:numel(opts.VariableNames)
    %opts.VariableNames{v}
    try
        opts.VariableTypes{v} = var_type(opts.VariableNames{v});
    catch ME
        fprintf(opts.VariableNames{v})
        fprintf('\n')
        
    end
end
%%
%opts=setvaropts(opts, 'datadate', 'InputFormat', 'yyyy-MM-dd');
clc
fx_daily = readtable(fpath, opts);
%fx_daily = table2timetable(fx_daily);
%fx_daily = retime(fx_daily, 'daily', 'previous');
%fx_daily = timetable2table(fx_daily);

%fx_daily = fillmissing
%{
% DONE only once.
fx_daily = table2timetable(removevars(fx_daily, {'datadate'}), 'RowTimes', fx_daily.datadate);
fx_daily = retime(fx_daily, min(fx_daily.Time):caldays(1):max(fx_daily.Time), 'previous');
fx_daily = timetable2table(fx_daily);
fx_daily.Properties.VariableNames{1} = 'datadate';

% flatten table
n = size(fx_daily, 1);
fx_daily=arrayfun(@(i) [...
    fx_daily(:, 1),...
    table(repmat(fx_daily.Properties.VariableNames(i), n, 1), 'VariableNames', {'curnam'}),...
    table(fx_daily.(fx_daily.Properties.VariableNames{i}), 'VariableNames', {'curval'}),...
    ], 3:size(fx_daily, 2), 'UniformOutput', false);
fx_daily=vertcat(fx_daily{:});

%}
% Note: all fx rates are given in per unit of USD

fpath = fullfile(fx_dir, 'ISOs.csv');
opts=detectImportOptions(fpath);
for v=1:numel(opts.VariableNames)
    opts.VariableTypes{v} = var_type(opts.VariableNames{v}); 
end
opts=setvaropts(opts, 'From', 'InputFormat', 'yyyy-MM-dd');
opts=setvaropts(opts, 'To', 'InputFormat', 'yyyy-MM-dd');
ISOs = readtable(fpath, opts);

currency = containers.Map(ISOs.Code, ISOs.Current_currency);


