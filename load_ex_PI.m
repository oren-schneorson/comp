
clc
smooth_ex_PI = false;
date_format = 'yyyy-MM-dd';

% expected inflation from capital markets (data from FAME)
fname = 'ex_PI_Cleveland.xlsx';
fpath_ex_PI_Cleveland = fullfile(lib_data, 'ex_PI', fname);

% Cleveland
series = 'Cleveland';
fname = sprintf('ex_PI_%s_metadata.xlsx', series);
fpath = fullfile(lib_data, 'ex_PI', fname);
metadata_Cleveland = readtable(fpath);

metadata_Cleveland.start_date = datetime(metadata_Cleveland.start_date, 'InputFormat', date_format);
metadata_Cleveland.end_date = datetime(metadata_Cleveland.end_date, 'InputFormat', date_format);

metadata_Cleveland.start_date = metadata_Cleveland.start_date-1;
metadata_Cleveland.end_date = metadata_Cleveland.end_date-1;

fname = sprintf('ex_PI_%s.xlsx', series);
fpath = fullfile(lib_data, 'ex_PI', fname);
opts = detectImportOptions(fpath);

opts.VariableNames(1) = {'date'};
opts.VariableTypes(1) = {'datetime'};
opts.VariableTypes(2:end) = repmat({'double'}, 1, numel(opts.VariableTypes)-1);


opts = setvaropts(opts,'date', ...
                       'DatetimeFormat',date_format);
opts = setvaropts(opts,{'date'}, ...
                       'InputFormat',date_format);

ex_PI_Cleveland = readtable(fpath, opts);
ex_PI_Cleveland.date = ex_PI_Cleveland.date-1;


% TIPS
series = 'TIPS';

fname = sprintf('ex_PI_%s_metadata.xlsx', series);
fpath = fullfile(lib_data, 'ex_PI', fname);
metadata_TIPS = readtable(fpath);

fname = sprintf('ex_PI_%s.xlsx', series);
fpath = fullfile(lib_data, 'ex_PI', fname);
opts = detectImportOptions(fpath);

opts.VariableNames(1) = {'date'};
opts.VariableTypes(1) = {'datetime'};
opts.VariableTypes(2:end) = repmat({'double'}, 1, numel(opts.VariableTypes)-1);


opts = setvaropts(opts,'date', ...
                       'DatetimeFormat',date_format);
opts = setvaropts(opts,{'date'}, ...
                       'InputFormat',date_format);


ex_PI_TIPS = readtable(fpath, opts);
ex_PI_TIPS{:, 2:end} = fillmissing(ex_PI_TIPS{:, 2:end}, 'previous');

%idx = all(isnan(ex_PI_TIPS{:,2:end}), 2);
%ex_PI_TIPS = ex_PI_TIPS(~idx, :);
%idx = all(isnan(ex_PI_Cleveland{:,2:end}), 2);
%ex_PI_Cleveland = ex_PI_Cleveland(~idx, :);


% create first difference variables: ex_PI_TIPS
for var = ex_PI_TIPS.Properties.VariableNames(2:end)
    var = char(var);
    newvar = ['d', var];
    
    if smooth_ex_PI
        idx = true(size(ex_PI_TIPS, 1), 1);
        idx = idx & ~isnan(ex_PI_TIPS{:, var});
        idx = idx & ~isinf(ex_PI_TIPS{:, var});        
        ex_PI_TIPS{idx, var} = hpfilter(ex_PI_TIPS{idx, var}, 1600);
        
    end

    ex_PI_TIPS.(newvar) = [NaN; diff(ex_PI_TIPS.(var))];

end



% create first difference variables
for var = ex_PI_Cleveland.Properties.VariableNames(2:end)
    var = char(var);
    newvar = ['d', var];
    
    if smooth_ex_PI        
        idx = true(size(ex_PI_Cleveland, 1), 1);
        idx = idx & ~isnan(ex_PI_Cleveland{:, var});
        idx = idx & ~isinf(ex_PI_Cleveland{:, var});
        ex_PI_Cleveland{idx, var} = hpfilter(ex_PI_Cleveland{idx, var}, 1600);

    end

    ex_PI_Cleveland.(newvar) = [NaN; diff(ex_PI_Cleveland.(var))];
end

ex_PI_TIPS = renamevars(ex_PI_TIPS, 'date', 'datadate');
ex_PI_Cleveland = renamevars(ex_PI_Cleveland, 'date', 'datadate');




