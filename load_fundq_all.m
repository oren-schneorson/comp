function [T] = load_fundq_all(gvkey, fx_daily, data_dir, var_type)
%LOAD_FUNDQ_ALL Loads the fundq of a company.
%   Input is gvkey, output is the relevant data from fundq, clean.

DateFormat = 'yyyy-MM-dd';            

% to be able to not input data_dir
if nargin < 3
    username = getenv("USER")
    data_dir = fullfile('/media', username, 'D', 'data', 'comp', 'fundq');
end

% to be able to input just data_dir='fundq' or 'g_fundq'
if isempty(fileparts(data_dir))
    username = getenv("USER")
    data_dir = fullfile('/media', username, 'D', 'data', 'comp', data_dir);
end
    
fname = [char(gvkey), '.csv'];
fpath = fullfile(data_dir, fname);
opts = gen_opts(fpath, {'datadate', 'rdq', 'fdateq', 'pdateq'}, DateFormat, var_type);


T = readtable(fpath, opts);

if isempty(T)
    return
end


% 21-03-2023: fx_daily is not up to date, losing some data in the end
if ~all(ismember(T.curcdq, {'USD'}))
    classes = varfun(@class, T, 'OutputFormat', 'cell');
    idx_classes = find(strcmp(classes, 'double'));

    T = innerjoin(...
        T,...
        fx_daily(:, {'datadate', 'curcdq', 'currency_val'}),...
        'Keys', {'datadate', 'curcdq'});
    T{:, idx_classes} = T{:, idx_classes}./T.currency_val;
end










end

