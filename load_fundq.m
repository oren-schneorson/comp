function [T] = load_fundq(gvkey, fx_daily, data_dir)
%LOAD_FUNDQ Loads the fundq of a company for computing Merton DD
%   Input is gvkey, output is the relevant data from fundq, clean.

% to be able to not input data_dir
if nargin < 3
    data_dir = '/media/u70o/D/data/comp/fundq';
end

% to be able to input just data_dir='fundq' or 'g_fundq'
if isempty(fileparts(data_dir))
    data_dir = fullfile('/media/u70o/D/data/comp', data_dir);
end
    
fname = [char(gvkey), '.csv'];
fpath = fullfile(data_dir, fname);
T = readtable(fpath);

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


components_lctq = {'apq', 'lcoq', 'dlcq', 'txpq'};
components_ltq = {'lctq', 'loq', 'dlttq', 'txditcq'};
components_ltq_ = {'loq', 'dlttq', 'txditcq'};
% g_fundq doesn't have txpq, its included in lco
if contains(data_dir, 'g_fundq')
    components_lctq = components_lctq(1:end-1);
    components_ltq = components_ltq(1:end-1);
    components_ltq_ = components_ltq_(1:end-1);
end

% compute lctq when missing
%{
% for a subset of companies, compute lctq and ltq from components even when
% it is available.
subset = {'006385', '008381', '015399', '020629', '020877', '065840', '066072', '101317', '105460'};

if ismember(gvkey, subset)
    vals = sum(T{:, components_lctq}, 2, 'omitnan');
    vals(vals==0) = NaN;
    T{:, 'lctq'} = vals;
    vals = sum(T{:, components_ltq}, 2, 'omitnan');
    vals(vals==0) = NaN;
    T{:, 'ltq'} = vals;    
end
%}

idx = isnan(T.lctq) | T.lctq <= 0;
% first, via lctq components
vals = sum(T{idx, components_lctq}, 2, 'omitnan');
vals(vals==0) = NaN;
if numel(vals) > 0
    T{idx, 'lctq'} = vals;
end

idx = isnan(T.ltq) | T.ltq <= 0;
% second, ltq components
vals = sum(T{idx, components_ltq}, 2, 'omitnan');
vals(vals==0) = NaN;
if numel(vals) > 0
    T{idx, 'ltq'} = vals;
end


% third, if still missing data, via other components of ltq
idx = isnan(T.lctq);
vals = T.ltq(idx)-sum(T{idx, components_ltq_}, 2, 'omitnan');
vals(vals==0) = NaN;
if numel(vals) > 0
    T{idx, 'lctq'} = vals;
end

% omit data if one component is <= than zero
idx = T.lctq <= 0 | T.ltq <= 0;
T = T(~idx, :);

% crop
%T = T(:, {'datadate', 'ltq', 'lctq', 'loq', 'dlttq', 'txditcq', 'apq', 'lcoq', 'dlcq', 'txpq'});
if ~ismember({'rdq'}, T.Properties.VariableNames)
    T.rdq = T.pdateq;
end

if strcmp(class(T.rdq), 'double')
    if all(isnan(T.rdq))
        T.rdq = NaT(size(T.rdq));
    else
        error('rdq double, but not all is data missing.')
    end
end


idx = isnat(T.rdq);
T.rdq(idx) = T.datadate(idx) + caldays(1); 
T.rdq(idx) = T.rdq(idx) + calmonths(3); 
T.rdq(idx) = T.rdq(idx) - caldays(1); 

T = sortrows(T, {'datadate', 'rdq'});


[~, idx] = unique(T(:, {'datadate'}), 'rows');
T = T(idx, :);

T.datadate = T.rdq; % use rdq (reported date)
T = T(:, {'datadate', 'ltq', 'lctq', 'atq', 'actq'});


idx_col = varfun(@class, T);
idx_col = varfun(@(c) strcmp(c, 'double'), idx_col);
idx_col = idx_col{:, :};

idx =...
    find(~all(isnan(T{:, idx_col}), 2), 1, 'first'):...
    find(~all(isnan(T{:, idx_col}), 2), 1, 'last');
T = T(idx, :);

idx = ~all(isnan(T{:, idx_col}), 2);
T = T(idx, :);


%{
% to quarterly data, using previous observation.
[~, idx] = unique(T.datadate);
T = sortrows(T(idx,:), {'datadate'});

%fundq{'ltq', 'lctq'} = fillmissing(fundq{'ltq', 'lctq'}, 'previous');

% TODO: change this, it makes the dates wrong.
% better work with rdq, reported date... or pdateq for global firms g_fundq
T = table2timetable(T(idx,:));
T.datadate = dateshift(T.datadate, 'start', 'quarter');
T = retime(T, 'quarterly', 'previous');
T.datadate = dateshift(T.datadate, 'end', 'quarter');
T = timetable2table(T);
%}







end

