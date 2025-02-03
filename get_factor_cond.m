function [ret_idx_ew, ret_idx, xret_idx_ew, xret_idx, GroupCount, prctiles, T, mcaps, T_InputVariables, HML] = get_factor(T, InputVariables, bp, params)
%GET_FACTOR Finds the factor/s (robust for multiple dimensions
%   Uses InputVariable as the indexing variable. GroupingVariable is
%   usually date. Notes: the alogrithm assumes daily frequency. Firms with
%   less that 100 observations are dropped.
%
%   Inputs:
%   T               table of firm-level data
%   InputVariables  variables to factorize
%   bp              cell containing break-points for each InputVarible
%
%   Outputs:
%
%   ret_idx_ew      equally weighted returns for each
%   ret_idx         value weighted returns for each group
%   GroupCount      number of firms in each portfolio
%   prctiles        membership of ptile group at each date for each 
%                   InputVariables


%{
Quality check: 15 Dec 2022
Quality check: 30 May 2023: Group membership based on retimed table
It's simpler. Didn't resolve market 

1 Jun 2023: doesn't work for monthly data. grp empty
retime monthly gives 1st of the month...
%}


if numel(InputVariables) > 1 && numel(bp)==1
    bp = repmat(bp, 1, numel(InputVariables));
elseif numel(InputVariables) ~= numel(bp)
    error('numel(InputVariables) ~= numel(bp) ~= 1')
end


nptiles_i = cellfun(@numel, bp);
nptiles = prod(nptiles_i);

if ischar(InputVariables)
    InputVariables = {InputVariables};
end

wins = 0; % winsorize away extreme observations at each date;

% name of percentile variables
prctile_vars = strcat('prctile_', InputVariables);


for InputVariable = InputVariables
% create prctile membership at time t
if strcmp(InputVariable, 'b2m')
% Note: dropping negative book to market-equity observations

    ub = inf; % upper bound for book to market-equity ratio
    lb = -inf; % lower bound for book to market-equity ratio

elseif strcmp(InputVariable, 'beta_1') 

    ub = inf; % upper bound for betas
    lb = -inf; % lower bound for betas

elseif strcmp(InputVariable, 'beta_m')

    ub = inf; % upper bound for betas
    lb = -inf; % lower bound for betas

elseif strcmp(InputVariable, 'beta_Epi')

    ub = inf; % upper bound for betas
    lb = -inf; % lower bound for betas

else
   
    ub = inf;
    lb = -inf;
    
end

% removes NaN....
idx = (T{:, InputVariable} >= lb) & (T{:, InputVariable} <= ub);

% remove obs for using conditions on ub and lb, for each InputVariable
T = T(idx, :);

end

T_InputVariables = T(:, [{'datadate', 'gvkey'}, InputVariables]);


% set gvkey membership in ptile group at each date, for each InputVariables

% make sure group membership assigned even if on particular date a firm
% might have NaN value in any InputVariable
% use T_InputVariables_i to complete values based on previous non-NaN val

% TODO: can this be done more compactly with varfun + retime, 
% grouping on gvkey?
if strcmp(params('freq'), 'D')
    T_InputVariables = table;
    [g, g_] = findgroups(T.gvkey);
    for i_g = 1:size(g_,1)
        clc
        fprintf('Complete values for InputVariables in daily freq, %.2f%%\n',...
        i_g/size(g_,1)*100)
        idx_g = g == i_g;
        T_InputVariables_i = table2timetable(T(idx_g, [{'datadate', 'gvkey'}, InputVariables]));
        T_InputVariables_i = retime(T_InputVariables_i, 'daily', 'previous');
        T_InputVariables_i = timetable2table(T_InputVariables_i);
        T_InputVariables = [T_InputVariables; T_InputVariables_i];
        clear aux_
        clear idx_g
    end
elseif strcmp(params('freq'), 'M')

    T_InputVariables = table;
    [g, g_] = findgroups(T.gvkey);
    for i_g = 1:size(g_,1)
        clc
        fprintf('Complete values for InputVariables in monthly freq, %.2f%%\n',...
        i_g/size(g_,1)*100)
        idx_g = g == i_g;
        T_InputVariables_i = table2timetable(T(idx_g, [{'datadate', 'gvkey'}, InputVariables]));
        T_InputVariables_i = retime(T_InputVariables_i, 'monthly', 'previous');
        T_InputVariables_i = timetable2table(T_InputVariables_i);
        T_InputVariables_i = T_InputVariables_i(~isnan(T_InputVariables_i.gvkey), :);
        
        % TODO: for 15 of the month this is incorrect
        T_InputVariables_i.datadate = T_InputVariables_i.datadate-1;
        
        T_InputVariables = [T_InputVariables; T_InputVariables_i];
        clear aux_
        clear idx_g
    end
else
    error('freq must be M or D.')
end


% GroupMembership: collects results of group membership for each InputVariable
GroupMembership = zeros(size(T_InputVariables, 1), numel(InputVariables));

% prctiles: collects percentiles for each InputVariable, conditional on
% membership in former groups of InputVariables. It's rows are unique dates
prctiles = table;
idx_wf0 = true(size(T_InputVariables, 1), 1);

params_wf = containers.Map();
params_wf('InputVariables') = InputVariables;
params_wf('bp') = bp;
params_wf('prctile_vars') = prctile_vars;
params_wf('wins') = wins;

% waterfall of group membership: conditional on being a member of former
% groups, what is the group a firm is in in a given time.
[prctiles, GroupMembership] = get_factor_wf(1, idx_wf0, 1, T_InputVariables, prctiles, GroupMembership, params_wf);
%sum(GroupMembership) == [0,0]


%{
% set NaN to 0 (or simply GroupMembership delcated zeros)
[ia, ib] = find(isnan(GroupMembership));
GroupMembership(ia, ib) = 0;
%}

% add thresholds of group membership to T for each date
T_InputVariables = outerjoin(T_InputVariables, prctiles,...
    'Keys', 'datadate', 'MergeKeys', true, 'Type', 'left');



T_InputVariables = sortrows(T_InputVariables, {'gvkey', 'datadate'});

% collect overall group membership (effective if numel(InputVariables)>1)
GroupMembership_overall = zeros(size(GroupMembership, 1), 1);

%{
% similar to before, simply declare it as zeros
% index of time-firm obs with no membership in at least one InputVariable
idx = any(GroupMembership == 0, 2);
% group membership=0 if firm-time obs has at least one dimesnion of 0
GroupMembership_overall(idx) = 0;
%}

idx = any(GroupMembership == 0, 2);
GroupMembership_ = GroupMembership(~idx,:); % remove membership=0 from GroupMembership
% break down to cell containing membership for each InputVariable
GroupMembership_ = mat2cell(GroupMembership_, size(GroupMembership_, 1), ones(1, size(GroupMembership_,2)));

% determine overall group membership
GroupMembership_overall(~idx) = sub2ind(nptiles_i+1, GroupMembership_{:});
clear GroupMembership_

T_InputVariables.grp = GroupMembership_overall; % add to T_InputVariables
clear aux_


unq_dates = unique(T_InputVariables.datadate);

if strcmp(params('freq'), 'D')
% for daily freq, end of month values
    idx = ...
        ismember(...
        T_InputVariables.datadate,...
        unique(eomdate(unq_dates)));
elseif strcmp(params('freq'), 'M')
% for monthly freq, end of year values
    idx = ...
        T_InputVariables.datadate.Month == 12 &...
        T_InputVariables.datadate.Day == 31;
end

% if firm doesn't have a valid group at eom... not sure if to include this
%T_InputVariables.grp = fillmissing(T_InputVariables.grp, 'previous');

aux = T_InputVariables(idx, {'gvkey', 'datadate', 'grp'});
[g, g_] = findgroups(aux.gvkey);
grp = table;
%dataloss = table; % debug
for i_g=1:size(g_, 1)
    idx_g = g == i_g;

    grp_i = aux(idx_g, :);
    grp_i = table2timetable(grp_i);
    grp_i = retime(grp_i, 'daily', 'previous');
    grp_i = timetable2table(grp_i);
    grp = [grp; grp_i];
    
    %{
    % debug, dataloss
    gvkey = g_(i_g);
    idx_gvkey = T_InputVariables.gvkey == gvkey;
    dates_lost_i = setdiff(T_InputVariables.datadate(idx_gvkey), grp_i.datadate);
    time_lost_i = calmonths(between(min(dates_lost_i), max(dates_lost_i)));
    
    dataloss_i = numel(dates_lost_i)/sum(idx_gvkey);
    dataloss_i = [gvkey, time_lost_i, sum(idx_gvkey), dataloss_i];
    dataloss_i = array2table(dataloss_i, 'VariableNames',...
        {'gvkey', 'months_lost', 'nobs', 'loss'});
    
    dataloss = [dataloss; dataloss_i];
    clear idx_gvkey dataloss_i
    %}
    clear idx_g grp_i
end
clear aux

% sum(grp.grp == 0) == 0



T_InputVariables = removevars(T_InputVariables, 'grp');

% this is where the loss of data occurs -->
T_InputVariables = outerjoin(T_InputVariables, grp,...
    'Type', 'left', 'Keys', {'gvkey', 'datadate'},...
    'MergeKeys', true);
T_InputVariables = sortrows(T_InputVariables, {'gvkey', 'datadate'});


%{
    2023-06-28
********************
    NOTE ON DATALOSS
********************
up to here there is 0.9% loss, meaning an average firm has 111
months of data, which is close to 160 months but not close enough.
I'm losing 37,518 obs because grp table does not have matching
Keys (gvkey, date) that T_InputVariables does have.

This precisely equals sum(isnan(T_InputVariables.grp))...

Observe that:
1. grp is built in a way to preserve eomdate,
2. T_InputVariables is retimed to daily in the top...

I'm losing twice as much data as I need to lose, because loss occurs
both at the *end* of firm life and at the beginning.

Fill T_InputVariables.grp with previous by gvkey.
this roughly about halve the dataloss, gaining obs lost in the end 
of data that didn't end at eom... This is done just below.
%}
    
% make sure rows are sorted.
T_InputVariables = sortrows(T_InputVariables, {'gvkey', 'datadate'});

missing_vals=varfun(@(x) fillmissing(x, 'previous'), T_InputVariables,...
    'InputVariables', {'datadate', 'grp'}, 'GroupingVariables', 'gvkey');
missing_vals = removevars(missing_vals, 'GroupCount');
missing_vals = renamevars(missing_vals, 'Fun_datadate', 'datadate');
missing_vals = renamevars(missing_vals, 'Fun_grp', 'grp');

% check that missing_vals has exactly the same (gvkey, date) pairs. V
% this should be true because T_InputVariables is sorted by gvkey, date
if ~all(...
    T_InputVariables.gvkey == missing_vals.gvkey & ...
    T_InputVariables.datadate == missing_vals.datadate)

    error('missing_vals and T_InputVariables do not have the same rows.')

end

% fillmissing grp values
T_InputVariables.grp = missing_vals.grp;
% the rest of the observations that are NaN filled with 0.
T_InputVariables.grp(isnan(T_InputVariables.grp))=0;



T = outerjoin(T, T_InputVariables(:, {'datadate', 'gvkey', 'grp'}),...
    'Type', 'left', 'Keys', {'datadate', 'gvkey'},...
    'MergeKeys', true);
T.grp(isnan(T.grp))=0;


unq_i = unique(T.grp(T.grp>0))'; % unique groups

ret_idx = table; % declare table for index of firms by prctile group
ret_idx_ew = table; % declare table for equally weighted index
xret_idx = table; % declare table for index of firms by prctile group
xret_idx_ew = table; % declare table for equally weighted index
GroupCount = table;

for i = unq_i
    clc

    fprintf('Computing prctile group %d of %d', i, numel(unq_i))
    T_grp_i = T(T.grp == i, :); % data for each group
        
    [~, ret_idx_ew_grp_i,...
        ret_idx_grp_i,...
        xret_idx_ew_grp_i,...
        xret_idx_grp_i,...
        ~] = ...
        get_indices(T_grp_i, 'adj', params);
    
    GroupCount_prctile_i = ret_idx_ew_grp_i(:, {'datadate', 'GroupCount'});
    ret_idx_ew_grp_i = removevars(ret_idx_ew_grp_i, 'GroupCount');
    ret_idx_grp_i = removevars(ret_idx_grp_i, 'GroupCount');
    xret_idx_ew_grp_i = removevars(xret_idx_ew_grp_i, 'GroupCount');
    xret_idx_grp_i = removevars(xret_idx_grp_i, 'GroupCount');
    
    % adjust variable names
    ret_idx_ew_grp_i.Properties.VariableNames{end} = [...
        ret_idx_ew_grp_i.Properties.VariableNames{end},...
        sprintf('_%d', i)];
    ret_idx_grp_i.Properties.VariableNames{end} = [...
        ret_idx_grp_i.Properties.VariableNames{end},...
        sprintf('_%d', i)];
    xret_idx_ew_grp_i.Properties.VariableNames{end} = [...
        xret_idx_ew_grp_i.Properties.VariableNames{end},...
        sprintf('_%d', i)];
    xret_idx_grp_i.Properties.VariableNames{end} = [...
        xret_idx_grp_i.Properties.VariableNames{end},...
        sprintf('_%d', i)];
    GroupCount_prctile_i.Properties.VariableNames{end} = [...
        GroupCount_prctile_i.Properties.VariableNames{end},...
        sprintf('_%d', i)];
    

    if isempty(ret_idx_ew)
        ret_idx_ew = ret_idx_ew_grp_i;
        ret_idx = ret_idx_grp_i;
        xret_idx_ew = xret_idx_ew_grp_i;
        xret_idx = xret_idx_grp_i;
        GroupCount = GroupCount_prctile_i;
    else
        ret_idx_ew = outerjoin(ret_idx_ew, ret_idx_ew_grp_i, 'Keys', 'datadate', 'MergeKeys', true);
        ret_idx = outerjoin(ret_idx, ret_idx_grp_i, 'Keys', 'datadate', 'MergeKeys', true);
        xret_idx_ew = outerjoin(xret_idx_ew, xret_idx_ew_grp_i, 'Keys', 'datadate', 'MergeKeys', true);
        xret_idx = outerjoin(xret_idx, xret_idx_grp_i, 'Keys', 'datadate', 'MergeKeys', true);
        GroupCount = outerjoin(GroupCount, GroupCount_prctile_i, 'Keys', 'datadate', 'MergeKeys', true);
    end
    
    

end

% compute total market cap. of each group
% smooth mcap aggregate
mcaps = table;
[g, g_] = findgroups(T.gvkey);
for i_g = 1:size(g_,1)
    idx = g == i_g;
    mcaps_i = T(idx, {'datadate', 'gvkey', 'grp', 'mcap'});
    mcaps_i = table2timetable(mcaps_i);
    mcaps_i = retime(mcaps_i, 'daily', 'previous');
    mcaps_i = timetable2table(mcaps_i);
    mcaps = [mcaps; mcaps_i];
    clear mcaps_i
end


% compute mcap by group
mcaps = mcaps(mcaps.grp ~= 0, :);
mcaps = varfun(@(x) sum(x, 'omitnan'), mcaps,...
    'InputVariable', 'mcap',...
    'GroupingVariable', {'datadate', 'grp'});
mcaps = removevars(mcaps, 'GroupCount');

mcaps = unstack(mcaps, 'Fun_mcap', 'grp');
mcaps.Properties.VariableNames(2:end) = ...
    arrayfun(@(i) sprintf('mcap_%d', i), 1:size(mcaps,2)-1,...
    'UniformOutput', false);


HML = ret_idx_ew(:, 1);
HML.ret_idx = ret_idx{:, end}-ret_idx{:, 2};
HML.ret_idx_ew = ret_idx_ew{:, end}-ret_idx_ew{:, 2};


% preserve group membership by InputVariable
GroupMembership = array2table(GroupMembership, 'VariableNames', strcat('grp_', InputVariables));
T_InputVariables = [T_InputVariables, GroupMembership];


end