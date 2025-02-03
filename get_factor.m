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


Quality check: 28 Jun 2023
TODO!
Must check that all join keep rows sorted in a given way, e.g. gvkey, date
if this gets mixed up, all results are crap.
rerun all results after this. check if changed... if they did, throw away
all and redo all results.

Follow-up: 29 Jun 2023
I've added a couple of sortrows after every join of T_InputVariables
rerunning results

%}

if numel(InputVariables) > 1 && numel(bp)==1
    bp = repmat(bp, 1, numel(InputVariables));
elseif numel(InputVariables) ~= numel(bp)
    error('numel(InputVariables) ~= numel(bp) ~= 1')
end


nptiles_i = cellfun(@numel, bp);
nptiles_i = [nptiles_i, 1];
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

% set gvkey membership in ptile group at each date, for each InputVariables

% make sure group membership assigned even if on particular date a firm
% might have NaN value in any InputVariable
% use T_InputVariables_i to complete values based on previous non-NaN val
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
        
        % TODO: for 15 of the month this is incorrect
        T_InputVariables_i.datadate = T_InputVariables_i.datadate-1;
        T_InputVariables = [T_InputVariables; T_InputVariables_i];
        clear aux_
        clear idx_g
    end
else
    error('freq must be M or D.')
end

prctiles = table;
for i_InputVariables = 1:numel(InputVariables)
    InputVariable = InputVariables(i_InputVariables);
% prctiles gets a (1, nptiles+1) vector of percentile values for each date
prctiles_i=varfun(@(x)...
    prctile(x, [wins, bp{i_InputVariables}, 100-wins])-...
    [1e-10, zeros(1, numel(bp{i_InputVariables})+1)],... % adjustment for lowest entry
    T_InputVariables(:, [{'datadate'}, InputVariable]),...
    'InputVariables', InputVariable,...
    'GroupingVariables', 'datadate');
prctiles_i.Properties.VariableNames(end) = ...
     prctile_vars(i_InputVariables); % change name of variables
prctiles_i = removevars(prctiles_i, 'GroupCount');

if isempty(prctiles)
    prctiles = prctiles_i;
else
    prctiles = outerjoin(prctiles, prctiles_i, 'Keys', 'datadate', 'MergeKeys', true);
end
end

% add thresholds of group membership to T for each date
T_InputVariables = outerjoin(T_InputVariables, prctiles, 'Keys', 'datadate', 'MergeKeys', true, 'Type', 'left');
T_InputVariables = sortrows(T_InputVariables, {'gvkey', 'datadate'});

% GroupMembership collects results: group membership for each InputVariable
GroupMembership = NaN(size(T_InputVariables, 1), numel(InputVariables));
for i = 1:numel(InputVariables)

aux_ = NaN(size(T_InputVariables, 1), 1); % auxliary for each InputVariable

% index of group membership for each observation, 
% size(aux) = [size(T_InputVariables,1), nptiles]
aux = ...
    T_InputVariables{:, InputVariables{i}} > T_InputVariables{:, prctile_vars(i)}(:, 1:end-1) &...
    T_InputVariables{:, InputVariables{i}} <= T_InputVariables{:, prctile_vars(i)}(:, 2:end);

% remove obs not in any group
aux_(~any(aux, 2)) = 0;
% CHECK: potentially obselete, since I removed obs >(<) ub(lb)
% sum(~any(aux, 2)) == 0

[ia, ib] = find(aux); % group membership number (single dimension)
aux_(ia) = ib; % ib is the group, ia is the time-firm index
GroupMembership(:, i) = aux_; % add result to GroupMembership

clear aux aux_ ia ib
end

% collect overall group membership (effective if numel(InputVariables)>1)
GroupMembership_overall = NaN(size(GroupMembership, 1), 1);

% index of time-firm obs with no membership in at least one InputVariable
idx = any(GroupMembership == 0, 2);
% group membership=0 if firm-time obs has at least one dimesnion of 0
GroupMembership_overall(idx) = 0;

GroupMembership = GroupMembership(~idx,:); % remove membership=0 from GroupMembership
% break down to cell containing membership for each InputVariable
GroupMembership = mat2cell(GroupMembership, size(GroupMembership, 1), ones(1, size(GroupMembership,2)));

% determine overall group membership

%sub2ind(nptiles_i + 1, unique(GroupMembership{1}), unique(GroupMembership{end}));
%error('dd')
GroupMembership_overall(~idx) = sub2ind(nptiles_i+1, GroupMembership{:});



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


T_InputVariables = removevars(T_InputVariables, prctile_vars);
aux = T_InputVariables(idx, {'datadate', 'gvkey', 'grp'});
[g, g_] = findgroups(aux.gvkey);
grp = table;
for i_g=1:size(g_, 1)
    idx_g = g == i_g;
    grp_i = aux(idx_g, :);
    grp_i = table2timetable(grp_i);
    grp_i = retime(grp_i, 'daily', 'previous');
    grp_i = timetable2table(grp_i);
    
    grp = [grp; grp_i];
    clear idx_g grp_i
end
clear aux


T_InputVariables = removevars(T_InputVariables, 'grp');
T_InputVariables = outerjoin(T_InputVariables, grp,...
    'Type', 'left', 'Keys', {'datadate', 'gvkey'},...
    'MergeKeys', true);
T_InputVariables = sortrows(T_InputVariables, {'gvkey', 'datadate'});

T = outerjoin(T, T_InputVariables(:, {'datadate', 'gvkey', 'grp'}),...
    'Type', 'left', 'Keys', {'datadate', 'gvkey'},...
    'MergeKeys', true);

T_InputVariables.grp(isnan(T_InputVariables.grp))=0;
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
        get_indices(T_grp_i, 'reg', params);
    
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



end