function [T, ret_idx_ew, ret_idx, xret_idx_ew, xret_idx, mcap_idx] = get_indices(T, InputVariable, params)
%GET_INDICES Computes return index
%   Uses InputVariable as the indexing variable. GroupingVariable is
%   usually date



if strcmp(params('freq'), 'D')
    var_ret = 'retd';
    var_xret = 'xretd';
    var_retr = 'retrd';
    idx = abs(T.retd) < 20;
elseif strcmp(params('freq'), 'M')
    idx = abs(T.retm) < 240;
    var_ret = 'retm';
    var_xret = 'xretm';
    var_retr = 'retrm';
else
    error('freq must be D or M.')
end

T = T(idx, :);

GroupingVariable = {'datadate'};

if strcmp(InputVariable, 'reg')
    % regular value weight, no constraints
    InputVariables = {'mcap'};

elseif strcmp(InputVariable, 'adj')
    
    % adjusted value weight, constraints
    InputVariables = {'mcap_adj'};
    T.mcap_adj = T.cshoc_index .* T.prcbd .* T.pubow/100;
    T = get_weight_index(T, InputVariables, GroupingVariable);

end

T = get_weight_index(T, InputVariables, GroupingVariable);

T.(['w', var_ret]) = T.w .* T.(var_ret);
T.(['w', var_xret]) = T.w .* T.(var_xret);
T.(['w', var_retr]) = T.w .* T.(var_retr);

T.wmcap = T.w .* T.mcap;
%T.wretd = T.w .* T.retd;
%T.wretrd = T.w .* T.retrd;
%T.wxretd = T.w .* T.xretd;

InputVariables = {'wmcap'};
aux = varfun(@(x) sum(x, 'omitnan'), T,...
    'InputVariables', InputVariables,...
    'GroupingVariables', GroupingVariable);
aux.Properties.VariableNames(end) =...
    {'mcap_idx'};

aux = removevars(aux, 'GroupCount');
mcap_idx = aux;

vars_T = ~ismember(T.Properties.VariableNames,...
    setdiff(aux.Properties.VariableNames, GroupingVariable));

T = outerjoin(T(:, vars_T), aux, 'Keys', 'datadate',...
    'MergeKeys', true, 'Type', 'left');

InputVariables = {['w', var_ret]};
aux = varfun(@(x) sum(x, 'omitnan'), T,...
    'InputVariables', InputVariables,...
    'GroupingVariables', GroupingVariable);
aux.Properties.VariableNames(end) =...
    {'ret_idx'};
ret_idx = aux;

aux = removevars(aux, 'GroupCount');
vars_T = ~ismember(T.Properties.VariableNames,...
    setdiff(aux.Properties.VariableNames, GroupingVariable));

T = outerjoin(T(:, vars_T), aux, 'Keys', 'datadate',...
    'MergeKeys', true, 'Type', 'left');


InputVariables = {['w', var_xret]};
aux = varfun(@(x) sum(x, 'omitnan'), T,...
    'InputVariables', InputVariables,...
    'GroupingVariables', GroupingVariable);
aux.Properties.VariableNames(end) =...
    {'xret_idx'};
xret_idx = aux;

aux = removevars(aux, 'GroupCount');
vars_T = ~ismember(T.Properties.VariableNames,...
    setdiff(aux.Properties.VariableNames, GroupingVariable));

T = outerjoin(T(:, vars_T), aux, 'Keys', 'datadate',...
    'MergeKeys', true, 'Type', 'left');



InputVariables = {['w', var_retr]};
aux = varfun(@(x) sum(x, 'omitnan'), T,...
    'InputVariables', InputVariables,...
    'GroupingVariables', GroupingVariable);
aux.Properties.VariableNames(end) =...
    {'retr_idx'};
%wretr_idx = aux;

aux = removevars(aux, 'GroupCount');
vars_T = ~ismember(T.Properties.VariableNames,...
    setdiff(aux.Properties.VariableNames, GroupingVariable));

T = outerjoin(T(:, vars_T), aux, 'Keys', 'datadate',...
    'MergeKeys', true, 'Type', 'left');

InputVariables = {var_ret};
aux = varfun(@(x) mean(x, 'omitnan'), T,...
    'InputVariables', InputVariables,...
    'GroupingVariables', GroupingVariable);
aux.Properties.VariableNames(end) =...
    {'ret_idx_ew'};
ret_idx_ew = aux;

aux = removevars(aux, 'GroupCount');
vars_T = ~ismember(T.Properties.VariableNames,...
    setdiff(aux.Properties.VariableNames, GroupingVariable));

T = outerjoin(T(:, vars_T), aux, 'Keys', 'datadate',...
    'MergeKeys', true, 'Type', 'left');


InputVariables = {var_xret};
aux = varfun(@(x) mean(x, 'omitnan'), T,...
    'InputVariables', InputVariables,...
    'GroupingVariables', GroupingVariable);
aux.Properties.VariableNames(end) =...
    {'xret_idx_ew'};
xret_idx_ew = aux;

aux = removevars(aux, 'GroupCount');
vars_T = ~ismember(T.Properties.VariableNames,...
    setdiff(aux.Properties.VariableNames, GroupingVariable));

T = outerjoin(T(:, vars_T), aux, 'Keys', 'datadate',...
    'MergeKeys', true, 'Type', 'left');



InputVariables = {['w', var_retr]};
aux = varfun(@(x) mean(x, 'omitnan'), T,...
    'InputVariables', InputVariables,...
    'GroupingVariables', GroupingVariable);
aux.Properties.VariableNames(end) =...
    {'retr_idx_ew'};

aux = removevars(aux, 'GroupCount');
vars_T = ~ismember(T.Properties.VariableNames,...
    setdiff(aux.Properties.VariableNames, GroupingVariable));

T = outerjoin(T(:, vars_T), aux, 'Keys', 'datadate',...
    'MergeKeys', true, 'Type', 'left');

end