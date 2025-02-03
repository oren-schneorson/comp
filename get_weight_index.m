function [T] = get_weight_index(T, InputVariable, GroupingVariable)
%GET_WEIGHT_INDEX Finds the weight of each observation
%   Uses InputVariable as the indexing variable. GroupingVariable is
%   usually date

cap = .005; % cap of weight at 0.5% of index.
%w = ['w_', char(InputVariable)];
% denominator of index, name of variable
denom = {['nansum_', char(InputVariable)]};

if ismember(InputVariable, {'equal_weight', 'ew'})
    T.ew = ones(size(T,1), 1);
end

if ischar(InputVariable)
    InputVariable = {InputVariable};
end

if ischar(GroupingVariable)
    InputVariable = {InputVariable};
end


aux = varfun(@(x) sum(x, 'omitnan'), T,...
    'InputVariables', InputVariable,...
    'GroupingVariables', GroupingVariable);    
aux = removevars(aux, 'GroupCount');
aux.Properties.VariableNames(end) = denom;

vars_T = ~ismember(T.Properties.VariableNames,...
    setdiff(aux.Properties.VariableNames, GroupingVariable));
T = join(T(:, vars_T), aux, 'Keys', GroupingVariable);

T.w = T{:, InputVariable}./T{:, denom};
T = removevars(T, denom);

if strcmp(InputVariable, 'ew')
    T = removevars(T, InputVariable);
end

% not exactly what TASE does but close
T.w = min(T.w, cap);

aux = varfun(@(x) sum(x, 'omitnan'), T,...
    'InputVariables', 'w',...
    'GroupingVariables', GroupingVariable);    
aux.Properties.VariableNames(end) = denom;
aux = removevars(aux, 'GroupCount');

% join the capped sum of weights
vars_T = ~ismember(T.Properties.VariableNames,...
    setdiff(aux.Properties.VariableNames, GroupingVariable));
T = join(T(:, vars_T), aux, 'Keys', GroupingVariable);

T.w = T.w./T{:, denom}; % adjust weight to cap
T = removevars(T, denom);



end