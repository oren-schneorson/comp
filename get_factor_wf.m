function [prctiles, GroupMembership] = get_factor_wf(i_InputVariables, idx_wf0, i_wf0, T_InputVariables, prctiles, GroupMembership, params_wf)
%GET_FACTOR_WF Dynamic loop for get_factor_cond
%   Finds group membership conditional on former group membership.


InputVariables = params_wf('InputVariables');
[prctiles_i, GroupMembership] = get_factor_wf_operation(i_InputVariables, T_InputVariables, idx_wf0, i_wf0, GroupMembership, params_wf);
[grp_wf, grps_wf] = findgroups(GroupMembership(:, i_InputVariables));

if isempty(prctiles)
    prctiles = prctiles_i;
else
    prctiles = outerjoin(prctiles, prctiles_i,...
        'Keys', 'datadate', 'MergeKeys', true);
end

prctiles = sortrows(prctiles, 'datadate');

if i_InputVariables >= numel(InputVariables)
    return
end

for i_wf1 = 1:numel(grps_wf)
    % group membership conditional on former groups
    idx_wf1 = idx_wf0 & grp_wf == i_wf1;
    [prctiles, GroupMembership] = get_factor_wf(i_InputVariables+1, idx_wf1, i_wf1, T_InputVariables, prctiles, GroupMembership, params_wf);
end






end

