function [prctiles_i, GroupMembership] = get_factor_wf_operation(i_InputVariables, T_InputVariables, idx_wf0, i_wf0, GroupMembership, params_wf)
%GET_FACTOR_WF_OPERATION Operation inside the loop of get_factor_wf
%   Conditional factoring requires computing group membership in each
%   InputVariable conditional on membership in former InputVariables. This
%   function computes group membership conditional on former
%   InputVariables. Combined with get_factor_wf, it provides overall group
%   membership. For example, grp 1: grp 11, 12, 13; grp 2: grp 21, 22, 23

    InputVariables = params_wf('InputVariables');
    bp = params_wf('bp');
    wins = params_wf('wins');
    InputVariable = InputVariables(i_InputVariables);

    % prctiles gets a (1, nptiles+1) vector of percentile values for each date
    prctiles_i=varfun(@(x)...
        prctile(x, [wins, bp{i_InputVariables}, 100-wins])-...
        [1e-10, zeros(1, numel(bp{i_InputVariables})+1)],... % adjustment for lowest entry
        T_InputVariables(idx_wf0, [{'datadate'}, InputVariable]),...
        'InputVariables', InputVariable,...
        'GroupingVariables', 'datadate');
    prctile_vname = strcat(InputVariable, num2str(i_wf0));
    prctiles_i.Properties.VariableNames(end) = prctile_vname;
    prctiles_i = removevars(prctiles_i, 'GroupCount');

    % add thresholds of group membership to T for each date
    T_InputVariables_ = outerjoin(T_InputVariables, prctiles_i,...
        'Keys', 'datadate', 'MergeKeys', true, 'Type', 'left');
    T_InputVariables_ = sortrows(T_InputVariables_, {'gvkey', 'datadate'});
    
    % index of group membership for each observation, 
    % size(aux) = [size(T_InputVariables,1), nptiles]
    aux = ...
        T_InputVariables_{:, InputVariables{i_InputVariables}} >  T_InputVariables_{:, prctile_vname}(:, 1:end-1) &...
        T_InputVariables_{:, InputVariables{i_InputVariables}} <= T_InputVariables_{:, prctile_vname}(:, 2:end);
    
    %{
    % T_InputVariables_ and T_InputVariables were not sorted in the same
    % way... must go over all joins...
    if i_InputVariables == 2
        date_ = datetime(2000,1,10);
        gvkey = 815019;
        prctile_vname
        idx_ = T_InputVariables_.gvkey == gvkey;
        idx_ = idx_ & T_InputVariables.datadate == date_;
        aux(idx_,:)
        T_InputVariables_{idx_, InputVariables{i_InputVariables}}
        T_InputVariables_{idx_, prctile_vname}(:, 1:end-1)
        T_InputVariables_{idx_, prctile_vname}(:, 2:end)
        prctiles_i{prctiles_i.datadate == date_, end}
        
        [wins, bp{i_InputVariables}, 100-wins]
        ptile = prctile(...
            T_InputVariables{...
            T_InputVariables.datadate==date_ & idx_wf0,...
            InputVariable}, [wins, bp{i_InputVariables}, 100-wins])
        ptile-[1e-10, zeros(1, numel(bp{i_InputVariables})+1)]
        error('dd')
    end
        %}
    
    
    % this makes sure I'm not overwriting previous group allocations
    % ignore non-members, conditional on previous groups
    aux(~idx_wf0) = false(sum(~idx_wf0), 1);
    
    [ia, ib] = find(aux); % group membership number (single dimension)
    GroupMembership(ia, i_InputVariables) = ib; % add result to GroupMembership




end

