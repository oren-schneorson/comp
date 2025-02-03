clc

% sample
params('min_date') = min_date;
params('max_date') = max_date;

% Fama MacBeth: stage 1
lead = 0; % lead X on Y, 0 by default.

% for printing factor names in figures
factors_nm = cellfun(@(nm) strrep(nm, '_', '\_'),...
    factors, 'UniformOutput', false);

redo_stage1 = false;

if exist('stage1', 'var') ~= 1 || redo_stage1
    stage1 = table; % declare table to collect results
end

% loop on cross-section dimension
for gvkey_i = 1:ngvkeys
        
    %{
    txt = ['Fama McBeth: stage 1; factors: ', sprintf('%s, ', factors{:}),...
        sprintf('\n%.2f%%', i/ngvkeys*100), '\n'];
    %}
    clc
    fprintf('Fama McBeth: stage 1; %.2f%%', gvkey_i/ngvkeys*100)
    
    gvkey = gvkeys(gvkey_i).name;
    fpath_i = fullfile(lib_gvkey_proc, gvkey);
    [~, gvkey, ~] = fileparts(char(gvkey));
    gvkey = str2double(gvkey);
    %{
    if any(stage1.gvkey == gvkey)
        continue
    end
    %}
    
    if ~isempty(T)
    
        idx = ismember(T.gvkey,  gvkey);
        T_ = T(idx, :);
    elseif exist(fpath_i, 'file') == 2
        T_ = readtable(fpath_i);
        T_.retm = (T_.retm-1)*100;
    else
        continue
    end
     
    [~, idx] = unique(T_.datadate);
    T_ = T_(idx, :);
    
    if isempty(T_)
        continue
    end
    
    for factor = factors

        factor = char(factor);
        Keys = {'datadate'}; % if more group specific, need to specify other Keys
        COMMAND = sprintf("T_ = outerjoin(T_, %s(:, [Keys, {factor}]), 'Type', 'left', 'Keys', Keys, 'MergeKeys', true);",...
            factors_tab(factor));
        eval(COMMAND)

        T_ = T_(:, ~ismember(T_.Properties.VariableNames, 'GroupCount'));

    end

    stage1_i = stage1_reg(T_, params);    
    stage1_i.gvkey = gvkey;
    stage1 = [stage1; stage1_i];
    clear stage1_

    

end

%% write table to spreadsheet file
writetable(stage1, fpath_stage1) 


return
%%
clc
head(sortrows(stage1, 'R_sq', 'ascend'))
size(stage1)
clf
histogram(stage1.beta_1)

%%
% join beta to T
idx_T = contains(T.Properties.VariableNames, 'beta');
idx_stage1 = contains(stage1.Properties.VariableNames, 'beta');
idx_stage1 = idx_stage1 & ~contains(stage1.Properties.VariableNames, 't_');
idx_stage1 = idx_stage1 & ~contains(stage1.Properties.VariableNames, '_int');
vars_T = T.Properties.VariableNames(~idx_T);
vars_stage1 = stage1.Properties.VariableNames(idx_stage1);
Keys = {'gvkey'};
T = innerjoin(T(:, vars_T), stage1(:, [Keys, vars_stage1]), 'Keys', Keys);
%}

%T = renamevars(T, 'beta_1', 'beta_Epi');
