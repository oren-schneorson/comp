clc



% for printing factor names in figures
factors_nm = cellfun(@(nm) strrep(nm, '_', '\_'),...
    factors, 'UniformOutput', false);

stage1 = table;
% loop on cross-section dimension
params('ngvkeys') = ngvkeys;
redo_all = false;

min_gvkey = 1;
%min_gvkey = round(ngvkeys/2);

for gvkey_i = min_gvkey:ngvkeys
%for gvkey_i = flip(1:ngvkeys)
%for gvkey_i = flip(1:min_gvkey)
    
    clc
    fprintf('Fama McBeth: stage 1; %.2f%%', gvkey_i/ngvkeys*100)
    params('gvkey_i') = gvkey_i;

    gvkey = gvkeys(gvkey_i).name;

    fpath_sec_proc = fullfile(lib_gvkey_proc, gvkey);
    fpath_stage1_i = fullfile(root, 'stage1', 'mw', freq, gvkeys(gvkey_i).name);

    [~, gvkey, ~] = fileparts(char(gvkey));
    gvkey = str2double(gvkey);
        
    if exist(fpath_stage1_i, 'file') == 2 && ~redo_all
        continue
    end

    if ~isempty(T)
        idx = T.gvkey == gvkey;    
        sec_proc_i = T(idx, :); % firm-level data
    elseif exist(fpath_sec_proc, 'file') == 2
        sec_proc_i = readtable(fpath_sec_proc);
    else
        continue
    end
    
    if any(~ismember(sec_proc_i.exchg, [11, 12, 14]))
        continue
    end

    sec_proc_i = sortrows(sec_proc_i, 'datadate', 'ascend'); % sorted by date
    [~, idx] = unique(sec_proc_i.datadate);
    sec_proc_i = sec_proc_i(idx, :);
        
    if isempty(sec_proc_i)
        continue
    end

    sec_proc_i = outerjoin(sec_proc_i, R_f_,...
        'Keys', 'datadate', 'Type', 'left', 'MergeKeys', true);
    
    if strcmp(freq, 'D')
        sec_proc_i.xretd = sec_proc_i.retd-sec_proc_i.(...
            R_f_.Properties.VariableNames{end});
    elseif strcmp(freq, 'M')
        sec_proc_i.xretm = sec_proc_i.retm-sec_proc_i.(...
            R_f_.Properties.VariableNames{end});
    else
        error('freq must be D or M.')
    end
 
    
    for factor = factors

        factor = char(factor);
        Keys = {'datadate'}; % if more group specific, need to specify other Keys
        COMMAND = sprintf("sec_proc_i = outerjoin(sec_proc_i, %s, 'Type', 'left', 'Keys', Keys, 'MergeKeys', true);",...
            factors_tab(factor));
        eval(COMMAND)

        sec_proc_i = sec_proc_i(:, ~ismember(sec_proc_i.Properties.VariableNames, 'GroupCount'));

    end
    
    stage1_i = stage1_reg_mw(sec_proc_i, params);
    writetable(stage1_i, fpath_stage1_i)
    
end

%{
*********************************
% write table to spreadsheet file
%}
%save(fpath_stage1, 'stage1', '-v7.3');

return
%%
clc
stage1.date = stage1.max_date;

% join beta to T 
idx_T = contains(T.Properties.VariableNames, 'beta');
idx_stage1 = contains(stage1.Properties.VariableNames, 'beta');
idx_stage1 = idx_stage1 & ~contains(stage1.Properties.VariableNames, 't_');
idx_stage1 = idx_stage1 & ~contains(stage1.Properties.VariableNames, '_int');

vars_T = T.Properties.VariableNames(~idx_T);
vars_stage1 = stage1.Properties.VariableNames(idx_stage1);
Keys = {'iid', 'date'};
T = outerjoin(T(:, vars_T), stage1(:, [Keys, vars_stage1]),...
    'Keys', Keys, 'Type', 'left',...
    'MergeKeys', true);

T = renamevars(T, 'beta_1', 'beta_Epi');







return
%% plot charts of beta estimation, conditional beta
%{
2023-05-08
Results are unsatisfactory. The factor of inflation expectations
explains little variation in returns (single factor model). betas
are not significantly different from zero.

Need to check market factor and see if this is something that
occurs with standard factors as well.

%}


clc
if true
    [g, g_] = findgroups(stage1.iid);
    for i_g = 771:size(g_,1)
        clf
        subplot(1,2,1)
        hold on        
        idx = g == i_g;
        gvkey = g_(i_g);
        fname = [num2str(gvkey), '.png'];
        
        fpath = fullfile(...
            '/home/oren/Documents/MATLAB/TASE/stage1/gfx',...
            fname);

        plot(stage1.max_date(idx), stage1.beta_1(idx), 'b-')
        %plot(stage1.max_date(idx), stage1.beta_1_int1(idx), 'b--')
        %plot(stage1.max_date(idx), stage1.beta_1_int2(idx), 'b--')
        title(gvkey)
        
        a = [stage1.beta_1_int1(idx), stage1.beta_1_int2(idx)-stage1.beta_1_int1(idx)];
        a = area(stage1.max_date(idx), a);
        a(1).LineStyle = 'none';
        a(2).LineStyle = 'none';
        a(1).FaceColor = 'none';
        a(2).FaceColor = 'b';
        a(1).FaceAlpha = 0;
        a(2).FaceAlpha = .3;
        
        legend({'estimate of \beta', '', '95%  bounds'})
        subplot(1,2,2)
        plot(stage1.max_date(idx), stage1.R_sq(idx), 'b-')
        title('R^2')
        
        
        export_fig(fpath)
        
        
    end
    
    
end



%%


%{
if ~monthly && dates_CBS
    freq = 252;
else
    freq = 12;
end

hh = freq*5; % num busdays a year times 5 years
h = log(2)/freq;

weight_func = @(t_, tau_) exp(-abs(t_-tau_).*h)./sum(...
    exp(-abs(t_-[0:t_-2]).*h));

%{
% moving window
% check what weight_func actually outputs
% (for daily freq!)
clc; clf

if ~monthly && dates_CBS
    freq = 252;
else
    freq = 12;
end

hh = freq*5; % num busdays a year times 5 years
h = log(2)/freq;
t=300;
tau=0:t-1;

% not excluding time t
weight_func = @(t_, tau_) exp(-abs(t_-tau_).*h)./sum(...
    exp(-abs(t_-[1:t_-1]).*h));


res = weight_func(t, tau);
sum(weight_func(t, tau))
plot(tau, res./res(end))
title(t)

%weight_func(t, t)/weight_func(t, t-hh)
%}

%}

















