clear; clc

% set root
username = getenv("USER");
root = mfilename('fullpath');
if contains(root, 'LiveEditor') || numel(root) == 0
    % in case you're running via LiveEditor, default to root
    % change this line if running on a different system
    root = fullfile('/home', username, 'Documents/MATLAB/comp');
    cd(root)
    
else
    root = fileparts(root);
end

cd(root)


addpath(fullfile('/home', username, 'Documents', 'MATLAB', 'altmany-export_fig-410f0ad'))
addpath(fullfile('/home', username, 'Documents', 'MATLAB', 'my_functions'))


%{
Load data
%}
glob = '';

plot_gvkey = false;
vol = 'uncond_vol';
mu = 5;


lib_data = fullfile('/media', username, 'D', 'data');
comp_dir = fullfile(lib_data, 'comp');

secd_proc_dir = fullfile(comp_dir, [glob, 'secd_proc']);
merton_dir = fullfile(comp_dir, [glob, 'merton']);
merton_CDS_dir = fullfile(comp_dir, [glob, 'merton_CDS']);
merton_CDS_dir = fullfile(merton_CDS_dir,...
    sprintf('mu=%d%%__F_star&F__%s', mu, vol));

fpath = fullfile(comp_dir, 'company.csv');
company = readtable(fpath);

fpath = fullfile(comp_dir, 'gsector_to_Sector.csv');
gsector_to_Sector = readtable(fpath);

data_dir = merton_CDS_dir;

load handel.mat;
end_sound = y;
end_sound_Fs = Fs;


%{
******************************************
LOAD COMPANY MERTON INFORMATION BY GVKEY
%}

d = dir(data_dir);
d = d(~[d.isdir]);
d = d(~[cellfun(@(c) contains(c, '~'), {d.name})]);
d = char({d.name}');


d = d(:, 1:end-4);
gvks = cellstr(d);
clear d

T = table;

for gvk_i = 1:numel(gvks)
    
    clc
    fprintf('%.2f%%\n', gvk_i/numel(gvks)*100)
    gvkey = gvks{gvk_i};
    fname = [gvkey, '.csv'];
    fpath = fullfile(data_dir, fname);
    T_ = readtable(fpath);
    T_ = T_(~isnan(T_.V),:);
    T_ = T_(~isinf(T_.DD),:);
    
    T_.gvkey = str2double(gvkey)*ones(size(T_,1), 1);
    T_ = T_(:, [{'gvkey', 'datadate'}, T_.Properties.VariableNames(2:end-1)]);
    T = [T; T_];
    
    if plot_gvkey
        clf
        set(gcf, 'Color', 'w')
        
        company_i = company(company.gvkey ==  str2double(gvkey),:);
        
        hold on
        plot(T_.datadate, T_.mcap)
        plot(T_.datadate, T_.E)
                
        title(company_i.conm)
        leg = {'mcap', 'E^*'};
        legend(leg)
        ylabel('mm $')
        fname = [gvkey, '.png'];
        fpath = fullfile(root, 'gfx', 'merton_CDS', fname);
        export_fig(fpath)
        
        
        
        
    end
    
end

clear T_


%% add variables
%{
******************************************
ADD VARIABLES FROM company.csv
(with values varying only in the cross section)
%}

clc
vars = {'gsector'};
%vars = {'mcap'};
min_comp_per_day = 10;

if ismember(vars, {'mcap'})
    clc
    num_grps = 3;
    aux = (0:num_grps)/num_grps*100;
    unq_dates = unique(T.datadate);
    prctiles = NaN(size(unq_dates, 1), num_grps+1);
    prctile_grp = NaN(size(T, 1), 1);
    for t = 1:numel(unq_dates)
        clc
        fprintf('Computing group membership by size, daily, %.2f%%',...
            t/numel(unq_dates)*100)
        idx = T.datadate == unq_dates(t);
        prctiles(t, :) = prctile(T.mcap(idx), aux);
        %prctiles(t, :)
        grps = true(sum(idx), size(prctiles, 2));
        grps = grps & T.mcap(idx) >= prctiles(t, :);
        [ia, ib] = find(grps);
        i = [ia, ib];
        i = sortrows(i, [1,2], {'ascend', 'descend'});
        [~, idx_i] = unique(i(:,1));
        prctile_grp(idx) = i(idx_i, 2);
        
    end
    
    
    prctile_grp = array2table(prctile_grp);
    prctile_grp = renamevars(prctile_grp, 'prctile_grp', 'size');
    prctile_grp.gvkey = T.gvkey;
    prctile_grp.datadate = T.datadate;
    
    
    
    T_ = outerjoin(T,...
    prctile_grp,...
    'Type', 'left',...
    'MergeKeys', true,...
    'Keys', {'gvkey', 'datadate'});

%{
    T_ = outerjoin(T,...
    prctile_grp(:, ['gvkey', vars]),...
    'Type', 'left',...
    'MergeKeys', true,...
    'Keys', 'gvkey');
%}

    
    
    
    %T_.size = prctile_grp;
    %vars = {'size'};
    
else
idx = ~isnan(company{:, vars});
idx = any(idx, 2);
unq_gkvey = unique(company.gvkey(idx));

head(company)
%{
% all firms should have gsector
T_ = innerjoin(T(ismember(T.gvkey, unq_gkvey),:),...
    company(idx, {'gvkey', var}),...
    'Keys', 'gvkey');
%}
T_ = outerjoin(T,...
    company(idx, ['gvkey', vars]),...
    'Type', 'left',...
    'MergeKeys', true,...
    'Keys', 'gvkey');

end


%%
head(T_)

%%
%{
T_ = outerjoin(T,...
prctile_grp,...
'Type', 'left',...
'MergeKeys', true,...
'Keys', {'gvkey', 'datadate'});
%}

T_ = outerjoin(T,...
    company(:, ['gvkey', vars]),...
    'Type', 'left',...
    'MergeKeys', true,...
    'Keys', 'gvkey');

%vars = {'size'};

GroupingVariables = [{'datadate'}, vars];

sum_mcap = varfun(@(mcap) sum(mcap, 'omitnan'),...
    T_,...
    'InputVariables', 'mcap',...
    'GroupingVariables', GroupingVariables);
sum_mcap.Properties.VariableNames{end} = 'sum_mcap';



idx = sum_mcap.GroupCount >= min_comp_per_day;
sum_mcap = sum_mcap(idx,:);
unq_dates = sum_mcap.datadate;

T_ = T_(ismember(T_.datadate, unq_dates), :);




T_ = innerjoin(T_, sum_mcap, 'Keys', GroupingVariables);

T_.w = T_.mcap./T_.sum_mcap;
T_.ew = 1./T_.GroupCount;
T_ = removevars(T_, 'sum_mcap');

T_.ratio_F_star_F = log(T_.F_star./T_.F);
T_.ratio_F_star_F_normalized = T_.ratio_F_star_F./T_.sigma_V;


T_.w_DD = T_.w .* T_.DD;
T_.w_ratio_F_star_F = T_.w .* T_.ratio_F_star_F;
T_.w_ratio_F_star_F_normalized = T_.w .* T_.ratio_F_star_F_normalized;

T_.ew_DD = T_.ew .* T_.DD;
T_.ew_ratio_F_star_F = T_.ew .* T_.ratio_F_star_F;
T_.ew_ratio_F_star_F_normalized = T_.ew .* T_.ratio_F_star_F_normalized;
%%
w = 'w';
var = 'ratio_F_star_F';
var = [w, '_', var];

index_ = varfun(@(x) sum(x, 'omitnan'),...
    T_,...
    'InputVariables', var,...
    'GroupingVariables', GroupingVariables);
index_.Properties.VariableNames{end} = 'w';


index = index_;

w = 'ew';
var = 'ratio_F_star_F';
var = [w, '_', var];

index_ = varfun(@(x) sum(x, 'omitnan'),...
    T_,...
    'InputVariables', var,...
    'GroupingVariables', GroupingVariables);
index_.Properties.VariableNames{end} = 'ew';

index = join(index, index_, 'Keys', [GroupingVariables, {'GroupCount'}]);
clear index_

head(index)

clf
start_date = datetime(1999,1,1);
end_date = datetime(2023,5,28);

subplot(1,1,1)
leg = {};
[g, g_] = findgroups(index.gsector);
%for i_g=4:8%:size(g_,1)
%for i_g=[9, 10, 11, 7]
for i_g=1:size(g_,1)
    %gsector_to_Sector(gsector_to_Sector.gsector == gsector, :)
    gsector = g_(i_g);
    
    idx = true(size(index, 1), 1);
    idx = idx & index.datadate >= start_date;
    idx = idx & index.datadate <= end_date;
    idx = idx & g == i_g;
    val = g_(i_g);
    leg = [leg,...
        gsector_to_Sector.GIC(...
        gsector_to_Sector.gsector == val)];
    subplot(2,2,1)
    hold on
    a = scatter(index.datadate(idx), index.w(idx));
    a.SizeData = 6;
    a.Marker = '.';
    
    subplot(2,2,2)
    hold on
    idx = true(size(sum_mcap, 1), 1);
    idx = idx & sum_mcap.datadate >= start_date;
    idx = idx & sum_mcap.datadate <= end_date;
    idx = idx & sum_mcap.gsector == gsector;
    %idx = idx & sum_mcap.size == val;
    plot(sum_mcap.datadate(idx), sum_mcap.sum_mcap(idx))
    
    subplot(2,2,3)
    hold on
    plot(index.datadate(idx), index.GroupCount(idx))
    
    
end


var = strrep(var, '_', '\_');

subplot(2,2,1)
title('Merton CDS', var)
legend(leg)

subplot(2,2,2)
title('mcap')
legend(leg)

subplot(2,2,3)
title('GroupCount')
legend(leg)
% sprintf('by %s', GroupingVariables{2:end})
%leg = {'value weighted', 'equally weighted'};

%%
clc

%head(statarray)

%summary(DD_index)

InputVariables = {'mcap', 'V', 'sigma_V', 'DD', 'F', 'F_star', 'ratio_F_star_F'};
GroupBy = {'gvkey'};
%GroupBy = {'datadate'};
stats2compute = {'mean','median','mode','std','min','max','skewness','kurtosis'};
idx = T_.datadate > datetime(1999,1,1); 



for InputVariable = InputVariables

    InputVariable
    statarray = grpstats(T_(idx, [InputVariable, GroupBy]),GroupBy, stats2compute);

    fname = sprintf('merton_CDS_GroupBy=%s_InputVariable=%s.csv', char(GroupBy), char(InputVariable));
    fpath = fullfile(root, 'desc_stats', fname);
    writetable(statarray, fpath)
    
    
    for v=statarray.Properties.VariableNames(2:end)
        fname = sprintf('merton_CDS_GroupBy=%s_InputVariable=%s.png', char(GroupBy), char(v));
        fpath = fullfile(root, 'desc_stats', 'gfx', fname);
        
        clf
        if ismember(GroupBy, {'gvkey'})
            histogram(statarray{:, v},'Normalization','probability')
        elseif ismember(GroupBy, {'datadate'})
            plot(statarray.datadate, statarray{:, v})
        end
        % sum(h.Values) == 1
        
        tit = [char(v), sprintf('_GroupBy=%s', char(GroupBy))];
        tit = strrep(tit, '_', '\_');
        title(tit)
        set(gcf, 'Color', 'w')
        export_fig(fpath)
        
        
        
        
    end
    
    
    
end






















