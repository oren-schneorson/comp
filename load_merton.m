clear; clc


%{
Load data
%}
glob = '';

addpath('/home/oren/Documents/MATLAB/altmany-export_fig-410f0ad')
addpath('/home/oren/Documents/MATLAB/my_functions')
root = '/home/oren/Documents/MATLAB/comp';
cd(root)

comp_dir = '/media/oren/D/data/comp';
secd_proc_dir = fullfile(comp_dir, [glob, 'secd_proc']);
merton_dir = fullfile(comp_dir, [glob, 'merton']);
merton_CDS_dir = fullfile(comp_dir, [glob, 'merton_CDS']);


%data_dir = merton_dir;

load handel.mat;
end_sound = y;
end_sound_Fs = Fs;

% get gvkeys from merton_CDS_dir
d = dir(merton_CDS_dir);
d = d(~[d.isdir]);
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
    fpath = fullfile(merton_dir, fname);
    T_ = readtable(fpath);
    idx = T_.DD > prctile(T_.DD, 1) & T_.DD < prctile(T_.DD, 99);
    T_ = T_(idx,:);
    T_ = T_(~isnan(T_.DD),:);
    T_ = T_(22:end, :);
    
    T_.gvkey = str2double(gvkey)*ones(size(T_,1), 1);
    T_ = T_(:, [{'gvkey', 'datadate'}, T_.Properties.VariableNames(2:end-1)]);
    
    
    T = [T; T_];
    
end

%%

sum_mcap = varfun(@(mcap) sum(mcap, 'omitnan'),...
    T,...
    'InputVariables', 'mcap',...
    'GroupingVariables', 'datadate');
sum_mcap.Properties.VariableNames{end} = 'sum_mcap';

idx = sum_mcap.GroupCount >= 500;
unq_dates = sum_mcap.datadate(idx);
T_ = T(ismember(T.datadate, unq_dates), :);
sum_mcap = sum_mcap(idx,:);




T_ = innerjoin(T_, sum_mcap, 'Keys', 'datadate');

%%

T_.weight = T_.mcap./T_.sum_mcap;
T_.weighted_DD = T_.weight .* T_.DD;


DD_index = varfun(@(weighted_DD) sum(weighted_DD, 'omitnan'),...
    T_,...
    'InputVariables', 'weighted_DD',...
    'GroupingVariables', 'datadate');
DD_index.Properties.VariableNames{end} = 'DD_index';


%%
clc
head(T_)


idx = DD_index.datadate > datetime(1999,1,1); 

subplot(1,2,1)
hold on
plot(DD_index.datadate(idx), DD_index.DD_index(idx))
%%
title('Comparing Merton DD and DD from CDS')
legend({'Merton\_CDS', 'Merton, \mu=0'})

























