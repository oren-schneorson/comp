%{
script to replace compute_merton_vars.m

cleaning of data should be done in advance
all values should be in USD
all fundq vars should be in millions (already are...?)
mcap should be in millions

# fundq: holes sholud be filled in advance, or perhaps
I don't need to fill holes, because I'm using retime anyway

# fundq: edges cropped
# leave only unique observations
# save fundq as fundq_proc


secd_proc contains only the relevant variables

firms with missing data should be cleaned out in advance
and tracked



%}
clear; clc

project = 'comp';

matlab_dir = '/home/oren/Documents/MATLAB';
data_dir = '/media/oren/D/data';


root = fullfile(matlab_dir, project);
comp_dir = fullfile(data_dir, project);
fundq_dir = fullfile(comp_dir, 'fundq');

addpath(fullfile(matlab_dir, 'my_functions'))
addpath(fullfile(matlab_dir, 'altmany-export_fig-410f0ad'))

cd (root)



fx_daily_fpath = fullfile(data_dir, 'fx', 'fx_daily.csv');
fx_daily = readtable(fx_daily_fpath);

%%
clc

fpath = fullfile(comp_dir, 'liab_structure.xlsx');
liab_structure = readtable(fpath);
liab_structure = liab_structure(liab_structure.usgaap==1,:);
components_lctq = {'apq', 'lcoq', 'dlcq', 'txpq'};

%
fundq_vars = {'apq', 'lcoq', 'dlcq', 'txpq', 'dlttq', 'loq', 'txditcq', 'xaccq', 'dd1q', 'npq'};
%fundq_vars = {'apq', 'lcoq', 'dlcq', 'txpq', 'dlttq', 'loq', 'txditcq'};


d = dir(fundq_dir);
d = d(~[d.isdir]);
d = d(~cellfun(@(c) contains(c, '#'), {d.name}));

for i_d = 1:numel(d)

    fname = d(i_d).name;
    gvkey = strrep(fname, '.csv', '');
    
    if ~strcmp(gvkey, '001081')
        %continue
    end
    
    fpath = fullfile(fundq_dir, fname);

    fig_fpath = fullfile(root, 'gfx', 'fundq', [gvkey, '.png']);

    T = readtable(fpath);
    
    if isempty(T)
        continue
    end
    
    %{
    idx = isnan(T.lctq);
    % lct reported annualy, so will its components
    if sum(idx) > 0
        vec = sum(T{idx, components_lctq}, 2, 'omitnan');
        vec(vec==0) = NaN;
        T{idx, 'lctq'} = vec;
        idx = isnan(T.lctq);
        if sum(idx) > 0
            vec = T.ltq(idx)-sum(T{idx, {'loq', 'dlttq', 'txditcq'}}, 2, 'omitnan');
            vec(vec==0) = NaN;
            T{idx, 'lctq'} = vec;
        end
        T{:, {'ltq', 'lctq', 'apq', 'lcoq', 'dlcq', 'txpq'}} =...
            fillmissing(T{:, {'ltq', 'lctq', 'apq', 'lcoq', 'dlcq', 'txpq'}},...
            'previous');
    end
    
    idx = isnan(T.lctq);
    T{idx, 'lctq'}
    sum(idx)
    %}
    
    
    

    
    if ~all(ismember(T.curcdq, {'USD'}))
        T = innerjoin(...
            T,...
            fx_daily(:, {'datadate', 'curcdq', 'currency_val'}),...
            'Keys', {'datadate', 'curcdq'});
        T{:, fundq_vars} = T{:, fundq_vars}./T.currency_val;
    end

    idx = all(~isnan(T{:, fundq_vars}), 2);
    if sum(idx) < 2
        %continue
    end
    
    idx = true(size(T,1),1);
    clf
    subplot(2,2,1)
    area(T.datadate(idx), [T{idx, 'lctq'}, T.ltq(idx)-T{idx, 'lctq'}]./T.ltq(idx)*100)
    legend({'lctq', 'all else'}, 'Location', 'northwest') 
    title('ltq and lctq, relative')
    ylabel('pcent')
    ylim([0, 100])

    subplot(2,2,2)
    area(T.datadate(idx), [T{idx, 'lctq'}, T.ltq(idx)-T{idx, 'lctq'}])
    legend({'lctq', 'all else'}, 'Location', 'northwest') 
    title('ltq and lctq, absolute')
    ylabel('$mm')
    legend({'lctq', 'all else'}, 'Location', 'northwest') 

    subplot(2,2,3)
    area(T.datadate(idx), T{idx, components_lctq}./T.lctq(idx)*100)
    legend(components_lctq, 'Location', 'northwest') 
    title('lctq and its components, relative')
    ylabel('pcent')
    ylim([0, 100])

    subplot(2,2,4)
    area(T.datadate(idx), T{idx, components_lctq})
    legend(components_lctq, 'Location', 'northwest') 
    title('lctq and its components, relative')
    ylabel('$mm')

    export_fig(fig_fpath)

end
return

%% basic clean of csv files
clc


fundq_vars = {'acoq',...
    'actq',...
    'altoq',...
    'atq',...
    'bank',...
    'curcdq',...
    'curncdq',...
    'currtrq',...
    'lltq',...
    'lctq',...
    'ltq',...
    };

for v = fundq_vars

    v = char(v);
    fundq_dir = fullfile(comp_dir, lib,  'fundq');
    data_dir = fullfile(fundq_dir, v);

    d = dir(data_dir);
    d = d(~[d.isdir]);
    d = d(~cellfun(@(c) contains(c, '#'), {d.name}));

    loop_times = NaN(numel(d), 1);

    for i_d = 1:numel(d)
        tic
        clc
        i_d/numel(d)*100
        v
        %time_left(i_d, loop_times, true, '')
        fname = d(i_d).name;
        gvkey = strrep(fname, '.csv', '');
        fpath = fullfile(data_dir, fname);

        T = readtable(fpath);
        T = T(~isnan(T.(v)), :);
        idx = T.(v) == 0;
        T{idx, v} = NaN;
        T{idx, v} = fillmissing(T{idx, v}, 'previous');

        loop_times(i_d) = toc;
        writetable(T, fpath);
    
    
    
end

end





