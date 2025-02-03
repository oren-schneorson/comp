%{
load_iids.m, part of a package for TASE, sub-package, load_data

1. main.m
2. preamble.m
3. load_data.m --> load_iids
4. merge2T.m
5. get_indices.m

%}

clc
load_in_chunks = true;

% list all data files in lib_gvkeys
d = dir(lib_gvkey_proc);
d = d(~[d.isdir]);
d = d(~cellfun(@(c) contains(c, '#'), {d.name}));

tot_size = sum([d.bytes]);

T = table;
num_chunks = 20;
tempnames = cell(num_chunks, 1);
for chunk = 1:num_chunks
    tempnames{chunk} = tempname(fullfile(lib_comp, 'temp'));
    tempnames{chunk} = [tempnames{chunk}, '.csv'];
end

if strcmp(freq, 'D')
    vars_req = {'gvkey', 'datadate', 'retd', 'xretd', 'retrd', 'mcap', 'atq', 'ltq', 'gsector', 'ggroup', 'gind'};
elseif strcmp(freq, 'M')
    vars_req = {'gvkey', 'datadate', 'retm', 'xretm', 'retrm', 'mcap', 'atq', 'ltq', 'gsector', 'ggroup', 'gind'};
else
    error('freq must be D or M.')
end
    

aux = cellfun(@(c) strrep(c, '.csv', ''), {d.name},...
    'UniformOutput', false);
aux = cellfun(@str2double, aux,...
    'UniformOutput', false);
[d.('gvkey')] = aux{:};
date_format = 'yyyy-MM-dd';


R_f_ = R_f;
R_f_.DTB3 = R_f_.DTB3/100+1;
R_f_.DTB3 = R_f_.DTB3.^(1/252);
R_f_.DTB3 = (R_f_.DTB3-1)*100;

ls = [];
%%
chunk_ = 1;
chunk_size = 0;

for gvkey_i = 1:ngvkeys

    clc
    fprintf('%d: Loading data, %.2f%%\n', gvkey_i, gvkey_i/ngvkeys*100)
    
    chunk_size = chunk_size + d(gvkey_i).bytes;
    %chunk = ceil(chunk_size/tot_size*num_chunks);
    chunk = ceil(gvkey_i/ngvkeys*num_chunks);

    if load_in_chunks && chunk ~= chunk_
        writetable(T, tempnames{chunk_})
        T = table;
        chunk_size = 0;
        chunk_ = chunk;
        
    end

    
    fname = char(d(gvkey_i).name);
    gvkey = str2double(strrep(fname, '.csv', ''));
    % get cid, gvkey_typ, sector, subsector etc.
    
    fpath = fullfile(lib_gvkey_proc, fname);
    idx_company = company.gvkey == gvkey;
    %company(idx_company,:)
    gsector = company.gsector(idx_company);
    ggroup = company.ggroup(idx_company);
    gind = company.gind(idx_company);
        
    % skip empty files
    if d(gvkey_i).bytes <= 566
        continue
    end
        
    opts = gen_opts(fpath, 'datadate', date_format, var_type);
    T_ = readtable(fpath, opts);
    
    [~, idx] = unique(T_.datadate);
    T_ = T_(idx, :);
    
    T_ = sortrows(T_, 'datadate', 'ascend');
    % join CPI
    Keys = {'datadate'};
    T_ = outerjoin(T_, CPI_D(:, [Keys, {'CPI_SA', 'PI_SA'}]),...
        'Keys', Keys, 'Type', 'left', 'MergeKeys', true);
    T_ = outerjoin(T_, R_f_(:, [Keys, {'DTB3'}]),...
        'Keys', Keys, 'Type', 'left', 'MergeKeys', true);

    if isempty(T_)
        % drop iid with data only after last obs of CPI
        continue
    end
    
    T_min_date = min(T_.datadate);
    T_max_date = max(T_.datadate);
    
    if strcmp(freq, 'D')
        T_.retrd = [NaN; T_.retd(2:end)./(T_.CPI_SA(2:end)./T_.CPI_SA(1:end-1))];
        T_.retd = (T_.retd-1)*100;
        T_.retrd = (T_.retrd-1)*100;
        T_.xretd = T_.retd-T_.DTB3;
    elseif strcmp(freq, 'M')
        T_.retrm = [NaN; T_.retm(2:end)./(T_.CPI_SA(2:end)./T_.CPI_SA(1:end-1))];
        T_.retm = (T_.retm-1)*100;
        T_.retrm = (T_.retrm-1)*100;
        T_.xretm = T_.retm-T_.DTB3;
    end
    
    T_.gvkey = repmat(gvkey, size(T_,1), 1);
    T_.gsector = repmat(gsector, size(T_,1), 1);
    T_.ggroup = repmat(ggroup, size(T_,1), 1);
    T_.gind = repmat(gind, size(T_,1), 1);
        

    % reorder variables
    T_ = T_(:, vars_req);
        
    if ~isempty(T)
        T = [T; T_];
    else
        T = T_;
    end
    
    
       
end

if load_in_chunks
    writetable(T, tempnames{num_chunks})
    clear T

end


%%

if load_in_chunks
    T = table;
    for chunk = 1:num_chunks
        chunk/num_chunks
       
        T_ = readtable(tempnames{chunk});
        T = [T; T_];
        
    end
    
    for chunk = 1:num_chunks

        delete(tempnames{chunk});

    end

    clear T_
    
end




