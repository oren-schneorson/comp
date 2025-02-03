clc

load_in_chunks = true;
fundq = table;
num_chunks = 20;
tempnames = cell(num_chunks, 1);
for chunk = 1:num_chunks
    tempnames{chunk} = tempname(fullfile(lib_comp, 'temp'));
    tempnames{chunk} = [tempnames{chunk}, '.csv'];
end


lib_fundq = fullfile(lib_comp, 'fundq');

vars_fundq = {'datadate', 'rdq', 'ltq', 'atq'};

fname = [glob, 'fundq.csv'];
fpath_fundq = fullfile(lib_comp, fname);

if exist(fpath_fundq, 'file') == 2    
    fundq = readtable(fpath_fundq);
else
    chunk_ = 1;
    chunk_size = 0;
    for gvkey_i = 1:ngvkeys

        clc
        fprintf('%.2f%%\n', gvkey_i/ngvkeys*100)
        
        chunk_size = chunk_size + d(gvkey_i).bytes;
        %chunk = ceil(chunk_size/tot_size*num_chunks);
        chunk = ceil(gvkey_i/ngvkeys*num_chunks);
        
        if load_in_chunks && chunk ~= chunk_
            writetable(fundq, tempnames{chunk_})
            fundq = table;
            chunk_size = 0;
            chunk_ = chunk;

        end


        fname = gvkeys(gvkey_i).name;
        [~,gvkey,~] = fileparts(fname);

        fpath = fullfile(lib_fundq, fname);
        fundq_i = load_fundq_all(gvkey, fx_daily, lib_fundq, var_type);
        fundq_i = fundq_i(:, vars_fundq);
        fundq_i.gvkey = repmat({gvkey}, size(fundq_i,1), 1);
        
        if isempty(fundq_i)
            continue
        end

        fundq = [fundq; fundq_i];

    end
    
    if load_in_chunks
        writetable(fundq, tempnames{num_chunks})
        clear fundq
    end
    
    if load_in_chunks
        fundq = table;
        for chunk = 1:num_chunks
            chunk/num_chunks

            fundq_ = readtable(tempnames{chunk});
            fundq = [fundq; fundq_];

        end

        % delete tempnames
        for chunk = 1:num_chunks
            delete(tempnames{chunk});
        end

        clear fundq_

    end

    idx = isnat(fundq.rdq);
    % shift reports date backwards one quarter...
    fundq.rdq(idx) = eomdate(fundq.datadate(idx)+calmonths(1));
    fundq.gvkey = cellfun(@str2double, fundq.gvkey);


    writetable(fundq, fpath_fundq)

end





