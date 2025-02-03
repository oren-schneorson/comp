

clc
if strcmp(freq, 'D')

    fprintf('Loading risk free rate, 3M TBill, daily, 1954-2023\n')
    series = 'DTB3';
    lib_R_f = fullfile(lib_data, 'ir', 'DTB');
    [R_f, metadata_R_f] = load_data_files(series, lib_R_f, 'FRED');

    R_f.DTB3 = fillmissing(R_f.DTB3, 'previous');

    R_f = table2timetable(R_f);
    R_f = retime(R_f, 'daily', 'previous');
    R_f = timetable2table(R_f);

    R_f = renamevars(R_f, 'date', 'datadate');
elseif strcmp(freq, 'M')
    fprintf('Loading risk free rate, 1M TBill, monthly, 1926-2023\n')
    series = 'RF_FF';
    lib_R_f = fullfile(lib_data, 'ir');
    [R_f, metadata_R_f] = load_data_files(series, lib_R_f, 'FRED');
    % annualize
    R_f.RF_FF = (R_f.RF_FF/100+1).^12;
    R_f.RF_FF = R_f.RF_FF-1;
    R_f.RF_FF = R_f.RF_FF*100;

    R_f = renamevars(R_f, 'date', 'datadate');
else
    error('freq must be D or M.')    
end