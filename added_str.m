function [added_str_] = added_str(model_, factors, min_date, max_date, freq, monthly_avg, mw, w, dates_CBS)
%ADDED_STR Standard added string for robustness checks
%   Detailed explanation goes here


if dates_CBS
    dates_CBS = '_dates_CBS';
else
    dates_CBS = '';
end

if monthly_avg
    monthly_avg = '_avg';
else
    monthly_avg = '';
end

w = strrep(w, '_', '');

factors = cellstr(factors);
factors = strjoin(factors, '_');
date_format = 'yyyy-mm';

min_date = datestr(min_date, date_format);
max_date = datestr(max_date, date_format);

added_str_ = sprintf('%s_%s_%s_%s - %s__%s%s%s%s',...
    model_,...
    char(factors),...
    freq,...
    min_date,...
    max_date,...
    w,...
    mw,...
    dates_CBS,...
    monthly_avg);

added_str_ = strip(added_str_);

end

