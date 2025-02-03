
%{
preamble, part of a package for TASE.

1. main.m
2. preamble.m
3. load_data.m
4. merge2T.m
5. get_indices.m

For further details see main.m

Last quality check: 2 Jan 2023
TODO: there's a lot of overlap in metadata. Unite to one single file with
dictionary for heb
%}

% for using export_fig function to produce charts
addpath('/home/oren/Documents/MATLAB/altmany-export_fig-410f0ad')
addpath('/home/oren/Documents/MATLAB/my_functions')

FontName = 'Utopia'; % Font


gvkey_typs = {'Ordinary Share'}; % Ordinary Share by default
fname = 'equity';
% Note: change fname accordingly!

lib_data = '/media/oren/D/data'; % main directory of SM data
lib_comp = fullfile(lib_data, 'comp'); % main directory of SM data

if strcmp(freq, 'D')
    lib_gvkey = fullfile(lib_comp, 'secd'); % individual securities
    lib_gvkey_proc = fullfile(lib_comp, 'secd_proc'); % individual securities
elseif strcmp(freq, 'M')
    lib_gvkey = fullfile(lib_comp, 'secm'); % individual securities
    lib_gvkey_proc = fullfile(lib_comp, 'secm_proc'); % individual securities
else
    error('freq must be D or M.')
end

added_str_ = added_str('', {}, min_date, max_date, freq, monthly_avg, '', '', dates_CBS);
fname_full = [fname, added_str_, '.mat'];
fname = [fname, added_str_, '.csv'];

fpath_equity = fullfile(lib_comp, fname);
fpath_equity_full = fullfile(lib_comp, fname_full);

% note: for many applications, equity is needed, so I saved it to save time
% instead of loading it "from source" (csv file by gvkey).


%{
***************************
***** METADATA
***************************
Load metadata and auxiliary tables:

2. sectors: table connecting containing company level metadata

3. security level information
3.1. securitiesinfo: table containing security level metadata (listed)
3.2. delistedsecurities: table containing security level metadata (delisted)
3.3. companydelistedsecurities table containing security level metadata
with additional information on the company
(delisted)
3.4 attributes: table containing security level metadata (hebrew, listed)
?. iid_tase_detail: table containing metadata at security level (incl cid)
4. dates_CBS_inflation: table containing information about dates of the CBS
published inflation figures.
5. var_desc: dictionary of variables


Last quality check: 2 Jan 2023

%}


%{
***************************
***** company metadata (identifier: cid / corpn)

** sectors:

% identifiers and name
cid
corpn (corporation id, if company incorp abroad, the id abroad...)
Bloomberg (Bloomberg key for company/stock)
conm_heb (company name)
conml_heb (long company name)

% sector information
level%d: TASE sectoral categroizatin
level0
level1
level2
level3

% [compustat equivalent of sector information]
gsubind 
gsubind_name
gind
gind_name
ggroup
ggroup_name
gsector
gsector_name

% company contact details
Address
Phone
Fax
CompanyHomepage
E-Mail


Last quality check: 2 Jan 2023

** companydelistedsecurities_eng:

%}



%{
***************************
***** security metadata (identifier: gvkey / isin)

** securitiesinfo:
iid (TASE identifier)
isin (international identifier)
iidnm_heb
iidnm (name of security)
sym_heb
sym (trading symbol)
iid_typ_heb
iid_typ (Corporate bond, Ordinary Share...)
level1_heb
level1 (same as sectors)
level2_heb
level2 (same as sectors)
level3_heb
level3 (same as sectors)
idx_incl_heb
idx_incl (indices the iid is included in)
corpn (corporation number -- some missing values here)
cid (one missing obs, 1092709)
conm_heb
conm (company name, potential differences from sectors)
incorp_heb
incorp (country of incorporation)

Last quality check: 2 Jan 2023


** delistedsecurities:
data is from TASE delisted securities data.
English and Hebrew are different. "all" Sheet joins the two

iid
isin
iid_typ
iidnm_heb
sym_heb
iid_typ_heb
cid (many obs missing cid. A check in KOGNOS did not yield adequate results)
iidnm
sym

%}

run('var_types.m')
fpath = fullfile(lib_comp, 'company.csv');
company = readtable(fpath);
fpath = fullfile(lib_comp, 'security.csv');
security = readtable(fpath);


gvkeys = dir(lib_gvkey_proc);
gvkeys = gvkeys(~[gvkeys.isdir]);
gvkeys = gvkeys(~cellfun(@(c) contains(c, '#'), {gvkeys.name}));


ngvkeys = numel(gvkeys); % number of securities in data

Holidays = readtable('../holidays.xlsx');
%Closures = readtable('../nyseclosures.xlsx');

holidays = Holidays.holidays;
clear Holidays

 

