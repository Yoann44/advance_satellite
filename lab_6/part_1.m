close all
clear

% Load rover observations
load('datar.mat')
% Load master observations
load('datam.mat')

prn_number = 24;

min_tow = min(datam(:,1));
max_tow = max(datam(:,1));

np_of_sv = 8;

frequ_1 = 1575.42e6;
frequ_2 = 1227.6e6;
frequ_5 = 1176.45e6;

obs_master = zeros((max_tow - min_tow) * 4 * (np_of_sv - 1), 1);
obs_rover = zeros((max_tow - min_tow) * 4 * (np_of_sv - 1), 1);

i = 0;
for tow = min_tow:max_tow
    i = i + 1;
    
    ref_sat_m = select_prn(datam_tow, prn_number);
    ref_sat_r = select_prn(datar_tow, prn_number);
    
    datam_tow = select_tow(datam, tow);
    datam_tow = datam_tow(datam_tow(:, 2) ~= prn_number, :);
    
    datar_tow = select_tow(datar, tow);
    datar_tow = datar_tow(datar_tow(:, 2) ~= prn_number, :);
    
    j = 1;
    l(1) = ref_sat_m(3:6);
    l(2) = ref_sat_r(3:6);
    for data_prn_tow = datam_tow
        j = j + 2;
        
        l(j) = ref_sat_m(3:6);
        l(j+1) = ref_sat_r(3:6);
    end
end

function [data_tow] = select_tow(data, tow)
    data_tow = data(data(:, 1) == tow, :);
end

function [data_tow] = select_prn(data, prn)
    data_tow = data(data(:, 2) == prn, :);
end



