close all;
clear;
clc;

% Load rover observations
load('datar.mat')
% Load master observations
load('datam.mat')

%prn_number = 24;
prn_number = 8;

min_tow = min(datam(:,1));
max_tow = max(datam(:,1));

np_of_sv = 8;

%GPS band frequencies [Hz]
frequ_1 = 1575.42e6;
frequ_2 = 1227.6e6;
frequ_5 = 1176.45e6;

%Deviation [m] for [C1 C5 L1 L5]
std = [0.5 0.5 0.01 .01];

%Speed of light [m/s]
c = 299792458;

%Convert cycle to meters in L observations
datam(:,5) = datam(:,5) * c / frequ_1;
datam(:,6) = datam(:,6) * c / frequ_5;

datar(:,5) = datar(:,5) * c / frequ_1;
datar(:,6) = datar(:,6) * c / frequ_5;

obs_master = zeros((max_tow - min_tow), 4 * (np_of_sv - 1));
obs_rover = zeros((max_tow - min_tow), 4 * (np_of_sv - 1));

D = [
    1  -1   -1  1   0   0   0   0   0   0   0   0   0   0   0   0
    1  -1   0   0   -1  1   0   0   0   0   0   0   0   0   0   0
    1  -1   0   0   0   0   -1  1   0   0   0   0   0   0   0   0
    1  -1   0   0   0   0   0   0   -1  1   0   0   0   0   0   0
    1  -1   0   0   0   0   0   0   0   0   -1  1   0   0   0   0
    1  -1   0   0   0   0   0   0   0   0   0   0   -1  1   0   0
    1  -1   0   0   0   0   0   0   0   0   0   0   0   0   -1  1
    ];

A = [   
    1       0           0
    1       0           0
    1       c/frequ_1   0
    1       0           c/frequ_5
    ];

P = diag(1./std.^2);
N = A*A';
D = N(3:4,3:4);
N22 = (D - N(2:3,1)*N(1,2:3)/N(1,1));

i = 0;
for tow = min_tow:15:max_tow
    i = i + 1;
       
    %Get the all satellite for the divent tow
    datam_tow = select_tow(datam, tow);
    datar_tow = select_tow(datar, tow);
    
    %Get the ref satellite
    ref_sat_m = select_prn(datam_tow, prn_number);
    ref_sat_r = select_prn(datar_tow, prn_number);
    
    %Remove from the satellite list the referance sat
    datam_tow = datam_tow(datam_tow(:, 2) ~= prn_number, :);
    datar_tow = datar_tow(datar_tow(:, 2) ~= prn_number, :);
    
    %Create the l matrix with ref sat at the beginning
    l(1,:) = ref_sat_m(3:6);
    l(2,:) = ref_sat_r(3:6);
    l(3:2:(size(datam_tow, 1) * 2 + 1),:) = datam_tow(:,3:6);
    l(4:2:(size(datam_tow, 1) * 2 + 2),:) = datar_tow(:,3:6);
    
    ld = D*l;
    
    for le = ld'
        b = A'*P*le;
        b21 = (b(2)-N(2:3,1)*b(1)/N(1,1));
        x21 = N22^-1*b21;
    end
end

function [data_tow] = select_tow(data, tow)
    data_tow = data(data(:, 1) == tow, :);
end

function [data_tow] = select_prn(data, prn)
    data_tow = data(data(:, 2) == prn, :);
end



