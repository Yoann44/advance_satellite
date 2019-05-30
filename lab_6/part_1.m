close all;
clear;
clc;

%Chose the SV used for the double fifference :
%prn_number = 24;
prn_number = 8;

% Load master & rover observations
load('datar.mat');
load('datam.mat');

%Extract tow and number of sv from the data
tow = sort(unique(datam(:,1)));
prns = unique(datam(:,2));
nb_of_sv = size(prns, 1);

%GPS band frequencies [Hz]
frequ_1 = 1575.42e6;
frequ_2 = 1227.6e6;
frequ_5 = 1176.45e6;

%F1 and F2 for E1 respectively E5a frequencies
F1=154;
F2=115;

%Deviation [m] for [C1 C5 L1 L5]
std = [0.5 0.5 0.01 .01];

%Speed of light [m/s]
c = 299792458;

%Convert cycle to meters in L observations for master and rover
datam(:,5) = datam(:,5) * c / frequ_1;
datam(:,6) = datam(:,6) * c / frequ_5;

datar(:,5) = datar(:,5) * c / frequ_1;
datar(:,6) = datar(:,6) * c / frequ_5;

%Create the data structure to save the results
observations = zeros((nb_of_sv - 1), 2,size(tow, 1));
k = zeros((nb_of_sv - 1), 4, size(tow, 1));

%Create the D matrix used for the double difference
D = zeros((nb_of_sv-1), 2*nb_of_sv);
for i = 1:(nb_of_sv-1)
    D(i,1:2) = [1 -1];
    D(i,(i*2+1):(i*2+2)) = [-1 1];
end

%Create the A, P, N, N22 matrices
A = zeros(4, 3);
A(:, 1) = 1;
A(3, 2) = c/frequ_1;
A(4, 3) = c/frequ_5;

P = diag(1./std.^2);
N = A'*P*A;
N22_init = (N(2:3,2:3) - N(2:3,1)*N(1,2:3)/N(1,1));
N22 = N22_init;

%Temporary variable used to remember last iteration b
b21_last = zeros((nb_of_sv - 1), 2);

for i = 1:size(tow, 1)
    %Get the all satellite for the divent tow
    datam_tow = select_tow(datam, tow(i));
    datar_tow = select_tow(datar, tow(i));
    
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
    
    j = 0;
    for le = ld'
        j = j + 1;
        b = A'*P*le;
        b21 = b(2:3)-N(2:3,1)*b(1)/N(1,1);
        
        if i > 1
            b21 = b21 + b21_last(j,:)';
        end
        b21_last(j,:) = b21';
        
        %Evaluate and save the observation
        x21 = N22^-1*b21;
        observations(j, :, i) = x21;
    end
    
    k(:,1,i) = observations(:, 1, i) - observations(:, 2, i);
    k(:,2,i) = F2 .* observations(:, 1, i) - F1 .* observations(:, 2, i);
    k(:,3,i) = round(k(:,1,i));
    k(:,4,i) = round(k(:,2,i));
    
    %Prepare N22 for the next iteration
    N22 = N22 + N22_init;
end

%Show some results
sprintf("Epoch : %d", tow(size(tow, 1)));
tow_index = size(tow, 1);
for i = 1:(nb_of_sv-1)
    fprintf("DD(%d - %d)\t: N1=\t%2.1f(\t%d) N2=\t%2.2f(\t%d) Nw=\t%2.1f\n", ...
        prn_number, prns(i), ...
        observations(i,1,tow_index), round(observations(i,1,tow_index)), ...
        observations(i,2,tow_index), round(observations(i,2,tow_index)), ...
        k(i,1,tow_index));
end

function [data_tow] = select_tow(data, tow)
    data_tow = data(data(:, 1) == tow, :);
end

function [data_tow] = select_prn(data, prn)
    data_tow = data(data(:, 2) == prn, :);
end



