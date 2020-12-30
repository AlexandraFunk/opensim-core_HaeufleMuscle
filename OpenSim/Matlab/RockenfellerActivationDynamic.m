% This script is a Matlab implementation of the Hatze activation dynamic,
% based on the publication of Rockenfeller and Günther in 2018 and shall
% provide a verification of the implementation which was included in
% OpenSim.
% 
% 
% In general this file takes the fiberlength of the contractile element of 
% a muscle which is available in a .sto output file from an OpenSim 
% simulation as input data file, to generate based on the provided data the
% corresponding activation dynamics. 
%

%% USER INPUT
% define start values from M_elbow_flexor_r.dat
lceopt = 0.2127;
time_constant_hatze = 11.3; % [1/s]
nue = 3;
rho_0 = 5.27; % [L/mol*10^4]
gamma_C = 1.37; % [mol/L*10^-4] 
Kuh0 = 0.005; % [] Hatze
gamma0 = 0.01; % initla condition for gamma to meet OpenSim implementation
excitation = 1; % or 1 ???

% get also the data from the ***states.sto file
% read sto file:
relative_path_to_sto_file = "..\Examples\TestOwnRockenfellerMuscle\build\RockenfellerActivationDynamics_states.sto";
full_path_to_sto_file = fullfile(pwd, relative_path_to_sto_file);

%% Read .sto file to matlab:
% open file
stofile = fopen(full_path_to_sto_file,'r');
% skip first 12 lines
for k=1:9
    fgets(stofile);
end
% get the data for time / Rockenfeller / original model 
data = readmatrix(full_path_to_sto_file, 'FileType', 'text');
time = data(:,1);
gamma_rockenfeller_cpp = data(:,14); % TODO check which activation value this is!! seems to be calcium concentration!!
lce_rockenfeller_cpp = data(:,15);
% a_millard_cpp = data(:,16);
% lce_millard_cpp = data(:,17);

% delete unused big data file 
clearvars data;
fclose(stofile);

%% compute activation dynamics based on Hatze formula
a_rockenfeller_matlab = zeros(length(time),1);
gamma_rockenfeller_matlab = zeros(length(time),1);
gamma = gamma0;
for sim_stepper = 1:length(time)
    % calculate new activation with current activation and current lce   
    curent_lce = lce_rockenfeller_cpp(sim_stepper,1);
    rho = (gamma_C * rho_0 * curent_lce / lceopt);
    if (gamma < 0) 
        gamma = 0;
    end
    rhogam = (gamma * rho).^nue;
    a_rockenfeller_matlab(sim_stepper,1) = ((Kuh0 + rhogam)/(1.0 + rhogam));
    gamma_rockenfeller_matlab(sim_stepper,1) = gamma;
    gamma_dot = time_constant_hatze * (excitation - gamma);
    % last time step cant calculate dt
    if sim_stepper < length(time)
        dt = time(sim_stepper+1,1) - time(sim_stepper,1);
        gamma = gamma + dt * gamma_dot;
    end
end



%% output relevant data
figure;
title("Ca2+ concentration calculated with Matlab and C++");
plot(time,gamma_rockenfeller_matlab);
hold on;
plot(time, gamma_rockenfeller_cpp);
legend("Ca2+_{Matlab}", "Ca2+_{C++}");
ylabel("Normalized Ca2+ concentration [-]");
xlabel("Simulation time [s]");

