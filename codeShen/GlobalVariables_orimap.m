global  num_units; % number of units being simulated
global  dt; % time step
global  t ; % simulation time 
global save_path; % save path of figures
global sigma_ff sigma_ori
global Wee Wei Wie Wii G_Wee G_Wei G_Wie G_Wii
global z
global grid_sz;
global J sigma_J
global tau_E_m tau_I_m
global k_m n_m
global kk

grid_sz = 20;
num_units = grid_sz*grid_sz; % 75 x 75 grid
delta_x = 16/grid_sz;

kk = [0.1 ;0.5];

% ke = 0.1;
% ki = 0.5;

% Jee = 0.1;
% Jei = -0.089;
% Jie = 0.38;
% Jii=-0.096;

J = [0.1  0.089; 0.38 0.096]; 
sigma_J = [8 4;12 4];

k_m = 0.012;
n_m = [2;2.2];

sigma_ori = 45;
sigma_ff = 32;
sigma_rf = delta_x;

dt = 1e-3;
t = 0:dt:1;

tau_E_m = 20e-3;
tau_I_m = 10e-3;

save_path='/Users/shen/Google Drive/private/lab/SSN/Results/';
