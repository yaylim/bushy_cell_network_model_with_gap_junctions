function [cost] = finding_threshold_exct(mult_exct_supra,mult_exct_sub, num_bushy_cell,g_gap)

network_type = 'fully_connected';
na_type = 'nacn'; % nacn or nav11
temp_scaling = 'on';
temp = 34;
V_shift = 4.3; % 0 for RM03, 4.3 for XM13
% setting up simulation environment
response_type = 'synaptic'; 
input_level = 'threshold';
% arranging the ODE initial values
type = 'bs2'; % cell type
load(['rest_val_' type '_XM13_34C_nacn']); %v2 is temp scaling on

Fs = 100e3;  % sampling rate in Hz (must be 100, 200 or 500 kHz)
T = 0.05;
t = 0:1/Fs:T-1/Fs; % time vector
mxpts = length(t);
sim_time = T; 
step_size = 1/Fs;
step_size_syn = step_size/1e-3; % working in msec scale
t_sim = 0:step_size_syn:(sim_time*1e3)-step_size_syn;

pin = zeros(1,mxpts); % there is a problem with the rise and fall timer.
pin(1000) = 1;
g_syn_exct_BS = zeros(num_bushy_cell,length(t));
g_syn_inh_BS = zeros(num_bushy_cell,length(t));
g_syn_exct_BS(1,:) = mult_exct_supra * exp2syn(pin,'exct',type,step_size_syn);

for kk = 1:num_bushy_cell-1
    g_syn_exct_BS(1+kk,:) =  mult_exct_sub * exp2syn(pin,'exct',type,step_size_syn);
end

threshold_val = -20;

for init_indx = 1:num_bushy_cell
    init_cond_BS((init_indx-1)*12 +1 :(init_indx)*12 ) = rest_val;
end

weight_exct_BS = 1; 
weight_inh_BS = 0;

[~,d_out] = ode45(@(t,init_cond_BS) bushy_cell_ODE_with_gap_junctions(t,init_cond_BS,g_syn_exct_BS,g_syn_inh_BS,weight_exct_BS,weight_inh_BS,g_gap,num_bushy_cell,V_shift,temp,temp_scaling,step_size_syn,response_type,input_level,network_type,na_type),t_sim,init_cond_BS);
v_out = d_out(:,1); % second cell doesn't receive anything, we are checking the threshold for the first cell 

cost = abs(max(v_out) - (threshold_val)); % trying to find the minimum I_sub value that will give us a spike higher than threshold which is -20mV

end

