clear

addpath( genpath( [pwd '/BEZ2018_ANF_model/' ]), '-begin' );
addpath( genpath( [pwd '/Data/' ]), '-begin' );

% setting up simulation environment
BS_cell_type = 'SBC';
threshold = -20;
seed_number = 1000; % for reproducibility
rng(seed_number)
CF_BS = 340; % Characteristic frequency of the bushy cell
F0 = 340; % frequency of the stimulus
network_type = 'fully_connected';
num_bushy_cell = 5; % number of cells in a fully connected network
g_gap = 20; % nS, can be arranged to any arbitrary number
weight_inh_BS = 0; % nS, can be arranged to any arbitrary number
weight_exct_BS_threshold = finding_min_threshold_g_gap(g_gap,num_bushy_cell);

% range of frequencies, 1/32 octave steps, from 200Hz to 32kHz
pip_freq_range = zeros(1,236);
pip_sweep_start = 200; % the characteristic freq of the postsynaptic cell
pip_sweep_end = 32000;
pip_sweep_step = 1/32; % in octaves
pip_freq_range(1) = pip_sweep_start;
k = 2;
while pip_freq_range(k-1) < pip_sweep_end
    pip_freq_range(k) = pip_freq_range(k-1)*(2^(pip_sweep_step));
    k = k+1;
end

weights_range = 0:length(pip_freq_range)-1;

% fully connected bushy cell cluster CF range
[~, CF_close_index] = min(abs(pip_freq_range - CF_BS));
bushy_CF_vec  = pip_freq_range((CF_close_index-floor(num_bushy_cell/2)): (CF_close_index+ceil(num_bushy_cell/2)-1));


% input ANF SR types, 1 -> low, 2 -> mid, 3 -> high SR
input_fiber_type_DS = [1,2,3];
input_fiber_type_TV = [1,2];
input_fiber_type_BS = 3;

fiber_type_count_DS = length(input_fiber_type_DS);
fiber_type_count_TV = length(input_fiber_type_TV);
fiber_type_count_BS = length(input_fiber_type_BS);

% network connectivity parameters
if strcmp(BS_cell_type,'SBC')
    AN_BS_input_count = 3;
else
    AN_BS_input_count = 12;

end

AN_BS_input_sd = 0.05;
DS_BS_input_count = 7;
DS_BS_input_sd = 0.208;
TV_BS_input_count = 6;
TV_BS_input_sd = 0.069;

AN_DS_input_count = 12;
AN_DS_input_sd = 0.4;

AN_TV_input_count = 12;
AN_TV_input_sd = 0.1;

AN_DS_freq_chosen_indx = zeros(num_bushy_cell,DS_BS_input_count);
AN_TV_freq_chosen_indx = zeros(num_bushy_cell,TV_BS_input_count);
AN_BS_freq_chosen_indx = zeros(num_bushy_cell,AN_BS_input_count);

% arranging range of inputs for each cells in the bushy cell cluster
for weight_counter = 1: length(bushy_CF_vec)
    CF_post_syn = bushy_CF_vec(weight_counter);
    % creating the range of weight vector
    [~, CF_close_index] = min(abs(pip_freq_range - CF_post_syn));

    weights_DS_BS = lognpdf(weights_range,log(CF_close_index),DS_BS_input_sd);
    [s_DS] = RandStream.create('mlfg6331_64','Seed',seed_number); %for reproductibility
    DS_freq_chosen_indx = datasample(s_DS,1:length(weights_range),DS_BS_input_count,'Replace',false,...
        'Weights',weights_DS_BS); % draw unique samples according to vaules
    % defined in the probability matrix of weights
    AN_DS_freq_chosen_indx(weight_counter,:) = sort(DS_freq_chosen_indx);

    weights_TV_BS = lognpdf(weights_range,log(CF_close_index),TV_BS_input_sd);
    [s_TV] = RandStream.create('mlfg6331_64','Seed',seed_number); %for reproductibility
    TV_freq_chosen_indx = datasample(s_TV,1:length(weights_range),TV_BS_input_count,'Replace',false,...
        'Weights',weights_TV_BS); % draw unique samples according to vaules
    % defined in the probability matrix of weights
    AN_TV_freq_chosen_indx(weight_counter,:) = sort(TV_freq_chosen_indx);

    weights_AN_BS = lognpdf(weights_range,log(CF_close_index),AN_BS_input_sd);
    [s_AN] = RandStream.create('mlfg6331_64','Seed',seed_number); %for reproductibility
    AN_BS_freq_chosen = datasample(s_AN,1:length(weights_range),AN_BS_input_count,'Replace',false,...
        'Weights',weights_AN_BS); % draw unique samples according to vaules
    % defined in the probability matrix of weights
    AN_BS_freq_chosen_indx(weight_counter,:) = sort(AN_BS_freq_chosen);

end

% create a new AN population or use the existing one
if exist(['ANpopulation_' num2str(seed_number) '.mat'],'file')
    load (['ANpopulation_' num2str(seed_number) '.mat']);
    disp('Loading existing population of AN fibers saved in ANpopulation.mat')
else
    [sponts,tabss,trels] = generateANpopulation_seed(length(pip_freq_range), [1 1 1],seed_number);
end

% stimulus parameters
Fs = 100e3;  % sampling rate in Hz (must be 100, 200 or 500 kHz)
rt = 3.9e-3; % rise/fall time in seconds
irpts = rt*Fs;
ondelay = 0; % onset delay
onbin = ondelay*Fs;
sim_time = 20; % in seconds
step_size = 1/Fs;
step_size_syn = step_size/1e-3; % working in msec scale
t = 0:step_size:(sim_time)- step_size; % time vector for the AN model
t_sim = 0:step_size_syn:(sim_time*1e3)-step_size_syn; % time vector for the synaptic model
nrep_stim = 200; % repetition for the stimulus
stimpts = 25*Fs*1e-3; % 25 msec
T_stim = 1/F0;
stimdb = 50; % stimulus dB SPL
[pin] = create_pin(stimdb,onbin,stimpts,irpts,nrep_stim,F0,Fs,t);

%entering AN model parameters
cohc  = 1.0;    % normal ohc function
cihc  = 1.0;    % normal ihc function
species = 1;    % 1 for cat (2 for human with Shera et al. tuning; 3 for human with Glasberg & Moore tuning)
noiseType = 1;  % 1 for variable fGn (0 for fixed fGn)
implnt = 0;     % "0" for approximate or "1" for actual implementation of the power-law functions in the Synapse
nrep_AN = 1;       % number of stimulus repetitions (e.g., 50);


%% DS Cell Layer
% initializing model parameters
temp_scaling = 'on';
temp = 34;
V_shift = 0; % 0 for RM03, 4.3 for XM13
weight_exct_DS = 20;
cell_type_DS = 'type12'; % cell type
na_type_DS = 'na';
init_cond_DS = load(['rest_val_' cell_type_DS '_RM03_22C_nacn']);
init_cond_DS = init_cond_DS.rest_val;
response_type_DS = 'synaptic_exct';
input_level_DS = 'sub';

DS_layer_out = zeros(num_bushy_cell,DS_BS_input_count, length(t_sim));

for DS_cell_counter = 1: num_bushy_cell
    for cell_counter = 1:DS_BS_input_count
        % choosing the range ANFs for each DS cell
        weights_AN_DS = lognpdf(weights_range,log(AN_DS_freq_chosen_indx(DS_cell_counter,cell_counter)),AN_DS_input_sd);
        [s_AN] = RandStream.create('mlfg6331_64','Seed',seed_number); %for reproductibility
        AN_DS_freq_chosen = datasample(s_AN,1:length(pip_freq_range),AN_DS_input_count,'Replace',false,...
            'Weights',weights_AN_DS); % draw unique samples according to values defined in the probability matrix of weights
        AN_DS_freq_chosen = sort(AN_DS_freq_chosen);
        psth_sum = zeros(1,length(t_sim));
        for AN_counter = 1:AN_DS_input_count
            CF = pip_freq_range(AN_DS_freq_chosen(AN_counter));
            AN_param = weights_range(AN_DS_freq_chosen(AN_counter));
            sponts_concat = [sponts.LS(AN_param) sponts.MS(AN_param) sponts.HS(AN_param)];
            tabss_concat = [tabss.LS(AN_param) tabss.MS(AN_param) tabss.HS(AN_param)];
            trels_concat = [trels.LS(AN_param) trels.MS(AN_param) trels.HS(AN_param)];
            for fiber_type_count = 1:fiber_type_count_DS
                spont = sponts_concat(input_fiber_type_DS(fiber_type_count));
                tabs = tabss_concat(input_fiber_type_DS(fiber_type_count));
                trel = trels_concat(input_fiber_type_DS(fiber_type_count));

                vihc = model_IHC_BEZ2018(pin,CF,nrep_AN,1/Fs,sim_time,cohc,cihc,species);
                psth = model_Synapse_BEZ2018(vihc,CF,nrep_AN,1/Fs,noiseType,implnt,...
                    spont,tabs,trel);

                psth_sum = psth_sum + psth; % adds all the inputs from all of the fibers
            end
        end
        g_syn_exct_DS = exp2syn(psth_sum,'exct',cell_type_DS,step_size_syn); % 0.03 is good for type1c

        [~,d_out] = ode45(@(t,init_cond_DS) DS_TV_ODE(t,init_cond_DS,cell_type_DS,g_syn_exct_DS,response_type_DS,na_type_DS,input_level_DS,weight_exct_DS,temp,temp_scaling,V_shift,step_size_syn),t_sim,init_cond_DS);

        [~,l] = findpeaks(d_out(:,1),'MinPeakHeight',threshold);
        DS_layer_spike_train = zeros(1,length(t_sim));
        DS_layer_spike_train(l) = 1;
        DS_layer_out(DS_cell_counter,cell_counter,:) = DS_layer_spike_train;
    end
end


%% TV Cell Layer
% initializing model parameters
temp_scaling = 'off';
V_shift = 4.3; % 0 for RM03, 4.3 for XM13
weight_exct_TV = 20;
cell_type_TV = 'tv'; % cell type
na_type_TV = 'na';
init_cond_TV = load(['rest_val_' cell_type_TV '_XM13_34C_nacn']);
init_cond_TV = init_cond_TV.rest_val;
response_type_TV = 'synaptic_exct';
input_level_TV = 'sub';

TV_layer_out = zeros(num_bushy_cell,TV_BS_input_count, length(t_sim));
for TV_cell_counter = 1:num_bushy_cell
    for cell_counter = 1:TV_BS_input_count
        % choosing the range ANFs for each TV cell
        weights_AN_TV = lognpdf(weights_range,log(AN_TV_freq_chosen_indx(TV_cell_counter,cell_counter)),AN_TV_input_sd);
        [s_AN] = RandStream.create('mlfg6331_64','Seed',seed_number); %for reproductibility
        AN_TV_freq_chosen = datasample(s_AN,1:length(pip_freq_range),AN_TV_input_count,'Replace',false,...
            'Weights',weights_AN_TV); % draw unique samples according to values defined in the probability matrix of weights
        AN_TV_freq_chosen = sort(AN_TV_freq_chosen);
        psth_sum = zeros(1,length(t_sim));
        for AN_counter = 1:AN_TV_input_count
            CF = pip_freq_range(AN_TV_freq_chosen(AN_counter));
            AN_param = weights_range(AN_TV_freq_chosen(AN_counter));
            sponts_concat = [sponts.LS(AN_param) sponts.MS(AN_param) sponts.HS(AN_param)];
            tabss_concat = [tabss.LS(AN_param) tabss.MS(AN_param) tabss.HS(AN_param)];
            trels_concat = [trels.LS(AN_param) trels.MS(AN_param) trels.HS(AN_param)];
            for fiber_type_count = 1:fiber_type_count_TV
                spont = sponts_concat(input_fiber_type_TV(fiber_type_count));
                tabs = tabss_concat(input_fiber_type_TV(fiber_type_count));
                trel = trels_concat(input_fiber_type_TV(fiber_type_count));

                vihc = model_IHC_BEZ2018(pin,CF,nrep_AN,1/Fs,sim_time,cohc,cihc,species);
                psth = model_Synapse_BEZ2018(vihc,CF,nrep_AN,1/Fs,noiseType,implnt,...
                    spont,tabs,trel);

                psth_sum = psth_sum + psth; % adds all the inputs from all of the fibers
            end
        end
        g_syn_exct_TV = exp2syn(psth_sum,'exct',cell_type_TV,step_size_syn); % 0.03 is good for type1c

        [~,d_out] = ode45(@(t,init_cond_TV) DS_TV_ODE(t,init_cond_TV,cell_type_TV,g_syn_exct_TV,response_type_TV,na_type_TV,input_level_TV,weight_exct_TV,temp,temp_scaling,V_shift,step_size_syn),t_sim,init_cond_TV);

        [~,l] = findpeaks(d_out(:,1),'MinPeakHeight',threshold);
        TV_layer_spike_train = zeros(1,length(t_sim));
        TV_layer_spike_train(l) = 1;
        TV_layer_out(TV_cell_counter,cell_counter,:) = TV_layer_spike_train;
    end
end

%% Bushy Cell Layer
% initializing model parameters
temp_scaling = 'off';
V_shift = 4.3; % 0 for RM03, 4.3 for XM13
na_type_BS = 'nacn'; % nacn or nav11
cell_type_BS = 'bs2'; % cell type
rest_val_BS = load(['rest_val_' cell_type_BS '_XM13_34C_nacn']);
init_cond_BS_one_cell = rest_val_BS.rest_val;
for init_indx = 1:num_bushy_cell
    init_cond_BS((init_indx-1)*12 +1 :(init_indx)*12 ) = init_cond_BS_one_cell;
end
response_type_BS = 'synaptic';

if strcmp(BS_cell_type,'SBC')
    input_level_BS = 'supra';
else
    input_level_BS = 'sub';

end

g_syn_inh_BS = zeros(num_bushy_cell,length(t_sim));
g_syn_exct_BS = zeros(num_bushy_cell,length(t_sim));

for BS_cell_counter = 1:num_bushy_cell
    DS_layer_inh_prep = DS_layer_out(BS_cell_counter,:,:);
    TV_layer_inh_prep = TV_layer_out(BS_cell_counter,:,:);

    BS_input_from_DS = sum(DS_layer_inh_prep,2);
    BS_input_from_TV = sum(TV_layer_inh_prep,2);
    g_syn_inh_from_TV  = exp2syn(BS_input_from_TV,'inh','tv',step_size_syn);
    g_syn_inh_from_DS  = exp2syn(BS_input_from_DS,'inh','type12',step_size_syn);
    g_syn_inh_BS(BS_cell_counter,:) = g_syn_inh_from_TV + g_syn_inh_from_DS;

    for cell_counter = 1:AN_BS_input_count
        weights_AN_BS = lognpdf(weights_range,log(AN_BS_freq_chosen_indx(BS_cell_counter, cell_counter)),AN_BS_input_sd);
        [s_AN] = RandStream.create('mlfg6331_64','Seed',seed_number); %for reproductibility
        AN_BS_freq_chosen = datasample(s_AN,1:length(pip_freq_range),AN_BS_input_count,'Replace',false,...
            'Weights',weights_AN_BS); % draw unique samples according to values defined in the probability matrix of weights
        AN_BS_freq_chosen = sort(AN_BS_freq_chosen);
        psth_sum = zeros(1,length(t_sim));
        for AN_counter = 1:AN_BS_input_count
            CF = pip_freq_range(AN_BS_freq_chosen(AN_counter));
            AN_param = weights_range(AN_BS_freq_chosen(AN_counter));
            sponts_concat = [sponts.LS(AN_param) sponts.MS(AN_param) sponts.HS(AN_param)];
            tabss_concat = [tabss.LS(AN_param) tabss.MS(AN_param) tabss.HS(AN_param)];
            trels_concat = [trels.LS(AN_param) trels.MS(AN_param) trels.HS(AN_param)];
            for fiber_type_count = 1:fiber_type_count_BS
                spont = sponts_concat(input_fiber_type_BS(fiber_type_count));
                tabs = tabss_concat(input_fiber_type_BS(fiber_type_count));
                trel = trels_concat(input_fiber_type_BS(fiber_type_count));

                vihc = model_IHC_BEZ2018(pin,CF,nrep_AN,1/Fs,sim_time,cohc,cihc,species);
                psth = model_Synapse_BEZ2018(vihc,CF,nrep_AN,1/Fs,noiseType,implnt,...
                    spont,tabs,trel);

                psth_sum = psth_sum + psth; % adds all the inputs from all of the fibers
            end
        end

    end

    g_syn_exct_BS(BS_cell_counter,:) = exp2syn(psth_sum,'exct',cell_type_BS,step_size_syn);

end

[~,d_out] = ode45(@(t,init_cond_BS) bushy_cell_ODE_with_gap_junctions(t,init_cond_BS,g_syn_exct_BS,g_syn_inh_BS,weight_exct_BS_threshold,weight_inh_BS,g_gap,num_bushy_cell,V_shift,temp,temp_scaling,step_size_syn,response_type_BS,input_level_BS,network_type,na_type_BS),t_sim,init_cond_BS);

% Find the appropriate indexing for the cells and put them in a cell structure again
BS_spike_train = zeros(num_bushy_cell,length(t_sim));
d_out_vec = zeros(num_bushy_cell,length(t_sim));
fire_rate_sec = zeros(1,num_bushy_cell);
SI = zeros(1,num_bushy_cell);

for bushy_counter = 1: num_bushy_cell
    d_out_vec(bushy_counter,:) = d_out(:,(bushy_counter-1)*12 + 1);
    [~,l] = findpeaks(d_out_vec(bushy_counter,:),'MinPeakHeight',threshold,'MinPeakDistance',100);
    BS_spike_train(bushy_counter,l) = 1;

    d_out_vec_out = zeros(nrep_stim,10000);
    fire_rate = zeros(1,nrep_stim);
    for k = 1:nrep_stim
        d_out_vec_out(k,:) = BS_spike_train(bushy_counter,(k-1)*10000+ 1: k*10000);
        fire_rate(k) = sum(d_out_vec_out(k,1000:2500));

    end

    fire_rate_mean = mean(fire_rate);
    fire_rate_sec(bushy_counter) = (1000/15) * fire_rate_mean;

    % calculating the synchronization index scores from period histograms
    psth_one_T = sum(d_out_vec_out,1);
    psth_SI_prep = psth_one_T(1,1000:2500);

    num_psthbins = ceil(T_stim*Fs);
    periodic_psth_sum = zeros(1,num_psthbins);
    tpsth = (0:(length(psth_SI_prep)-1))/Fs;

    for lp = 1:length(psth_SI_prep)
        phbin = round(rem(2*pi*F0*tpsth(lp),2*pi)/(2*pi*F0)*Fs)+1;
        if phbin == num_psthbins+1
            phbin = 1;
        end
        periodic_psth_sum(phbin) = periodic_psth_sum(phbin)+psth_SI_prep(lp);
    end

    SI_sin = periodic_psth_sum * sin(2*pi*(1:length(periodic_psth_sum))...
        /length(periodic_psth_sum))';
    SI_cos= periodic_psth_sum * cos(2*pi*(1:length(periodic_psth_sum))...
        /length(periodic_psth_sum))' ;

    SI_sin = SI_sin/sum(periodic_psth_sum);
    SI_cos = SI_cos/sum(periodic_psth_sum);


    SI(bushy_counter) = sqrt(SI_sin^2 + SI_cos^2);
end

