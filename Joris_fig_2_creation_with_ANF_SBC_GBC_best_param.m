%% raster plot comparison
%% SI Plotting
clear,

load threshold_val.mat
g_gap_vec = logspace(0.01,2,10);
inh_vec = logspace(0.01,2,10);
g_gap_strength_vec_SBC = threshold_val*4;

bushy_CF = 340; % in Hz
db_vec = [-inf, -50:10:120];
db_vec_ticks = -60:10:120;

seed_number = 1000;

g_gap_index_SBC = 8; % the index of the gap junction strength level
inh_index_SBC = 7; % the index of the inhibition level
g_gap_index_GBC = 9; % the index of the gap junction strength level
inh_index_GBC = 1; % the index of the inhibition level
db_level = 50; % in dB
db_index = find(db_vec == db_level); % the index of the db level
cell_number = 3;


pip_sweep_start = 200; % the characteristic freq of the postsynaptic cell
pip_sweep_end = 32000;
pip_sweep_step = 1/32; % in octaves
pip_freq_range(1) = pip_sweep_start;
k = 2;
while pip_freq_range(k-1) < pip_sweep_end
    pip_freq_range(k) = pip_freq_range(k-1)*(2^(pip_sweep_step));
    k = k+1;
end
pip_count = length(pip_freq_range);
[~, freq_index] = min(abs(pip_freq_range - bushy_CF));

load (['Data/input_instances_CF' num2str(bushy_CF) '_seed_number_' num2str(seed_number) '.mat'])
load (['Data/gap_junc_network_init_output_' num2str(bushy_CF) 'Hz_seed_number_' num2str(seed_number) '.mat'])
load ('Data/gbc_best_param_grid_search_g_inh_1.0233_g_gap_60.102_weight_exct_mult_77.5598_CF_344Hz.mat')
GBC_stored = BS_stored;
SI_GBC = SI;
fire_rate_sec_GBC = fire_rate_sec;
load ('Data/gap_junct_network_bs_layer_out_340_g_gap_21.7103_inh_36.1225_exct_319.5559_high_SR_seed_1000_3high.mat')


%% calculating the SI_AN and fire_rate_AN
for db_vec_index = 1:19

    d_out_vec_spike = zeros(1,length(t));
    d_out_vec_spike_BS = zeros(1,length(t));
    d_out_vec_spike_GBC = zeros(1,length(t));

    d_out_vec_spike(psthstore{db_vec_index,freq_index,3}) = 1;
    d_out_vec_spike_BS(BS_stored{3,db_vec_index}) = 1;
    d_out_vec_spike_GBC(GBC_stored{3,db_vec_index}) = 1;
    
    BS_vec_spike(db_vec_index,:) = d_out_vec_spike_BS;
    GBC_vec_spike(db_vec_index,:) = d_out_vec_spike_GBC;
    AN_vec_spike(db_vec_index,:) = d_out_vec_spike;
    
    d_out_vec_out = zeros(200,10000);

    for k = 1:200
        d_out_vec_out(k,:) = d_out_vec_spike((k-1)*10000 + 1: k*10000);
        fire_rate(k) = sum(d_out_vec_out(k,1000:2500));
    end

    fire_rate_mean = mean(fire_rate);
    fire_rate_sec_AN(db_vec_index) = (1000/15) * fire_rate_mean;
    psth_one_T = sum(d_out_vec_out,1);
    psth_SI_prep = psth_one_T(1000:2500);

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


    SI_AN_calc(db_vec_index) = sqrt(SI_sin^2 + SI_cos^2);

end


%% Joris Fig_1

db_vec_plot = [12 13 14];

figure,
t_layout = tiledlayout(3,4);
nexttile(1)
title('ANF 340Hz Fire rate and SI')
yyaxis left
plot(db_vec_ticks,fire_rate_sec_AN,'*'), ylim([0 400]),ylabel('Firing rate/sec')
hold on, plot(db_vec_ticks,fire_rate_sec_AN)
yyaxis right
plot(db_vec_ticks,SI_AN_calc,'o'),ylabel('Synchrony'),hold on, plot(db_vec_ticks,SI_AN_calc,'--'),ylim([0 1])
xticks([-60 -40 0 50 100])
xticklabels({'-inf', '-40', '0', '50', '100'})
xlim([-60 120])

for m = 1:length(db_vec_plot)
    nexttile(m+1)
    kl = db_vec_plot(m);
    AN_spike_train = AN_vec_spike(kl,:);

    t_raster = 0:1/Fs: (100e-3)-1/Fs ;
    % raster plotting
    for k = 1:200
        if isempty(find(AN_spike_train((k-1)*10000 +1:k*10000)==1))
            plot([],'bo','markerfacecolor','b','markersize',2)
            hold on
        else
            plot(t_raster(find(AN_spike_train((k-1)*10000 +1:k*10000)==1)),k * ones(1,length(find(AN_spike_train((k-1)*10000 +1:(k)*10000 )==1))),'bo','markerfacecolor','b','markersize',2)
            hold on
        end
    end
    xlim([0 0.035000])
    title(['ANF = ' num2str(bushy_CF) 'Hz, db = ' num2str(db_vec(kl))])

end

%% SBC 
nexttile(5)
title('SBC 340Hz Fire rate and SI')

yyaxis left
plot(db_vec_ticks, fire_rate_sec(3,:),'*'), ylim([0 400]),ylabel('Firing rate/sec')
hold on, plot(db_vec_ticks,fire_rate_sec(3,:))
yyaxis right
plot(db_vec_ticks,SI(3,:),'o'),ylabel('Synchrony'),hold on, plot(db_vec_ticks,SI(3,:),'--'),ylim([0 1])
xticks([-60 -40 0 50 100])
xticklabels({'-inf', '-40', '0', '50', '100'})
xlim([-60 120])

for m = 1:length(db_vec_plot)
    nexttile(5+m)
    kl = db_vec_plot(m);
    for kk = 3
        psth_raster = BS_vec_spike(kl,:);

        for k = 1:nrep_stim
            if isempty(find(psth_raster((k-1)*10000 +1:k*10000)==1))
                plot([],'bo','markerfacecolor','b','markersize',2)
                hold on
            else
                plot(t_raster(find(psth_raster((k-1)*10000 +1:k*10000)==1)),k * ones(1,length(find(psth_raster((k-1)*10000 +1:(k)*10000 )==1))),'bo','markerfacecolor','b','markersize',2)
                hold on
            end
        end

        xlim([0 0.035000])
        title(['SBC CF = ' num2str(bushy_CF) 'Hz, db = ' num2str(db_vec(kl))])
    end

end

%% GBC 

nexttile(9)
title('GBC 340Hz Fire rate and SI')

yyaxis left
plot(db_vec_ticks, fire_rate_sec_GBC(3,:),'*'), ylim([0 400]),ylabel('Firing rate/sec')
hold on, plot(db_vec_ticks,fire_rate_sec_GBC(3,:))
yyaxis right
plot(db_vec_ticks,SI_GBC(3,:),'o'),ylabel('Synchrony'),hold on, plot(db_vec_ticks,SI_GBC(3,:),'--'),ylim([0 1])
xticks([-60 -40 0 50 100])
xticklabels({'-inf', '-40', '0', '50', '100'})
xlim([-60 120])
xlabel('dB SPL')

for m = 1:length(db_vec_plot)
    nexttile(9+m)
    kl = db_vec_plot(m);
    for kk = 3
        psth_raster = GBC_vec_spike(kl,:);

        for k = 1:nrep_stim
            if isempty(find(psth_raster((k-1)*10000 +1:k*10000)==1))
                plot([],'bo','markerfacecolor','b','markersize',2)
                hold on
            else
                plot(t_raster(find(psth_raster((k-1)*10000 +1:k*10000)==1)),k * ones(1,length(find(psth_raster((k-1)*10000 +1:(k)*10000 )==1))),'bo','markerfacecolor','b','markersize',2)
                hold on
            end
        end

        xlim([0 0.035000])
        xlabel('time (sec)')
        title(['GBC CF = ' num2str(bushy_CF) 'Hz, db = ' num2str(db_vec(kl))])
    end

end


title (t_layout,['ANF, SBC and GBC CF = ' num2str(bushy_CF) 'Hz, g\_gap = ' num2str(g_gap_vec(g_gap_index_SBC)), 'nS, g\_inh =  ' num2str(inh_vec(inh_index_SBC)) 'nS'])

