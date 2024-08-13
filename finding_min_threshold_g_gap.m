function [mult_exct_supra_out] = finding_min_threshold_g_gap(g_gap_in,num_bushy_cell)
% in this code we are trying to find the minimum threshold input level
% for a cell in the fully connected network to fire an action potential
% when the other cells in the cluster doesn't receive any inputs. The
% expectation is if we increase the number of fully connected cells, the
% level of subthreshold input should increase. We are also trying to figure
% out the level of gap junctions effect on this.

warning("off")

rng(1000)

lower_lim = 0;
upper_lim = 200;
mult_exct_sub = 0;

g_gap = g_gap_in;
mult_exct_supra = randi(100); % initialize this with a random value. This is the value that you are trying to find the minimum of

my_objective = @(mult_exct_supra) finding_threshold_exct(mult_exct_supra,mult_exct_sub,num_bushy_cell,g_gap);

[best_param ,~] = patternsearch(my_objective,1,[],[],[],[],lower_lim,upper_lim); %,[],[],options);
mult_exct_supra_out = best_param;

end
