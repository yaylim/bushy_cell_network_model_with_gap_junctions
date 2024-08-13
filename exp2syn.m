function [syn_out] = exp2syn(syn_in,connection_type,type,step_size_syn)

t_syn = 0:step_size_syn:100-step_size_syn;

if strcmp(connection_type, 'exct')
    if strcmp(type, 'bs2')
        tao_exct1 = 0.05;
        tao_exct2 = 0.4;
    elseif strcmp(type, 'tv')
        tao_exct1 = 0.05;
        tao_exct2 = 0.2;
    elseif strcmp(type, 'type12')
        tao_exct1 = 0.05;
        tao_exct2 = 0.2;
    end

    t_norm = log(tao_exct2/tao_exct1)*(tao_exct1*tao_exct2)/(tao_exct2-tao_exct1);
    norm_const = 1/(exp(-t_norm/tao_exct2)-exp(-t_norm/tao_exct1));
    g_exct = norm_const * (exp(-t_syn/tao_exct2) - exp(-t_syn/tao_exct1));
    syn_out = filter(g_exct,1,syn_in);

elseif strcmp(connection_type, 'inh')
    tao_inh1 = 0.5;
    tao_inh2 = 4.88;

    t_norm = log(tao_inh2/tao_inh1)*(tao_inh1*tao_inh2)/(tao_inh2-tao_inh1);
    norm_const = 1/(exp(-t_norm/tao_inh2)-exp(-t_norm/tao_inh1));
    g_inh = norm_const * (exp(-t_syn/tao_inh2) - exp(-t_syn/tao_inh1));
    syn_out = filter(g_inh,1,syn_in);
end

end