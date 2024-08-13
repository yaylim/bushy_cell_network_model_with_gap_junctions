function [d_out] = bushy_cell_ODE_with_gap_junctions(t,init_cond,g_syn_exct,g_syn_inh,weight_exct,weight_inh,g_gap,num_bushy_cell, V_shift,temp,temp_scaling,step_size_syn,response_type,input_level,network_type,na_type)

global I_ext_val
global input_onset
global input_length
global cluster_out_mat
global num_cell_in_cluster

if strcmp(temp_scaling,'on')
    
    P = (temp-22)/10;
    Q10_time = 3;
    Q10_cond = 2;
    Q_time = Q10_time^P;
    Q_cond = Q10_cond^P;
    temp_cor_time = 1/Q_time;
else
    Q_cond = 1;
    temp_cor_time =1;
end


% defining the internal cell model parameters
% bushy cell model, MC'18
C_m = 26;
if strcmp(na_type,'nav11')
    g_na = 1000;
elseif strcmp(na_type,'nacn')
    g_na = 2300;
end
g_ht = 58;
g_lt = 80;
g_A = 0;
g_h = 30;
g_lk = 2;
V_K = -84;
V_na = 50;
V_h = -43;
V_lk = -65;

% initializing the output vector as a column vector
d_out = zeros(12*num_bushy_cell,1);

for bushy_indx = 1:num_bushy_cell
    
    V(bushy_indx) = init_cond((bushy_indx-1)*12 + 1);
    a(bushy_indx) = init_cond((bushy_indx-1)*12 + 2);
    b(bushy_indx) = init_cond((bushy_indx-1)*12 + 3);
    c(bushy_indx) = init_cond((bushy_indx-1)*12 + 4);
    w(bushy_indx) = init_cond((bushy_indx-1)*12 + 5);
    z(bushy_indx) = init_cond((bushy_indx-1)*12 + 6);
    n(bushy_indx) = init_cond((bushy_indx-1)*12 + 7);
    p(bushy_indx) = init_cond((bushy_indx-1)*12 + 8);
    m(bushy_indx) = init_cond((bushy_indx-1)*12 + 9);
    h(bushy_indx) = init_cond((bushy_indx-1)*12 + 10);
    s(bushy_indx) = init_cond((bushy_indx-1)*12 + 11);
    r(bushy_indx) = init_cond((bushy_indx-1)*12 + 12);
    
    
    
    % fast transient K+ current
    tao_a(bushy_indx) = temp_cor_time*((100*(7*exp((V(bushy_indx)+60)/14) + 29*exp(-(V(bushy_indx)+60)/24))^(-1)) + 0.1);
    tao_b(bushy_indx) = temp_cor_time*((1000*(14*exp((V(bushy_indx)+60)/27) + 29*exp(-(V(bushy_indx)+60)/24))^(-1)) + 1);
    tao_c(bushy_indx) = temp_cor_time*((90*(1 + exp(-(V(bushy_indx)+66)/17))^(-1)) + 10);
    
    a_inf(bushy_indx) = (1 + exp(-(V(bushy_indx)+31)/6))^(-1/4);
    b_inf(bushy_indx) = (1 + exp((V(bushy_indx)+66)/7))^(-1/2);
    c_inf(bushy_indx) = b_inf(bushy_indx);
    
    da(bushy_indx) = (a_inf(bushy_indx) - a(bushy_indx))/tao_a(bushy_indx);
    db(bushy_indx) = (b_inf(bushy_indx) - b(bushy_indx))/tao_b(bushy_indx);
    dc(bushy_indx) = (c_inf(bushy_indx) - c(bushy_indx))/tao_c(bushy_indx);
    
    I_A(bushy_indx)  = g_A * (a(bushy_indx)^4) * b(bushy_indx) * c(bushy_indx)...
        *(V(bushy_indx)-V_K) ;
    
    % low threshold K+ current
    w_inf(bushy_indx) = (1 + exp(-(V(bushy_indx)+48)/6))^(-1/4);
    zeta = 0.5;
    z_inf(bushy_indx) = (1-zeta) * ((1 + exp((V(bushy_indx)+71)/10))^(-1)) + zeta;
    
    tao_w(bushy_indx) = temp_cor_time*((100*(6*exp((V(bushy_indx)+60)/6) + 16*exp(-(V(bushy_indx)+60)/45))^(-1)) + 1.5);
    tao_z(bushy_indx) = temp_cor_time*((1000*(exp((V(bushy_indx)+60)/20) + exp(-(V(bushy_indx)+60)/8))^(-1)) + 50);
    
    dw(bushy_indx) = (w_inf(bushy_indx) - w(bushy_indx))/tao_w(bushy_indx);
    dz(bushy_indx) = (z_inf(bushy_indx) - z(bushy_indx))/tao_z(bushy_indx);
    
    I_lt(bushy_indx) = Q_cond*g_lt * (w(bushy_indx)^4) * z(bushy_indx) * (V(bushy_indx)-V_K);
    
    
    % high threshold K+ current
    sgm = 0.85;
    n_inf(bushy_indx) = (1 + exp(-(V(bushy_indx)+15)/5))^(-1/2);
    p_inf(bushy_indx) = (1 + exp(-(V(bushy_indx)+23)/6))^(-1);
    
    tao_n(bushy_indx) = temp_cor_time*((100*(11*exp((V(bushy_indx)+60+V_shift)/24) + 21*exp(-(V(bushy_indx)+60+V_shift)/23))^(-1)) + 0.7);
    tao_p(bushy_indx) = temp_cor_time*((100*(4*exp((V(bushy_indx)+60+V_shift)/32) + 5*exp(-(V(bushy_indx)+60+V_shift)/22))^(-1)) + 5);
    
    dn(bushy_indx) = (n_inf(bushy_indx) - n(bushy_indx))/tao_n(bushy_indx);
    dp(bushy_indx) = (p_inf(bushy_indx) - p(bushy_indx))/tao_p(bushy_indx);
    
    I_ht(bushy_indx) = Q_cond*g_ht * ((sgm * (n(bushy_indx)^2)) + (1-sgm)*p(bushy_indx) )* (V(bushy_indx)-V_K);
    
    % Fast Na+ Current
    if strcmp(na_type,'nacn')
        % cell 1
        m_inf(bushy_indx) = (1 + exp(-(V(bushy_indx)+38)/7))^(-1);
        h_inf(bushy_indx) = (1 + exp((V(bushy_indx)+65)/6))^(-1);
        
        tao_m(bushy_indx) = temp_cor_time*((10*(5*exp((V(bushy_indx)+60)/18) + 36*exp(-(V(bushy_indx)+60)/25))^(-1)) + 0.04);
        tao_h(bushy_indx) = temp_cor_time*((100*(7*exp((V(bushy_indx)+60)/11) + 10*exp(-(V(bushy_indx)+60)/25))^(-1)) + 0.6);
        
        dm(bushy_indx) = (m_inf(bushy_indx) - m(bushy_indx))/tao_m(bushy_indx);
        dh(bushy_indx) = (h_inf(bushy_indx) - h(bushy_indx))/tao_h(bushy_indx);
        ds(bushy_indx) = 0;
        
        I_na(bushy_indx) = g_na * (m(bushy_indx)^3) * h(bushy_indx) * (V(bushy_indx)-V_na);
        
        
    elseif strcmp(na_type,'na')
        % cell 1
        m_inf(bushy_indx) = (1 + exp(-(V(bushy_indx)+38)/7))^(-1);
        h_inf(bushy_indx) = (1 + exp((V(bushy_indx)+65)/6))^(-1);
        
        tao_m(bushy_indx) = temp_cor_time*((10*(5*exp((V(bushy_indx)+60)/18) + 36*exp(-(V(bushy_indx)+60)/25))^(-1)) + 0.04);
        tao_h(bushy_indx) = temp_cor_time*((100*(7*exp((V(bushy_indx)+60)/11) + 10*exp(-(V(bushy_indx)+60)/25))^(-1)) + 0.6);
        
        dm(bushy_indx) = (m_inf(bushy_indx) - m(bushy_indx))/tao_m(bushy_indx);
        dh(bushy_indx) = (h_inf(bushy_indx) - h(bushy_indx))/tao_h(bushy_indx);
        ds(bushy_indx) = 0;
        
        I_nav = Q_cond*g_na * (m(bushy_indx)^3) * h(bushy_indx) * (V(bushy_indx)-V_na);
        
        
    elseif strcmp(na_type,'nav11')
        % cell 1
        m_inf(bushy_indx) = (1 + exp(-(V(bushy_indx)+27.4+V_shift)*4.7*0.03937))^(-1);
        h_inf(bushy_indx) = (1 + exp((V(bushy_indx)+41.9+V_shift)/6.7))^(-1);
        s_inf(bushy_indx) = (1 + exp((V(bushy_indx)+46+V_shift)/6.6))^(-1);
        
        tao_m(bushy_indx) = temp_cor_time*0.15;
        tao_h(bushy_indx) = temp_cor_time*(23.12*exp(-0.5*((V1+77.58+V_shift)/43.92)^2));
        tao_s(bushy_indx) = temp_cor_time*(1000*140.4*exp(-0.5*((V1+71.3+V_shift)/30.9)^2));
        
        dm(bushy_indx) = (m_inf(bushy_indx) - m(bushy_indx))/tao_m(bushy_indx);
        dh(bushy_indx) = (h_inf(bushy_indx) - h(bushy_indx))/tao_h(bushy_indx);
        ds(bushy_indx) = (s_inf(bushy_indx) - s(bushy_indx))/tao_s(bushy_indx);
        
        I_na(bushy_indx) = g_na * (m(bushy_indx)^3) * h(bushy_indx) * s(bushy_indx) * (V(bushy_indx)-V_na);
        
        
    end
    
    
    % Hyperpolarization activated cation current 
    r_inf(bushy_indx) = (1 + exp((V(bushy_indx)+76)/7))^(-1);
    tao_r(bushy_indx) = temp_cor_time*((100000/(237*exp((V(bushy_indx)+60)/12) + 17*exp(-(V(bushy_indx)+60)/14))) + 25);
    dr(bushy_indx) = (r_inf(bushy_indx) - r(bushy_indx))/tao_r(bushy_indx);
    
    I_h(bushy_indx)  = g_h * r(bushy_indx) * (V(bushy_indx)-V_h) ;
    
    
    % leakage current
    I_lk(bushy_indx) = g_lk * (V(bushy_indx)-V_lk);
    
    if strcmp(response_type,'synaptic')
        
        if strcmp(input_level,'sub')
            g_reg = 0.5;
        elseif strcmp(input_level,'supra')
            g_reg = 3;
        elseif strcmp(input_level,'threshold')
            g_reg = 1;
        end

        V_exct = 0; %0
        V_inh = -75; %potassium equilib

        I_ext(bushy_indx) = 0;

        I_e_exct(bushy_indx) = g_reg * weight_exct * g_syn_exct(bushy_indx,floor(t/step_size_syn)+1) * (V(bushy_indx)-V_exct);
        I_e_inh(bushy_indx) = g_reg * weight_inh * g_syn_inh(bushy_indx,floor(t/step_size_syn)+1) * (V(bushy_indx)-V_inh);
        I_e(bushy_indx) = I_e_exct(bushy_indx) + I_e_inh(bushy_indx);


    elseif strcmp(response_type,'ext')

        I_ext(bushy_indx) = 0;
        if (t> input_onset(bushy_indx)) && (t<input_length + input_onset(bushy_indx))
            I_ext(bushy_indx) = I_ext_val(bushy_indx);
        else
            I_ext(bushy_indx) = 0;
        end
        I_e(bushy_indx) = 0;


    end
end

for bushy_indx = 1:num_bushy_cell
    if strcmp(network_type, 'chain')
        if bushy_indx == 1 % beggining of the chain
            I_gap(bushy_indx) = g_gap * (V(bushy_indx) - V(bushy_indx+1));
        elseif bushy_indx == num_bushy_cell % end of the chain
            I_gap(bushy_indx) = g_gap * (V(bushy_indx) - V(bushy_indx-1));
        else % cells in the middle which got effected by both sides
            I_gap(bushy_indx) = g_gap * (2*V(bushy_indx) - V(bushy_indx-1) - V(bushy_indx+1));
            % over here you have to take the difference between 2 cells twice
            % v3 - v2 and v3 - v1 for example, so it's gonna be 2*V_(bushy_indx)
        end
    elseif strcmp(network_type, 'fully_connected')
        cell_indx = setdiff(1:num_bushy_cell, bushy_indx);
        I_gap_presum = 0;
        I_gap(bushy_indx) = 0;
        for k = 1: num_bushy_cell-1
            I_gap_presum = g_gap * (V(bushy_indx) -  V(cell_indx(k)));
            I_gap(bushy_indx) = I_gap_presum + I_gap(bushy_indx);
        end
    elseif strcmp(network_type, 'cluster_connection')
        [idx_row,~] = find(cluster_out_mat==bushy_indx);
        I_gap(bushy_indx) = 0;

        for cluster_counter = 1:length(idx_row)
            cell_indx = setdiff(cluster_out_mat(idx_row(cluster_counter),:), bushy_indx);
            I_gap_presum = 0;
            for k = 1: num_cell_in_cluster-1
                I_gap_presum = g_gap * (V(bushy_indx) -  V(cell_indx(k)));
                I_gap(bushy_indx) = I_gap_presum + I_gap(bushy_indx);
            end

        end

    end

    %  a chain of bushy cells connected to each other via gap junctions.
    %  cell_1  cell_2   cell_3        cell_(k-1)  cell_k
    %  /\ ------>/\------>/\------> .... /\------> /\
    %  \/<-------\/<------\/<------ .... \/<------ \/


    d_out((bushy_indx-1)*12 + 1) = (1/C_m) * (-I_A(bushy_indx) - I_lt(bushy_indx) - I_ht(bushy_indx)...
        - I_na(bushy_indx) - I_h(bushy_indx) - I_lk(bushy_indx) - I_e(bushy_indx)...
        - I_gap(bushy_indx) + I_ext(bushy_indx) ) ;
    d_out((bushy_indx-1)*12 + 2) = da(bushy_indx);
    d_out((bushy_indx-1)*12 + 3) = db(bushy_indx);
    d_out((bushy_indx-1)*12 + 4) = dc(bushy_indx);
    d_out((bushy_indx-1)*12 + 5) = dw(bushy_indx);
    d_out((bushy_indx-1)*12 + 6) = dz(bushy_indx);
    d_out((bushy_indx-1)*12 + 7) = dn(bushy_indx);
    d_out((bushy_indx-1)*12 + 8) = dp(bushy_indx);
    d_out((bushy_indx-1)*12 + 9) = dm(bushy_indx);
    d_out((bushy_indx-1)*12 + 10) =     dh(bushy_indx);
    d_out((bushy_indx-1)*12 + 11) = ds(bushy_indx);
    d_out((bushy_indx-1)*12 + 12) = dr(bushy_indx);

end
end