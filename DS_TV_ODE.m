function [d_out] = DS_TV_ODE(t,init_cond,type,g_syn_exct,response_type,na_type,input_level,weight_exct,temp, temp_scaling, V_shift,step_size_syn)

global I_ext_val
global input_onset
global input_length

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

if strcmp(type,'tv')
    % tuberculoventral cell model, MC'18
    C_m = 35;
    g_na = 5800;
    g_ht = 400;
    g_lt = 0;
    % in the data table provided by Manis & Campagnola 2018 model,
    % there is no implication for klt channel, but since it has been said that
    % the model is based on tstellate cell model, we assumed g_lt would be zero
    g_A = 65;
    g_h = 2.5;
    g_lk = 4.5;
    V_K = -81.5;
    V_na = 50;
    V_h = -43;
    V_lk = -72;

elseif strcmp(type,'type12')
    % d stellate cell model, RM'03
    C_m = 12;
    g_na = 1000;
    g_ht = 150;
    g_lt = 20;
    g_A = 0;
    g_h = 2;
    g_lk = 2;
    V_K = -70;
    V_na = 55;
    V_h = -43;
    V_lk = -65;

end

d_out = zeros(12,1);    % a column vector

V = init_cond(1);
a = init_cond(2);
b = init_cond(3);
c = init_cond(4);
w = init_cond(5);
z = init_cond(6);
n = init_cond(7);
p = init_cond(8);
m = init_cond(9);
h = init_cond(10);
s = init_cond(11);
r = init_cond(12);

% fast transient K+ current
tao_a = temp_cor_time*((100*(7*exp((V+60)/14) + 29*exp(-(V+60)/24))^(-1)) + 0.1);
tao_b = temp_cor_time*((1000*(14*exp((V+60)/27) + 29*exp(-(V+60)/24))^(-1)) + 1);
tao_c = temp_cor_time*((90*(1 + exp(-(V+66)/17))^(-1)) + 10);

a_inf = (1 + exp(-(V+31)/6))^(-1/4);
b_inf = (1 + exp((V+66)/7))^(-1/2);
c_inf = b_inf;

da = (a_inf - a)/tao_a;
db = (b_inf - b)/tao_b;
dc = (c_inf - c)/tao_c;

I_A  = Q_cond * g_A * (a^4) * b * c *(V-V_K) ;

% low threshold K+ current
w_inf = (1 + exp(-(V+48)/6))^(-1/4);
zeta = 0.5;
z_inf = (1-zeta) * ((1 + exp((V+71)/10))^(-1)) + zeta;

tao_w = temp_cor_time*((100*(6*exp((V+60)/6) + 16*exp(-(V+60)/45))^(-1)) + 1.5);
tao_z = temp_cor_time*((1000*(exp((V+60)/20) + exp(-(V+60)/8))^(-1)) + 50);

dw = (w_inf - w)/tao_w;
dz = (z_inf - z)/tao_z;

I_lt = Q_cond * g_lt * (w^4) * z * (V-V_K);


% high threshold K+ current
sgm = 0.85;
n_inf = (1 + exp(-(V+15)/5))^(-1/2);
p_inf = (1 + exp(-(V+23)/6))^(-1);

tao_n = temp_cor_time*((100*(11*exp((V+60+V_shift)/24) + 21*exp(-(V+60+V_shift)/23))^(-1)) + 0.7);
tao_p = temp_cor_time*((100*(4*exp((V+60+V_shift)/32) + 5*exp(-(V+60+V_shift)/22))^(-1)) + 5);

dn = (n_inf - n)/tao_n;
dp = (p_inf - p)/tao_p;

I_ht = Q_cond * g_ht * ((sgm * (n^2)) + (1-sgm)*p )* (V-V_K);

% Fast Na+ Current
if strcmp(na_type,'nacn')
    m_inf = (1 + exp(-(V+38)/7))^(-1);
    h_inf = (1 + exp((V+65)/6))^(-1);

    tao_m = temp_cor_time*((10*(5*exp((V+60)/18) + 36*exp(-(V+60)/25))^(-1)) + 0.04);
    tao_h = temp_cor_time*((100*(7*exp((V+60)/11) + 10*exp(-(V+60)/25))^(-1)) + 0.6);

    dm = (m_inf - m)/tao_m;
    dh = (h_inf - h)/tao_h;
    ds = 0;

    I_na = Q_cond * g_na * (m^3) * h * (V-V_na);


elseif strcmp(na_type,'na')
    m_inf = (1 + exp(-(V+38)/7))^(-1);
    h_inf = (1 + exp((V+65)/6))^(-1);

    tao_m = temp_cor_time*((10*(5*exp((V+60)/18) + 36*exp(-(V+60)/25))^(-1)) + 0.04);
    tao_h = temp_cor_time*((100*(7*exp((V+60)/11) + 10*exp(-(V+60)/25))^(-1)) + 0.6);

    dm = (m_inf - m)/tao_m;
    dh = (h_inf - h)/tao_h;
    ds = 0;

    I_na = Q_cond * g_na * (m^3) * h * (V-V_na);

elseif strcmp(na_type,'nav11')
    m_inf = (1 + exp(-(V+27.4+V_shift)*4.7*0.03937))^(-1);
    h_inf = (1 + exp((V+41.9+V_shift)/6.7))^(-1);
    s_inf = (1 + exp((V+46+V_shift)/6.6))^(-1);

    tao_m = temp_cor_time*0.15;
    tao_h = temp_cor_time*(23.12*exp(-0.5*((V+77.58+V_shift)/43.92)^2));
    tao_s = temp_cor_time*(1000*140.4*exp(-0.5*((V+71.3+V_shift)/30.9)^2));

    dm = (m_inf - m)/tao_m;
    dh = (h_inf - h)/tao_h;
    ds = (s_inf - s)/tao_s;

    I_na = Q_cond *  g_na * (m^3) * h * s * (V-V_na);

end


% Hyperpolarization activated cation current
r_inf = (1 + exp((V+76)/7))^(-1);
tao_r = temp_cor_time * ((100000/(237*exp((V+60)/12) + 17*exp(-(V+60)/14))) + 25);
dr = (r_inf - r)/tao_r;

I_h  = Q_cond * g_h * r * (V-V_h) ;

% leakage current
I_lk = g_lk * (V-V_lk);

if strcmp(response_type,'synaptic_exct')

    if strcmp(input_level,'sub')
        g_reg = 0.5;
    else
        g_reg = 3;
    end

    V_exct = 0; %0
    V_inh = -75; %potassium equilib
    I_ext = 0;
    I_e_exct = g_reg * weight_exct * g_syn_exct(floor(t/step_size_syn)+1) * (V-V_exct);
    I_e = I_e_exct;

elseif strcmp(response_type,'ext')

    I_ext = 0;
    if (t> input_onset) && (t<input_length + input_onset)
        I_ext = I_ext_val;
    else
        I_ext = 0;
    end
    I_e = 0;

end

d_out(1) = (1/C_m) * (-I_A - I_lt - I_ht - I_na - I_h - I_lk - I_e + I_ext);
d_out(2) = da;
d_out(3) = db;
d_out(4) = dc;
d_out(5) = dw;
d_out(6) = dz;
d_out(7) = dn;
d_out(8) = dp;
d_out(9) = dm;
d_out(10) = dh;
d_out(11) = ds;
d_out(12) = dr;

end