clear all;
clc;

rx = zeros(1,181);
tx = zeros(1,181);
video = zeros(1,181);
for tau = 0:180
    %given parameters for ground station
    theta = 50;
    n = 0.6;
    phi = 13;
    d = 35786e3;
    %d = distance(r,theta);
    f_up = 14e9;
    L_out_ground = 1;
    P_out_W = 750;              %power output (W)
    P_out = 10*log10(P_out_W);  %power output (dBW)

    %given parameters for satellite TCR
    GT = -5;
    squelch = -115;
    L_out_sat = 1;
    G_sat = 6;  %Satellite telemetry/control gain (dB)
    k = -228.6; %Boltzmann constant (dBW/HzK)

    %given parmeters for video downlink
    L_out_video = 1;
    d_video = 1.5;
    f_down = f_up - 1748e6;

    %call functions to find loss and gain
    G_ground = gain_parabolic(n,phi,f_up);
    G_video = gain_parabolic(n,d_video,f_down);
    L_FS = path_loss(d);
    L_A = -effective_area(f_up);
    L_R = rain_loss_up(theta,tau);
    L_pol = loss_polarisation(tau);

    %calculate link performance
    L_up = L_FS + L_A + L_R + L_pol;
    EIRP = P_out - L_out_ground + G_ground;
    link_performance = EIRP - L_up + GT - k;

    %calculate bit error rate
    EN_required_BPSK = 12; %E_b/N_0 required (dB)
    R = 1e3;
    L_im = 1.5; %implementation loss (dB)
    EN_rx = link_performance - 10*log10(R) - L_im;

    %calculate link performance for downlink
    GT_down = 38;
    Tx_down_P = 10*log10(30);
    EIRP_down = Tx_down_P - L_out_sat + G_sat;
    L_A_down = -effective_area(f_down);
    L_R_down = rain_loss_down(theta,tau);
    L_down = L_FS + L_A_down + L_R_down + L_pol;
    link_performance_down = EIRP_down - L_down + GT_down - k;

    Tx_down_V = 10*log10(300);
    EIRP_video = Tx_down_V - L_out_video + G_video;
    link_performance_video = EIRP_video - L_down + GT_down - k;

    %calculate BER for downlink
    R_down = 16e3;
    R_video = 90e6;
    EN_required_QPSK = 16;
    EN_tx = link_performance_down - 10*log10(R_down) - L_im;
    EN_video = link_performance_video - 10*log10(R_video) - L_im;
    rx(tau+1) = EN_rx;
    tx(tau+1) = EN_tx;
    video(tau+1) = EN_video;
end

h(1) = plot(0:180,rx,'color','#17A589');
hold on;
h(2) = plot(0:180,tx,'color','#2874A6');
h(3) = plot(0:180,video,'color','#6C3483');
xlabel('Polarisation Misalignment (\circ)');
ylabel('Signal/Noise (dB)');
yline(16,'--');
yline(12,'--');
text(140,18,'QPSK Required');
text(140,14,'BPSK Required');
legend(h(1:3),'Command','Telemetry','Video');

