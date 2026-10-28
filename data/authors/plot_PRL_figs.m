% Plot Figure 3
clear
clf
load prl_fig3_data
% do a least square fit (constant slope line)
pp = polyfit(ratio_SS_CS,ratio_det_latency,1);
xi=1:0.01:4;
yi=polyval(pp,xi);
pp_bpsk = polyfit(ratio_SS_CS,ratio_det_latency_bpsk,1);
yi_bpsk=polyval(pp_bpsk,xi);

% scatter plot of ratio of latencies vs ratio of CRE (homodyne) for squeezed / classical
scatter(ratio_SS_CS,ratio_det_latency,"b","filled")
hold on
scatter(ratio_SS_CS,ratio_det_latency_bpsk,"r","filled")
plot(xi,yi,'color',"b")
plot(xi,yi_bpsk,'color',"r")

grid("on")
title(sprintf('N=%d, N_a in [%.2f,%.1f], loss change at t=%d (%d samples), ARL_0=%.1e',N,sinh(r_min)^2,N_a,cp,n,ARL0),'Fontsize',9)
xlabel('Ratio of relative entropies (increases with N_a)','Fontsize',10)
ylabel('Ratio of change-detection latencies','Fontsize',10)
legend("Unmodulated",  "BPSK modulation", "Linear least squares fit (Unmodulated)","Linear least squares fit (BPSK)", "location", "northwest")
print("PRL_Guha_Cherian-John_Basu_Main_Fig3.png")

% Plot Figure 7
clear
clf
load prl_fig7_data
plot(Na_range, mean_clat, 'linewidth', 2, 'marker', 'o', ";Classical;");
set(gca, "linewidth", 2, "fontsize", 12)
hold on;
styles = ['+','*','x','s'];
for m=1:length(blocksizes)
  plot(Na_range, mean_sqentlat{m}, 'linewidth', 2, 'marker', styles(m), sprintf(";Entangled (blocksize=%d);",blocksizes(m)));
endfor
legend('location','east');
title(sprintf('N=%d, loss change at t=%d (%d samples), ARL_0=%.1e',N,cp,n,ARL0),'Fontsize',9);
xlabel('Number of squeezed (or entangled) photons per pulse (N_a)','Fontsize',10);
ylabel('Mean detection latency (time in number of pulses)','Fontsize',10);
print("PRL_Guha_Cherian-John_Basu_Main_Fig7.png")


% Plot Figure 4(b) in Supplemental Material
% plot the latency vs dB squeezing
clear
clf
load prl_supplem_fig4b_data
dv_lat = mean(lat);
dBsqueeze = 10 * log10(exp(2.*r_range ));
scatter(dBsqueeze,c_lat_pos,"b","*", "linewidth", 1)
hold on
scatter(dBsqueeze,sq_lat_pos,"r","o", "linewidth", 1)
scatter(dBsqueeze,c_lat_bpsk_pos,"c","s", "linewidth", 1)
scatter(dBsqueeze,sq_lat_bpsk_pos,"g","x", "linewidth", 1)
xlim([0 8])
line ("xdata",[0.8,7], "ydata",[dv_lat,dv_lat], "linestyle", "--", "linewidth", 1.25)
line ("xdata",[1.83,1.83], "ydata",[0,55], "color", "magenta", "linestyle", "-.", "linewidth", 0.8)
text(0.2,51.5,"DV, N_a = 1","Fontsize",11)
text(1.9,54,"CV: needs steady stream of (deterministic) 1.83dB squeezed pulses", "Fontsize", 9)
text(1.9,51.5,"DV: needs steady stream of (deterministic) source of 1 photons", "Fontsize", 9)
grid("on")
set(gca,'GridLineStyle','--')
%title(sprintf('N=%d, N_a in [%.2f,%.1f], loss change at t=%d (%d samples), ARL_0=%.1e slots',N,sinh(r_min)^2,N_a,cp,n,ARL0),'Fontsize',12,'FontWeight','Normal')
xlabel('(Squeezing per pulse, in dB) 10log_{10}(e^{2r})', 'Fontsize', 11)
ylabel('(Change detection latency) \tau', 'Fontsize', 11)
legend("CV: Weak coherent pulses (Unmodulated)", "CV: Squeezing-augmented pulses (Unmodulated)", "CV: Weak coherent pulses (BPSK)", "CV: Squeezing-augmented pulses (BPSK)", "DV: single photon pulses (Unmodulated)", "location", "southwest", "Fontsize", 8, "box", "off")
print("PRL_Guha_Cherian-John_Basu_Supplement_Fig4b.png")
