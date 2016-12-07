%%
clear all
close all
clc


clear all
close all
clc
set(0,'defaulttextinterpreter','latex')
set(0,'DefaultTextFontname', 'CMU Serif')
set(0,'DefaultAxesFontName', 'CMU Serif')


Pt = 0; %or 15, 10 (neg values, dBm)
Gtx = 6; %dBi
Grx = 24; %dBi
Gt = 6; %dBi
r_th = [20:5:55];
f = 5.8e9; %Hz
lam = 3e8/f;


P_on_tag = Pt + Gtx + Gt + 20*log10((lam)./(4*pi.*r_th)); %dBm

P_r = Pt + Gtx + 2*Gt + Grx + 40*log10((lam)./(4*pi.*r_th)); %No tag amplification, dBm





% fileID = fopen('testio.bin','rb');
% A = fread(fileID,'float');
quadrant = 0;
count = 0;
average = [];
output = [];
all_outputs = [];

output_25m_60mv_250khz = read_complex_binary('25m_60mv_250khz.bin',16072);
output_25m_62mv_250khz = read_complex_binary('25m_62mv_250khz.bin',16072);
output_25m_64mv_250khz =   read_complex_binary('25m_64mv_250khz.bin',16072);
output_25m_66mv_250khz = read_complex_binary('25m_66mv_250khz.bin',16072);
output_25m_68mv_250khz = read_complex_binary('25m_68mv_250khz.bin',16072);
output_25m_60mv_1Mhz = read_complex_binary('25m_60mv_1Mhz.bin',16072);

output_30m_60mv_250khz = read_complex_binary('30m_60mv_250khz.bin',16072);
output_30m_62mv_250khz = read_complex_binary('30m_62mv_250khz.bin',16072);
output_30m_64mv_250khz =   read_complex_binary('30m_64mv_250khz.bin',16072);
output_30m_66mv_250khz = read_complex_binary('30m_66mv_250khz.bin',16072);
output_30m_68mv_250khz = read_complex_binary('30m_68mv_250khz.bin',16072);
output_30m_60mv_1Mhz = read_complex_binary('30m_60mv_1Mhz.bin',16072);

output_35m_60mv_250khz = read_complex_binary('35m_60mv_250khz.bin',16072);
output_35m_62mv_250khz = read_complex_binary('35m_62mv_250khz.bin',16072);
output_35m_64mv_250khz =   read_complex_binary('35m_64mv_250khz.bin',16072);
output_35m_66mv_250khz = read_complex_binary('35m_66mv_250khz.bin',16072);
output_35m_68mv_250khz = read_complex_binary('35m_68mv_250khz.bin',16072);
output_35m_60mv_1Mhz = read_complex_binary('35m_60mv_1Mhz.bin',16072);

output_40m_60mv_250khz = read_complex_binary('40m_60mv_250khz.bin',16072);
output_40m_62mv_250khz = read_complex_binary('40m_62mv_250khz.bin',16072);
output_40m_64mv_250khz =   read_complex_binary('40m_64mv_250khz.bin',16072);
output_40m_66mv_250khz = read_complex_binary('40m_66mv_250khz.bin',16072);
output_40m_68mv_250khz = read_complex_binary('40m_68mv_250khz.bin',16072);
output_40m_60mv_1Mhz = read_complex_binary('40m_60mv_1Mhz.bin',16072);

output_45m_60mv_250khz = read_complex_binary('45m_60mv_250khz.bin',16072);
output_45m_62mv_250khz = read_complex_binary('45m_62mv_250khz.bin',16072);
output_45m_64mv_250khz =   read_complex_binary('45m_64mv_250khz.bin',16072);
output_45m_66mv_250khz = read_complex_binary('45m_66mv_250khz.bin',16072);
output_45m_68mv_250khz = read_complex_binary('45m_68mv_250khz.bin',16072);
output_45m_60mv_1Mhz = read_complex_binary('45m_60mv_1Mhz.bin',16072);

output_50m_60mv_250khz = read_complex_binary('50m_60mv_250khz.bin',16072);
output_50m_62mv_250khz = read_complex_binary('50m_62mv_250khz.bin',16072);
output_50m_64mv_250khz =   read_complex_binary('50m_64mv_250khz.bin',16072);
output_50m_66mv_250khz = read_complex_binary('50m_66mv_250khz.bin',16072);
output_50m_68mv_250khz = read_complex_binary('50m_68mv_250khz.bin',16072);
output_50m_60mv_1Mhz = read_complex_binary('50m_60mv_1Mhz.bin',16072);


[max_value_output_25m_60mv_250khz,max_index_output_25m_60mv_250khz] = max(abs(output_25m_60mv_250khz));
[max_value_output_25m_62mv_250khz,max_index_output_25m_62mv_250khz] = max(abs(output_25m_62mv_250khz));
[max_value_output_25m_64mv_250khz,  max_index_output_25m_64mv_250khz] = max(abs(output_25m_64mv_250khz));
[max_value_output_25m_66mv_250khz,max_index_output_25m_66mv_250khz] = max(abs(output_25m_66mv_250khz));
[max_value_output_25m_68mv_250khz,  max_index_output_25m_68mv_250khz] = max(abs(output_25m_68mv_250khz));
[max_value_output_25m_60mv_1Mhz,max_index_output_25m_60mv_1Mhz] = max(abs(output_25m_60mv_1Mhz));

[max_value_output_30m_60mv_250khz,max_index_output_30m_60mv_250khz] = max(abs(output_30m_60mv_250khz));
[max_value_output_30m_62mv_250khz,max_index_output_30m_62mv_250khz] = max(abs(output_30m_62mv_250khz));
[max_value_output_30m_64mv_250khz,  max_index_output_30m_64mv_250khz] = max(abs(output_30m_64mv_250khz));
[max_value_output_30m_66mv_250khz,max_index_output_30m_66mv_250khz] = max(abs(output_30m_66mv_250khz));
[max_value_output_30m_68mv_250khz,  max_index_output_30m_68mv_250khz] = max(abs(output_30m_68mv_250khz));
[max_value_output_30m_60mv_1Mhz,max_index_output_30m_60mv_1Mhz] = max(abs(output_30m_60mv_1Mhz));

[max_value_output_35m_60mv_250khz,max_index_output_35m_60mv_250khz] = max(abs(output_35m_60mv_250khz));
[max_value_output_35m_62mv_250khz,max_index_output_35m_62mv_250khz] = max(abs(output_35m_62mv_250khz));
[max_value_output_35m_64mv_250khz,  max_index_output_35m_64mv_250khz] = max(abs(output_35m_64mv_250khz));
[max_value_output_35m_66mv_250khz,max_index_output_35m_66mv_250khz] = max(abs(output_35m_66mv_250khz));
[max_value_output_35m_68mv_250khz,  max_index_output_35m_68mv_250khz] = max(abs(output_35m_68mv_250khz));
[max_value_output_35m_60mv_1Mhz,max_index_output_35m_60mv_1Mhz] = max(abs(output_35m_60mv_1Mhz));

[max_value_output_40m_60mv_250khz,max_index_output_40m_60mv_250khz] = max(abs(output_40m_60mv_250khz));
[max_value_output_40m_62mv_250khz,max_index_output_40m_62mv_250khz] = max(abs(output_40m_62mv_250khz));
[max_value_output_40m_64mv_250khz,  max_index_output_40m_64mv_250khz] = max(abs(output_40m_64mv_250khz));
[max_value_output_40m_66mv_250khz,max_index_output_40m_66mv_250khz] = max(abs(output_40m_66mv_250khz));
[max_value_output_40m_68mv_250khz,  max_index_output_40m_68mv_250khz] = max(abs(output_40m_68mv_250khz));
[max_value_output_40m_60mv_1Mhz,max_index_output_40m_60mv_1Mhz] = max(abs(output_40m_60mv_1Mhz));

[max_value_output_45m_60mv_250khz,max_index_output_45m_60mv_250khz] = max(abs(output_45m_60mv_250khz));
[max_value_output_45m_62mv_250khz,max_index_output_45m_62mv_250khz] = max(abs(output_45m_62mv_250khz));
[max_value_output_45m_64mv_250khz,  max_index_output_45m_64mv_250khz] = max(abs(output_45m_64mv_250khz));
[max_value_output_45m_66mv_250khz,max_index_output_45m_66mv_250khz] = max(abs(output_45m_66mv_250khz));
[max_value_output_45m_68mv_250khz,  max_index_output_45m_68mv_250khz] = max(abs(output_45m_68mv_250khz));
[max_value_output_45m_60mv_1Mhz,max_index_output_45m_60mv_1Mhz] = max(abs(output_45m_60mv_1Mhz));

[max_value_output_50m_60mv_250khz,max_index_output_50m_60mv_250khz] = max(abs(output_50m_60mv_250khz));
[max_value_output_50m_62mv_250khz,max_index_output_50m_62mv_250khz] = max(abs(output_50m_62mv_250khz));
[max_value_output_50m_64mv_250khz,  max_index_output_50m_64mv_250khz] = max(abs(output_50m_64mv_250khz));
[max_value_output_50m_66mv_250khz,max_index_output_50m_66mv_250khz] = max(abs(output_50m_66mv_250khz));
[max_value_output_50m_68mv_250khz,  max_index_output_50m_68mv_250khz] = max(abs(output_50m_68mv_250khz));
[max_value_output_50m_60mv_1Mhz,max_index_output_50m_60mv_1Mhz] = max(abs(output_50m_60mv_1Mhz));



max_complex_values = zeros(1,36);
max_complex_values(1) = output_25m_60mv_250khz(max_index_output_25m_60mv_250khz);
max_complex_values(2) = output_25m_62mv_250khz(max_index_output_25m_62mv_250khz);
max_complex_values(3) = output_25m_64mv_250khz(max_index_output_25m_64mv_250khz);
max_complex_values(4) = output_25m_66mv_250khz(max_index_output_25m_66mv_250khz);
max_complex_values(5) = output_25m_68mv_250khz(max_index_output_25m_68mv_250khz);
max_complex_values(6) = output_25m_60mv_1Mhz(max_index_output_25m_60mv_1Mhz);

max_complex_values(7) = output_30m_60mv_250khz(max_index_output_30m_60mv_250khz);
max_complex_values(8) = output_30m_62mv_250khz(max_index_output_30m_62mv_250khz);
max_complex_values(9) = output_30m_64mv_250khz(max_index_output_30m_64mv_250khz);
max_complex_values(10) = output_30m_66mv_250khz(max_index_output_30m_66mv_250khz);
max_complex_values(11) = output_30m_68mv_250khz(max_index_output_30m_68mv_250khz);
max_complex_values(12) = output_30m_60mv_1Mhz(max_index_output_30m_60mv_1Mhz);

max_complex_values(13) = output_35m_60mv_250khz(max_index_output_35m_60mv_250khz);
max_complex_values(14) = output_35m_62mv_250khz(max_index_output_35m_62mv_250khz);
max_complex_values(15) = output_35m_64mv_250khz(max_index_output_35m_64mv_250khz);
max_complex_values(16) = output_35m_66mv_250khz(max_index_output_35m_66mv_250khz);
max_complex_values(17) = output_35m_68mv_250khz(max_index_output_35m_68mv_250khz);
max_complex_values(18) = output_35m_60mv_1Mhz(max_index_output_35m_60mv_1Mhz);

max_complex_values(19) = output_40m_60mv_250khz(max_index_output_40m_60mv_250khz);
max_complex_values(20) = output_40m_62mv_250khz(max_index_output_40m_62mv_250khz);
max_complex_values(21) = output_40m_64mv_250khz(max_index_output_40m_64mv_250khz);
max_complex_values(22) = output_40m_66mv_250khz(max_index_output_40m_66mv_250khz);
max_complex_values(23) = output_40m_68mv_250khz(max_index_output_40m_68mv_250khz);
max_complex_values(24) = output_40m_60mv_1Mhz(max_index_output_40m_60mv_1Mhz);

max_complex_values(25) = output_45m_60mv_250khz(max_index_output_45m_60mv_250khz);
max_complex_values(26) = output_45m_62mv_250khz(max_index_output_45m_62mv_250khz);
max_complex_values(27) = output_45m_64mv_250khz(max_index_output_45m_64mv_250khz);
max_complex_values(28) = output_45m_66mv_250khz(max_index_output_45m_66mv_250khz);
max_complex_values(29) = output_45m_68mv_250khz(max_index_output_45m_68mv_250khz);
max_complex_values(30) = output_45m_60mv_1Mhz(max_index_output_45m_60mv_1Mhz);

max_complex_values(31) = output_50m_60mv_250khz(max_index_output_50m_60mv_250khz);
max_complex_values(32) = output_50m_62mv_250khz(max_index_output_50m_62mv_250khz);
max_complex_values(33) = output_50m_64mv_250khz(max_index_output_50m_64mv_250khz);
max_complex_values(34) = output_50m_66mv_250khz(max_index_output_50m_66mv_250khz);
max_complex_values(35) = output_50m_68mv_250khz(max_index_output_50m_68mv_250khz);
max_complex_values(36) = output_50m_60mv_1Mhz(max_index_output_50m_60mv_1Mhz);



all_outputs(1,:) = output_25m_60mv_250khz;
all_outputs(2,:) = output_25m_62mv_250khz;
all_outputs(3,:) = output_25m_64mv_250khz;
all_outputs(4,:) = output_25m_66mv_250khz;
all_outputs(5,:) = output_25m_68mv_250khz;
all_outputs(6,:) = output_25m_60mv_1Mhz;

all_outputs(7,:) = output_30m_60mv_250khz;
all_outputs(8,:) = output_30m_62mv_250khz;
all_outputs(9,:) = output_30m_64mv_250khz;
all_outputs(10,:) = output_30m_66mv_250khz;
all_outputs(11,:) = output_30m_68mv_250khz;
all_outputs(12,:) = output_30m_60mv_1Mhz;

all_outputs(13,:) = output_35m_60mv_250khz;
all_outputs(14,:) = output_35m_62mv_250khz;
all_outputs(15,:) = output_35m_64mv_250khz;
all_outputs(16,:) = output_35m_66mv_250khz;
all_outputs(17,:) = output_35m_68mv_250khz;
all_outputs(18,:) = output_35m_60mv_1Mhz;

all_outputs(19,:) = output_40m_60mv_250khz;
all_outputs(20,:) = output_40m_62mv_250khz;
all_outputs(21,:) = output_40m_64mv_250khz;
all_outputs(22,:) = output_40m_66mv_250khz;
all_outputs(23,:) = output_40m_68mv_250khz;
all_outputs(24,:) = output_40m_60mv_1Mhz;

all_outputs(25,:) = output_45m_60mv_250khz;
all_outputs(26,:) = output_45m_62mv_250khz;
all_outputs(27,:) = output_45m_64mv_250khz;
all_outputs(28,:) = output_45m_66mv_250khz;
all_outputs(29,:) = output_45m_68mv_250khz;
all_outputs(30,:) = output_45m_60mv_1Mhz;

all_outputs(31,:) = output_50m_60mv_250khz;
all_outputs(32,:) = output_50m_62mv_250khz;
all_outputs(33,:) = output_50m_64mv_250khz;
all_outputs(34,:) = output_50m_66mv_250khz;
all_outputs(35,:) = output_50m_68mv_250khz;
all_outputs(36,:) = output_50m_60mv_1Mhz;


%%
Vpk = [];
Pdbm = [];
SNR = [];

for p = 1:36
        average = [];
        count = 0;
        output = all_outputs(p,:);
        max_complex = max_complex_values(p);
        phase_of_max = angle(max_complex);
        if (phase_of_max > 0 && phase_of_max <= pi/2)
            quadrant = 1;
            for b = 1:length(output)
                if (real(output(b)) > 0 && imag(output(b)) > 0)
                    count = count+1;
                    average(b) = output(b);
                end
            end
            
            
        elseif (phase_of_max > pi/2 && phase_of_max <= pi)
            quadrant = 2;
            for b = 1:length(output)
                if (real(output(b)) < 0 && imag(output(b)) > 0)
                    count = count+1;
                    average(b) = output(b);
                 end
            end

        elseif (phase_of_max > -pi/2 && phase_of_max <= 0)
            quadrant = 4;
            for b = 1:length(output)
                if (real(output(b)) > 0 && imag(output(b)) < 0)
                    count = count+1;
                    average(b) = output(b);
                 end
            end
        else
            quadrant = 3;
                for b = 1:length(output)
                if (real(output(b)) < 0 && imag(output(b)) < 0)
                    count = count+1;
                    average(b) = output(b);
                 end
            end
        end
        average_single(p) = sum(average)/count;
        EVM(p) = sqrt(sum(abs(average-average_single(p)*ones(1,length(average))).^2))/count;
        amplitude(p) = abs(average_single(p));
        Vpk(p) = abs(average_single(p));
        Pdbm(p) = 20*log10(Vpk(p))-5;
end
Pdbm = transpose(Pdbm);
Vpk = transpose(Vpk);
EVM = transpose(EVM);
SNR = 10*log10((Vpk.^2)./(EVM.^2));

%%
Pdbm_60mv_250k = [];
Pdbm_62mv_250k = [];
Pdbm_64mv_250k = [];
Pdbm_66mv_250k = [];
Pdbm_68mv_250k = [];
Pdbm_60mv_1M = [];

Vpk_60mv_250k = [];
Vpk_62mv_250k = [];
Vpk_64mv_250k = [];
Vpk_66mv_250k = [];
Vpk_68mv_250k = [];
Vpk_60mv_1M = [];

for a = 1:6
    Pdbm_60mv_250k(a) = Pdbm(6*a-5);
    Pdbm_62mv_250k(a) = Pdbm(6*a-4);
    Pdbm_64mv_250k(a) = Pdbm(6*a-3);
    Pdbm_66mv_250k(a) = Pdbm(6*a-2);
    Pdbm_68mv_250k(a) = Pdbm(6*a-1);
    Pdbm_60mv_1M(a)=   Pdbm(6*a);
 
    Vpk_60mv_250k(a) = Vpk(6*a-5);
    Vpk_62mv_250k(a) = Vpk(6*a-5);
    Vpk_64mv_250k(a) = Vpk(6*a-5);
    Vpk_66mv_250k(a) = Vpk(6*a-5);
    Vpk_68mv_250k(a) = Vpk(6*a-5);
    Vpk_60mv_1M(a)   = Vpk(6*a-5);
end
% distance = linspace(25,50,6);
% scatter(distance,Pdbm_60mv_250k)
% hold on
% scatter(distance,Pdbm_62mv_250k)
% scatter(distance,Pdbm_64mv_250k)
% scatter(distance,Pdbm_66mv_250k)
% scatter(distance,Pdbm_68mv_250k)
% % semilogy(Pdbm_70mv_250k)
% legend('60mv','62mv','64mv','66mv','68mv')
% axis([20 55 -120 -90])
scatter(real(all_outputs(31,:)),imag(all_outputs(31,:)))
hold on
scatter(real(average_single(31)),imag(average_single(31)))



r = [25:5:50];

figure

theory = plot(r_th, P_r ,'LineWidth', 2.5);
hold on
mv60 = scatter(r, Pdbm_60mv_250k, 200, 'g');
set(mv60, ...
    'Marker', 'o', ...
    'MarkerFaceColor', 'g');
hold on
mv62 = scatter(r, Pdbm_62mv_250k, 200, 'r');
set(mv62, ...
    'Marker', 'd', ...
    'MarkerFaceColor', 'r');
hold on
mv64 = scatter(r, Pdbm_64mv_250k, 200, 'b');
set(mv64, ...
    'Marker', 'v', ...
    'MarkerFaceColor', 'b');
hold on
mv66 = scatter(r, Pdbm_66mv_250k, 200, 'k');
set(mv66, ...
    'Marker', '^', ...
    'MarkerFaceColor', 'k');
hold on
mv68 = scatter(r, Pdbm_68mv_250k, 200, 'r');
set(mv68, ...
    'Marker', 's', ...
    'MarkerFaceColor', 'r');
hold on
mv601m = scatter(r, Pdbm_60mv_1M, 200, 'm');
set(mv601m, ...
    'Marker', 'p', ...
    'MarkerFaceColor', 'm');



axis([20, 55, -130, -90])


Legend = legend('Theoretical Link Budget', '$V_{DC}$ = 60.0 mV, $f_m$ = 250 kHz','$V_{DC}$ = 62.0 mV, $f_m$ = 250 kHz','$V_{DC}$ = 64.0 mV, $f_m$ = 250 kHz','$V_{DC}$ = 66.0 mV, $f_m$ = 250 kHz','$V_{DC}$ = 68.0 mV, $f_m$ = 250 kHz','$V_{DC}$ = 60.0 mV, $f_m$ = 1 MHz');
Xlabel = xlabel('Reader-to-Tag Distance (m)');
Ylabel = ylabel('Received Signal Strength (dBm)');
 
%Set Fonts and Axes Properties
set(gca, ...
     'FontName'  , 'Helvetica');
 
 set(gca, ...
     'FontSize'  , 22);
% 
set([Xlabel, Ylabel], ...
      'interpreter', 'latex', ...
      'FontSize', 28);
  
 set(Legend, ...
     'FontSize'  , 22, ...
     'interpreter', 'latex');
% 
%
set(gca, ...
  'FontSize'    , 28        , ...
  'Color'       , 'white'   , ...
  'Box'         , 'on'     , ...
  'TickDir'     , 'out'     , ...
  'TickLength'  , [.02 .02] , ...
  'XMinorTick'  , 'on'      , ...
  'YMinorTick'  , 'on'      , ...
  'XGrid'       , 'on'      , ...
  'YGrid'       , 'on'      , ...
  'XColor'      , [0 0 0], ...
  'YColor'      , [0 0 0], ...
  'YTick'       , -130:10:-90, ...
  'LineWidth'   , 1, ...
  'FontUnits', 'points', ...
   'FontWeight', 'normal', ...
   'FontSize', 28, ...
   'FontName', 'Times');

set(gcf, ...
  'Color'       , 'White'   );





% %export
% set(gcf, 'PaperPositionMode', 'auto');
% print -depsc2 ranges.eps
% close;
% 
% %postprocess
% fixPSlinestyle('ranges.eps', 'Fig_ranges.eps');







