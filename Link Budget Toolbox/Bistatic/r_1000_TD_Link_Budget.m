clear all
close all
clc

% inputs

%frequency
f = 5.8e9; %Hz
x = linspace(-100,-50,100);
% distance between reader and tag
r = 1000; %m

% miniumum required power at the receiver
Prmin = -110; %dBm

% minimum required SNR
SNRmin = 10; %dB

% noise figure
NF = 2; %dB

%Bandwidth of the receiver
BW = 1e6; %1 Hz

% Receiver Noise power
Pn = -111; %dBm

% Processing gain
Gpr = 10; %dB

%%

%Sensitivity of Receiver Calculation
N = 30;
% Transmitted power ranges
PT = linspace(-20, 30, N); %dBm
% Transmitting antenna gain
Gtx = linspace(1.8, 7.415, N); %dBi
% Receiving antenna gain
Grx = linspace(1, 35, N); %dBi
% Tag antenna gain
Gt = linspace(9, 13.07, N); %dBi

%%
%gain of the reflection amplifier fitted to the experimental results

% p1 =      0.9694;
% p2 =       2.946;
% p3 =       -3.38;
% p4 =      -16.22;
% p5 =        24.2;
p1 =    0.001327;%  (-0.008781, 0.01144)
p2 =      0.3925;%  (-2.5, 3.285)
p3 =       43.28;%  (-266.4, 353)
p4 =        2107;%  (-1.26e+04, 1.681e+04)
p5 =   3.821e+04;%  (-2.232e+05, 2.996e+05)
a1 =       34.63;%  (31.84, 42.52)
b1 =      -78.24;%  (-79.86, -74.42)
c1 =       93.67;%  (7.169, 13.86)
Gqtr = @(x) a1*exp(-((x-b1).^2/c1))+4; % x is incident tag power in dBm

% plot(x,Gqtr(x))
%%
temp_result = zeros(N^4,6);
lambda = 3e8/(f);%m
path = 10*log10((lambda./(4*pi.*r)).^2);
count = 1;
possible_rows = zeros(N^4,1);
for a = 1:length(PT)
    for b = 1:length(Gtx)
        for c = 1:length(Grx)
            for d = 1:length(Gt)
                Ptag = PT(a) + Gtx(b) + Gt(d) + path;
                G_tag = Gqtr(Ptag);
                P_backscattered = Ptag + G_tag + Gt(d) + Grx(c) + path;
                index = (N^3*a-N^3) + (N^2*b-N^2) + (N*c-N) + d;
                temp_result(index,4) = Gt(d);
                temp_result(index,3) = Grx(c);
                temp_result(index,2) = Gtx(b);
                temp_result(index,1) = PT(a);
                temp_result(index,5) = P_backscattered;
                temp_result(index,6) = Ptag;
                
                if((Ptag > -80 && Ptag < -60) && (P_backscattered > Prmin))
                    possible_rows(count) = index;
                    count = count + 1;
                end
            end
        end
    end
end
possible_rows = possible_rows(possible_rows~=0);
result = zeros(length(possible_rows),6);
for a = 1:(length(possible_rows))
    result(a,:) = temp_result( possible_rows(a)   ,:);
end

[a,min_tag_row] = min(result(:,4));
min_Gtag_indices = find(result(:,4) == min(result(:,4)));
min_Gtag_results = zeros(length(min_Gtag_indices),6);
for a = 1:length(min_Gtag_indices)
    min_Gtag_results(a,:) = result(min_Gtag_indices(a),:);
end

%TAKE 300
r = 1:1:1000;
path_2 = 10*log10((lambda./(4*pi.*r)).^2);

a = 300;
PT = min_Gtag_results(a,1);Gtx = min_Gtag_results(a,2); Grx = min_Gtag_results(a,3); Gt = min_Gtag_results(a,4);
Ptag = PT + Gtx + Gt + path_2;
G_tag = Gqtr(Ptag);
P_backscattered = Ptag + G_tag + Gt + Grx + path_2;

semilogx(r,P_backscattered)