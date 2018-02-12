clear all
close all
clc
set(0,'defaulttextinterpreter','latex')
set(0,'DefaultTextFontname', 'CMU Serif')
set(0,'DefaultAxesFontName', 'CMU Serif')


 
% inputs

%frequency
f = 5.8e9; %Hz
% miniumum required power at the receiver
Prmin = -110; %dBm
%Sensitivity of Receiver Calculation


%scenarios 
%Listed parameters are in the following order:
%Transmitted power ranges (dBm); Tx antenna gain (dBi); Rx antenna gain(dBi); tag antenna gain (dBi)

%Town house
PT = -9.65; Gtx = 1.8; Grx = 4.71; Gt = 1.8; r = 1:20;
%Football Field
%PT = -0.7; Gtx = 1.8; Grx = 12; Gt = 2.1; r = 1:50;
%Warehouse
%PT = 5.9; Gtx = 1.8; Grx = 18.3; Gt = 1.8; r = 10:100;
%Skyscraper
%PT = 15; Gtx = 1.8; Grx = 27; Gt = 2.1; r = 1:300;
%University Campus
%PT = 14.5; Gtx = 1.8; Grx = 30; Gt = 9; r = 10:700;
%Crop Field
%PT = 20; Gtx = 1.8; Grx = 35; Gt = 9; r = 10:1000;
%Airport
%PT = 23; Gtx = 1.8; Grx = 35; Gt = 10.7; r = 10:1500;
%City
%PT = 23; Gtx = 1.8; Grx = 35; Gt = 13.1; r = 10:2000;




 
%%
%gain of the reflection amplifier fitted to the experimental results

 
% p1 =      0.9694;
% p2 =       2.946;
% p3 =       -3.38;
% p4 =      -16.22;
% p5 =        24.2;
a1 =       37.18;%  (31.84, 42.52)
b1 =      -77.14;%  (-79.86, -74.42)
c1 =       10.52;%  (7.169, 13.86)
Gqtr = @(x) a1*exp(-((x-b1)/c1).^2); % x is incident tag power in dBm




%%

 
% Conversions to linear

lambda = 3e8/f;

PTLin = 10.^(PT/10); %mW
GtxLin = 10.^(Gtx/10); %lin
GrxLin = 10.^(Grx/10); %lin
GtLin = 10.^(Gt/10); %lin

Pt = PTLin*GtxLin*GtLin*(lambda./(4*pi*r)).^2; %power on tag (mW)

G = Gqtr(10*log10(Pt)); %dB
GLin = 10.^(G/10);

Pr = Pt.*GrxLin.*GtLin.*GLin.*(lambda./(4*pi*r)).^2;

figure(1)
semilogx(r, 10*log10(Pr), 'k:', 'LineWidth', 1.5);
hold on





%Football Field
PT = -10.67; Gtx = 13; Grx = 13; Gt = 1.8; r = 1:50;
 
% Conversions to linear

lambda = 3e8/f;

PTLin = 10.^(PT/10); %mW
GtxLin = 10.^(Gtx/10); %lin
GrxLin = 10.^(Grx/10); %lin
GtLin = 10.^(Gt/10); %lin

Pt = PTLin*GtxLin*GtLin*(lambda./(4*pi*r)).^2; %power on tag (mW)

G = Gqtr(10*log10(Pt)); %dB
GLin = 10.^(G/10);

Pr = Pt.*GrxLin.*GtLin.*GLin.*(lambda./(4*pi*r)).^2;

figure(1)
semilogx(r, 10*log10(Pr), 'b--', 'LineWidth', 1.5);
hold on





%Warehouse
PT = -10.67; Gtx = 18.1; Grx = 18.1; Gt = 1.8; r = 1:100;

 
% Conversions to linear

lambda = 3e8/f;

PTLin = 10.^(PT/10); %mW
GtxLin = 10.^(Gtx/10); %lin
GrxLin = 10.^(Grx/10); %lin
GtLin = 10.^(Gt/10); %lin

Pt = PTLin*GtxLin*GtLin*(lambda./(4*pi*r)).^2; %power on tag (mW)

G = Gqtr(10*log10(Pt)); %dB
GLin = 10.^(G/10);

Pr = Pt.*GrxLin.*GtLin.*GLin.*(lambda./(4*pi*r)).^2;

figure(1)
semilogx(r, 10*log10(Pr), 'r-.', 'LineWidth', 1.5);
hold on


%Skyscraper
PT = -10.67; Gtx = 27.67; Grx = 27.67; Gt = 1.8; r = 1:300;


lambda = 3e8/f;

PTLin = 10.^(PT/10); %mW
GtxLin = 10.^(Gtx/10); %lin
GrxLin = 10.^(Grx/10); %lin
GtLin = 10.^(Gt/10); %lin

Pt = PTLin*GtxLin*GtLin*(lambda./(4*pi*r)).^2; %power on tag (mW)

G = Gqtr(10*log10(Pt)); %dB
GLin = 10.^(G/10);

Pr = Pt.*GrxLin.*GtLin.*GLin.*(lambda./(4*pi*r)).^2;

figure(1)
semilogx(r, 10*log10(Pr), 'g-', 'LineWidth', 1.5);
hold on


%University Campus
PT = -10.67; Gtx = 35; Grx = 35; Gt = 1.8; r = 1:700;


PTLin = 10.^(PT/10); %mW
GtxLin = 10.^(Gtx/10); %lin
GrxLin = 10.^(Grx/10); %lin
GtLin = 10.^(Gt/10); %lin

Pt = PTLin*GtxLin*GtLin*(lambda./(4*pi*r)).^2; %power on tag (mW)

G = Gqtr(10*log10(Pt)); %dB
GLin = 10.^(G/10);

Pr = Pt.*GrxLin.*GtLin.*GLin.*(lambda./(4*pi*r)).^2;

figure(1)
semilogx(r, 10*log10(Pr), 'k:', 'LineWidth', 2.5);
hold on

%Crop Field
PT = -10.67; Gtx = 35; Grx = 35; Gt = 6; r = 1:1000;

PTLin = 10.^(PT/10); %mW
GtxLin = 10.^(Gtx/10); %lin
GrxLin = 10.^(Grx/10); %lin
GtLin = 10.^(Gt/10); %lin

Pt = PTLin*GtxLin*GtLin*(lambda./(4*pi*r)).^2; %power on tag (mW)

G = Gqtr(10*log10(Pt)); %dB
GLin = 10.^(G/10);

Pr = Pt.*GrxLin.*GtLin.*GLin.*(lambda./(4*pi*r)).^2;

figure(1)
semilogx(r, 10*log10(Pr), 'b--', 'LineWidth', 2.5);
hold on





%Airport
PT = -11.52; Gtx = 35; Grx = 35; Gt = 9; r = 1:1500;


% Conversions to linear

PTLin = 10.^(PT/10); %mW
GtxLin = 10.^(Gtx/10); %lin
GrxLin = 10.^(Grx/10); %lin
GtLin = 10.^(Gt/10); %lin

Pt = PTLin*GtxLin*GtLin*(lambda./(4*pi*r)).^2; %power on tag (mW)

G = Gqtr(10*log10(Pt)); %dB
GLin = 10.^(G/10);

Pr = Pt.*GrxLin.*GtLin.*GLin.*(lambda./(4*pi*r)).^2;

figure(1)
semilogx(r, 10*log10(Pr), 'r-.', 'LineWidth', 2.5);
hold on



%City
PT = -10.67; Gtx = 35; Grx = 35; Gt = 11; r = 1:2000;



% Conversions to linear

PTLin = 10.^(PT/10); %mW
GtxLin = 10.^(Gtx/10); %lin
GrxLin = 10.^(Grx/10); %lin
GtLin = 10.^(Gt/10); %lin

Pt = PTLin*GtxLin*GtLin*(lambda./(4*pi*r)).^2; %power on tag (mW)
G = Gqtr(10*log10(Pt)); %dB
GLin = 10.^(G/10);
Pr = Pt.*GrxLin.*GtLin.*GLin.*(lambda./(4*pi*r)).^2;

figure(1)
r = logspace(0,3.1);
Pt = PTLin*GtxLin*GtLin*(lambda./(4*pi*r)).^2; %power on tag (mW)
G = Gqtr(10*log10(Pt)); %dB
GLin = 10.^(G/10);
Pr = Pt.*GrxLin.*GtLin.*GLin.*(lambda./(4*pi*r)).^2;
semilogx(r, 10*log10(Pr), 'go-', 'LineWidth', 2.5);



Legend = legend('Town house','Football field','Warehouse','Skyscraper','University Campus','Crop field','Airport','City');
Xlabel = xlabel('Range $r$ (m)');
Ylabel = ylabel('Power on Reader $P_r$ (dBm)');

axis([0, 2000, -120, 0])


%Set Fonts and Axes Properties
set(gca, ...
    'FontName'  , 'Helvetica');

set(Legend, ...
    'FontName'  , 'AvantGarde');

set([Xlabel, Ylabel], ...
    'FontName'  , 'AvantGarde', ...
    'Color', 'k');

set(Legend, ...
     'FontSize'  , 22, ...
     'interpreter', 'latex');
     
 set([Xlabel, Ylabel], ...
     'FontSize'  , 28, ...
     'interpreter', 'latex');


       
set(gca, ...
  'FontSize'    , 26        , ...
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
  'YTick'       , [-120:10:0], ...
  'LineWidth'   , 1, ...
  'FontUnits', 'points', ...
   'FontWeight', 'normal', ...
   'FontSize', 22, ...
   'FontName', 'Times');

set(gcf, ...
  'Color'       , 'White'   );

set(gcf, ...
   'PaperPositionMode', 'auto');



