clear all;
close all;

%% polynomial fitting from Dr. Kang 
%% all of the units are in centimeters

span_wholewing=[0
    0.4475
    0.895
    1.3425
    1.79
    2.2375
    2.685
    3.1325
    3.58
    4.0275
    4.475
    4.9225
    5.37
    5.8175
    6.0997
    ];

chord_wholewing=[0.172
    0.6961
    1.326
    1.9723
    2.6186
    3.2319
    3.6628
    3.5136
    3.149
    2.6849
    1.8895
    1.5414
    1.2762
    0.746
    0.3319
    ];

cr_poly = polyfit(span_wholewing,chord_wholewing,6);
 
% definite integral of a polynomial over an interval
polyint_def=@(poly, interval) diff(polyval(polyint(poly), interval));

l=max(span_wholewing) % right wing length
S = polyint_def(cr_poly, [0 l]) % right wing area
AR = l^2/S
c_bar = S / l

cr_bar_poly = cr_poly.*(l.^(length(cr_poly)-1:-1:0)) /c_bar;

S_1 = polyint_def( conv(cr_poly,[1 0]), [0 l])
tilde_r_1 = S_1/S/l
disp('err');
disp(S_1/S/l- polyint_def( conv(cr_bar_poly,[1 0]), [0 1]))

S_2 = polyint_def( conv(cr_poly,[1 0 0]), [0 l])
tilde_r_2 = (S_2/S/l^2)^(1/2)
disp(tilde_r_2 - (polyint_def( conv(cr_bar_poly,[1 0 0]) ,[0 1]))^(1/2))

S_3 = polyint_def( conv(cr_poly,[1 0 0 0]), [0 l])
tilde_r_3 = (S_3/S/l^3)^(1/3)
disp(tilde_r_3 - (polyint_def( conv(cr_bar_poly,[1 0 0 0]) ,[0 1]))^(1/3))

tilde_v = polyint_def( conv(cr_bar_poly,cr_bar_poly) ,[0 1])
tilde_r_v_1 = 1/tilde_v* polyint_def( conv(conv(cr_bar_poly,cr_bar_poly),[1 0]), [0 1])
disp('err');
disp(polyint_def( conv(conv(cr_poly,cr_poly),[1 0]), [0 l]) /l^2 /c_bar^2 /tilde_v - tilde_r_v_1)
tilde_r_v_2 = (1/tilde_v* polyint_def( conv(conv(cr_bar_poly,cr_bar_poly),[1 0 0]), [0 1]))^(1/2)
disp('err');
disp(polyint_def( conv(conv(cr_poly,cr_poly),[1 0 0]), [0 l]) /l^3 /c_bar^2 /tilde_v - tilde_r_v_2^2)

%% plot chord
figure;
N=501;
r=linspace(0,l,N);
cr=polyval(cr_poly,r);
r_bar=linspace(0,1,N);
cr_bar=polyval(cr_bar_poly*c_bar,r_bar);
plot(r,cr,'r',span_wholewing,chord_wholewing,'ro');
hold on;
plot(r,cr_bar,'b:');
xlabel('$r\;(\mathrm{cm})$','interpreter','latex');
ylabel('$c(r)\;(\mathrm{cm})$','interpreter','latex');

%% save the computed morphological parameters into a structure variable

MONARCH.rho = 1.225; % air density (kg/m^3)
MONARCH.l = l*1E-2; % span of the right wing (m)
MONARCH.S = S*1E-4; % area of the right wing (m^2)
MONARCH.c_bar = MONARCH.S / MONARCH.l; % mean chord (m)
MONARCH.AR = AR;

MONARCH.tilde_r_1 = tilde_r_1; % non-dimensional radius of the first moment of wing area
MONARCH.tilde_r_2 = tilde_r_2; % non-dimensional radius of the second moment of wing area
MONARCH.tilde_r_3 = tilde_r_3; % non-dimensional radius of the thrid moment of wing area 
MONARCH.r_cp = MONARCH.l * MONARCH.tilde_r_3 / MONARCH.tilde_r_2; % CP of the translational force (lift and drag) (m)

MONARCH.tilde_v = tilde_v; % non-dimensional virtual mass 
MONARCH.tilde_r_v_1 = tilde_r_v_1; % non-dimensional radius of the first moment of wing volume 
MONARCH.tilde_r_v_2 = tilde_r_v_2; % non-dimensional radius of the second moment of wing volume 
MONARCH.r_rot = MONARCH.l * MONARCH.tilde_r_v_2^2 / MONARCH.tilde_r_v_1; % CP of the rotational force (m)
disp(MONARCH);

save('morp_MONARCH','MONARCH');
