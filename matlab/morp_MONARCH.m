clear all;
close all;

%% polynomial fitting from Dr. Kang 
%% all of the units are in centimeters
scale = 1e-2; %(certain data in centimeters)

DATA=[0	0.172	0.0329	-0.1391
0.4475	0.6961	0.2329	-0.4632
0.895	1.326	0.3974	-0.9286
1.3425	1.9723	0.5329	-1.4394
1.79	2.6186	0.6595	-1.9591
2.2375	3.2319	0.765	-2.4669
2.685	3.6628	0.8822	-2.7806
3.1325	3.5136	0.9702	-2.5434
3.58	3.149	1.0346	-2.1144
4.0275	2.6849	1.0727	-1.6122
4.475	1.8895	1.0815	-0.808
4.9225	1.5414	1.0845	-0.4569
5.37	1.2762	1.0933	-0.1829
5.8175	0.746	0.9614	0.2154
6.0997	0.3319	0.7621	0.4302];

m_B = 1.485e-4;
h_B = 1.463e-2;
w_B = 5.680e-3;
m_A = 1.092e-4;
h_A = 1.738e-2;
w_A = 3.710e-3;
m_R = 2.510e-5;

span_wholewing=DATA(:,1);
chord_wholewing=DATA(:,2);
chord_LE=DATA(:,3);
chord_TE=DATA(:,4);
cr_poly = polyfit(span_wholewing,chord_wholewing,6);
cr_LE_poly = polyfit(span_wholewing,chord_LE,6);
cr_TE_poly = polyfit(span_wholewing,chord_TE,6);

%% compute area / volume moment
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

%% compute momnet of inertia

J_R_xx = m_R/S * polyint_def( conv(cr_LE_poly-cr_TE_poly,[1 0 0]), [0 l]) * scale^2;
cr_LE_poly_pow2 = conv(cr_LE_poly,cr_LE_poly);
cr_TE_poly_pow2 = conv(cr_TE_poly,cr_TE_poly);
J_R_xy = - m_R/S * polyint_def( 0.5* conv(cr_LE_poly_pow2 - cr_TE_poly_pow2, [1 0]), [0 l])* scale^2;
cr_LE_poly_pow3 = conv(cr_LE_poly_pow2,cr_LE_poly);
cr_TE_poly_pow3 = conv(cr_TE_poly_pow2,cr_TE_poly);
J_R_yy = m_R/S * polyint_def( 1/3* (cr_LE_poly_pow3 - cr_TE_poly_pow3), [0 l])* scale^2;
J_R_zz = J_R_xx + J_R_yy;

nu_R_x = 1/S * polyint_def( 0.5* (cr_LE_poly_pow2 - cr_TE_poly_pow2), [0 l])* scale
nu_R_y = 1/S * polyint_def( 0.5* conv(cr_LE_poly - cr_TE_poly, [1 0]), [0 l])* scale
nu_R = [nu_R_x nu_R_y 0]';
nu_L = [nu_R_x -nu_R_y 0]';
nu_A = [-h_A/2, 0, 0]';

J_R = [J_R_xx J_R_xy 0;
    J_R_xy J_R_yy 0
    0 0 J_R_zz];

J_L = [J_R_xx -J_R_xy 0;
    -J_R_xy J_R_yy 0
    0 0 J_R_zz];

J_B_xx = 1/8*m_B*w_B^2;
J_B_yy = m_B*(1/16*w_B^2 + 1/12*h_B^2);
J_B_zz = J_B_yy;
J_B = diag([J_B_xx J_B_yy J_B_zz]);

J_A_xx = 1/8*m_A*w_A^2;
J_A_yy = m_A*(1/16*w_A^2 + 1/3*h_A^2);
J_A_zz = J_A_yy;
J_A = diag([J_A_xx J_A_yy J_A_zz]);

% N=5000;
% ys=linspace(0,l,N);
% Dy=ys(2)-ys(1);
% tmp=zeros(3,1);
% for k=1:N        
%     c_LE=polyval(cr_LE_poly,ys(k));
%     c_TE=polyval(cr_TE_poly,ys(k));
%     y=ys(k);
%     tmp(1)=tmp(1)+y^2*(c_LE-c_TE)*Dy;
%     tmp(2)=tmp(2)+-0.5*y*(c_LE^2-c_TE^2)*Dy;
%     tmp(3)=tmp(3)+1/3*(c_LE^3-c_TE^3)*Dy;
% end
% tmp*m_R*scale^2/S  

%% plot chord
figure;
N=501;
r=linspace(0,l,N);
cr=polyval(cr_poly,r);
c_LE=polyval(cr_LE_poly,r);
c_TE=polyval(cr_TE_poly,r);
r_bar=linspace(0,1,N);
cr_bar=polyval(cr_bar_poly*c_bar,r_bar);
plot(r,cr,'r',span_wholewing,chord_wholewing,'r.');
hold on;
plot(r,cr_bar,'b:');
xlabel('$r\;(\mathrm{cm})$','interpreter','latex');
ylabel('$c\;(\mathrm{cm})$','interpreter','latex');

figure;
plot(r,c_LE,'r', span_wholewing,chord_LE,'r.');
hold on;
plot(r,c_TE,'b', span_wholewing,chord_TE,'b.');
xlabel('$r\;(\mathrm{cm})$','interpreter','latex');
ylabel('$c_{LE}, c_{TE}\;(\mathrm{cm})$','interpreter','latex');
axis equal;

bool_print=1;
if bool_print
    figure(1);print('MONARCH_c','-depsc');
    figure(2);print('MONARCH_c_LE_TE','-depsc');
    !cp *.eps ../doc/Figs
end

%% save the computed morphological parameters into a structure variable

MONARCH.scale = scale;
MONARCH.rho = 1.225; % air density (kg/m^3)
MONARCH.g = 9.81; % air density (kg/m^3)
MONARCH.l = l*scale; % span of the right wing (m)
MONARCH.S = S*scale^2; % area of the right wing (m^2)
MONARCH.c_bar = MONARCH.S / MONARCH.l; % mean chord (m)
MONARCH.AR = AR;

MONARCH.m_B = m_B; % mass of body/thorax
MONARCH.h_B = h_B; % height of body/thorax
MONARCH.w_B = w_B; % width of body/thorax
MONARCH.m_A = m_A; % mass of abdomen
MONARCH.h_A = h_A; % height of abdomen
MONARCH.w_A = w_A; % width of abdomen
MONARCH.m_R = m_R; % mass of right wing
MONARCH.m_L = m_R; % mass of left wing
MONARCH.m = m_B + m_A + 2*m_R;

MONARCH.J_B = J_B;
MONARCH.J_R = J_R;
MONARCH.J_A = J_A;
MONARCH.J_L = J_L;

MONARCH.nu_R = nu_R;
MONARCH.nu_A = nu_A;
MONARCH.nu_L = nu_L;

MONARCH.mu_R = [0, w_B/2, 0]';
MONARCH.mu_L = [0, -w_B/2, 0]';
MONARCH.mu_A = [-h_B/2, 0, 0]';

MONARCH.tilde_r_1 = tilde_r_1; % non-dimensional radius of the first moment of wing area
MONARCH.tilde_r_2 = tilde_r_2; % non-dimensional radius of the second moment of wing area
MONARCH.tilde_r_3 = tilde_r_3; % non-dimensional radius of the thrid moment of wing area 
MONARCH.r_cp = MONARCH.l * MONARCH.tilde_r_3 / MONARCH.tilde_r_2; % CP of the translational force (lift and drag) (m)

MONARCH.tilde_v = tilde_v; % non-dimensional virtual mass 
MONARCH.tilde_r_v_1 = tilde_r_v_1; % non-dimensional radius of the first moment of wing volume 
MONARCH.tilde_r_v_2 = tilde_r_v_2; % non-dimensional radius of the second moment of wing volume 
MONARCH.r_rot = MONARCH.l * MONARCH.tilde_r_v_2^2 / MONARCH.tilde_r_v_1; % CP of the rotational force (m)

MONARCH.cr_poly = cr_poly;
MONARCH.cr_LE_poly = cr_LE_poly;
MONARCH.cr_TE_poly = cr_TE_poly;
disp(MONARCH);

save('morp_MONARCH','MONARCH');
