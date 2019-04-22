clear all;
close all;

MONARCH.rho = 1.225; % air density (kg/m^3)
MONARCH.l = 52.2E-3; % span of the right wing
MONARCH.S = 2.46E-3; % area of the right wing
MONARCH.c_bar = MONARCH.S / MONARCH.l; % mean chord

MONARCH.tilde_r_2 = 0.1; % non-dimensional radius of the second moment of wing area
MONARCH.tilde_r_3 = 0.1; % non-dimensional radius of the thrid moment of wing area (arbitrary value)
MONARCH.r_cp = MONARCH.l * MONARCH.tilde_r_3 / MONARCH.tilde_r_2;

MONARCH.tilde_v = 0.1; % non-dimensional virtual mass (arbitrary value)
MONARCH.tilde_r_v_1 =0.1; % non-dimensional radius of the first moment of wing volume (arbitrary value)
MONARCH.tilde_r_v_2 =0.1; % non-dimensional radius of the second moment of wing volume (arbitrary value)
MONARCH.r_rot = MONARCH.l * MONARCH.tilde_r_v_2^2 / MONARCH.tilde_r_v_1; % non-dimensional radius of the second moment of wing volume (arbitrary value)

save('morp_MONARCH');
