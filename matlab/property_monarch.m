%Initialize simulation model of insect type

% ===============================
% Monarch, with morphological parameters taken mostly from Sun and Xiong,
% ================================

%Edit this file to change properties of the insect

% clear all
% clc
% close all

%Global variables passed to the various sub-modules
global MASS_body MASS_wings G RHO_f WingInertia...
    IN_body L_BODY...
    IN_wing1 IN_wing2 ...
    R_CG2O R_CG2ACb ...
    R_O2WG R_O2AC CHORD_w R_w AREA_W R2MS R1S AxisRot...
    BODY_ANGLE_FREE          ...
    STROKE_AMP STROKE_F STROKE_AVG ...
    ALPHAD ALPHAU ALPHA_PHASE Cn SCP K BETA_0 DELTAUR...
    U Re C_W U_w_aveQSM U_w_maxQSM

%% MISC Constant
%Accel of gravity (ft/s^2)
G = [0; 0; 9.81];  %in inertial frame
%Air density (slug/ft^3)
RHO_f=1.225;

% This is a parameter that allows us to turn off wing inertia in the EOMs
WingInertia = 1;

%% BODY PROPERTIES
%TOTAL mass of the insect body in kg (head, thorax, abdomen, but not wings)
MASS_body = 0.2577E-3; %kg

% Total Body Length
L_BODY = 37.58E-3;        %


R_CG2O   = [5 ; 5 ; 0] * 1E-3;   % [m], Vector from cg to wing root
% R_CG2O   = [0 ; 0 ; 0] * 1E-3;
R_CG2ACb = [0 ; 0 ; 0] ;         % initially assume that the Aero Center of body is colocated with cg

%% Based on notes in "Morphological_Parameters_vs01"
IN_body =  [3.42E-9 0 0;
            0 2.5E-8 0;
            0 0 2.5E-8] ;   % in body frame, kg-m^2

% This parameter is not used
%  BODY_ANGLE_FREE = 0.0*pi/180;  % per Sun & Xiong, JEB, 2005,

%% WING PROPERTIES
% Mass of ALL wings [kg]
MASS_wings = [2*0.0251E-3;0]; % total mass of all four wings
R_w        = 52.2E-3;         % [m] span of 1 wing
AREA_W     = 2.46E-3;         % Measured Area of the overlapped wings and body, calculated using wingImageProcessor
CHORD_w    = AREA_W /(2*R_w); % [m]
R2MS       = 0.456*R_w;       % [m] radius to 2nd moment of area (S), calculated using wingImageProcessor
R1S        = 1.0;             % dimless, from Ellington, 1984 (this gives a distribution of chord)

% nondim, chordwise location of axis of rotation in %chord; wing root is located on this line.
AxisRot = 0.15;				% Ennos uses 0.15 (1988) ; Xiong and Sun use 0.3 (2005)

% distance from wing root to wing CG, in wing frame, where x is spanwise, -y is increasing chord and z is normal to wing
R_O2WG = [.456*R_w; -0.2*CHORD_w; 0];
R_O2AC = [R2MS ; 0; 0]; %in wing frame % T

% Wing Inertia
% Ixy and Iyx set to zero

IN_wing1 = [1.31  0    0;
            0    1.31  0;
            0    0     3.62] * 1E-9;   % in wing frame

%% Flapping Kinematics

%     STROKE_AMP = 70;    % [degress]
STROKE_F = 10.1;    % [Hz]
%     STROKE_AVG = 0;     % [degrees]


BETA_0 = 0; %initial stroke plane angle, need in degrees because ExtractKinematicsV2 uses this term: RotMatrixBodytoRW(wingstates*pi/180);

%% Nondimensionalizing Parameters
%U = 4*STROKE_AMP*pi/180*STROKE_F*R2MS;
% Madhu, changed U to mean forward velocity magnitude of the trajectory
U = 1.38;  % [m/s] based on average thorax velocity magnitude

