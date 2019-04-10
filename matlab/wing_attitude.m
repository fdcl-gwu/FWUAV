function [Q_R Q_L varargout]=wing_attitude(beta,Euler_R,varargin)
%wing_attitude compute wing attitude and wing angular velocity
%
%   [Q_R Q_L]=wing_attitude(beta, Euler) returns the attitudes of the right wing
%   and the left wing for a given stroke plane angle and a set of Euler angles, assuming that the
%   left wing is symmetric to the right wing
%
%   [Q_R Q_L]=wing_attitude(beta, Euler_R, Euler_L) returns the attitudes of the right wing
%   and the left wing 
%
%   [Q_R Q_L W_R W_L]=wing_attitude(beta, Euler_R, Euler_L, Euler_R_dot,
%   Euler_L_dot) returns the attitude and the angular velocity of both
%   wings for given Euler anges and their time-derivatives
%
%   [Q_R Q_L W_R W_L W_R_dot W_L_dot]=wing_attitude(beta, Euler_R, Euler_L, Euler_R_dot,
%   Euler_L_dot, Euler_R_ddot, Euler_L_ddot) returns the attitude, the angular velocity, and the angular acceleration of both
%   wings for given Euler anges and their time-derivatives

e1=[1 0 0]';
e2=[0 1 0]';
e3=[0 0 1]';

if nargin == 2
    Euler_L = Euler_R;
elseif nargin > 2
    Euler_L = varargin{1};
end

phi_R=Euler_R(1);
theta_R=Euler_R(2);
psi_R=Euler_R(3);
Q_R=expmso3(beta*e2)*expmso3(phi_R*e1)*expmso3(-psi_R*e3)*expmso3(theta_R*e2);

phi_L=Euler_L(1);
theta_L=Euler_L(2);
psi_L=Euler_L(3);
Q_L=expmso3(beta*e2)*expmso3(-phi_L*e1)*expmso3(psi_L*e3)*expmso3(theta_L*e2);

if nargin > 3
    Euler_R_dot = varargin{2};
    Euler_L_dot = varargin{3};

    T_R=[cos(psi_R)*cos(theta_R), 0, sin(theta_R);
        sin(psi_R), 1, 0;
        cos(psi_R)*sin(theta_R), 0, -cos(theta_R)];
    W_R=T_R*Euler_R_dot;
    
    T_L = [-cos(psi_L)*cos(theta_L), 0, -sin(theta_L);
        sin(psi_L), 1, 0;
        -cos(psi_L)*sin(theta_L), 0, cos(theta_L)];
    W_L=T_L*Euler_L_dot;
    varargout{1}=W_R;
    varargout{2}=W_L;
    
    if nargin > 5
        Euler_R_ddot = varargin{4};
        Euler_L_ddot = varargin{5};
        
        T_R=[cos(psi_R)*cos(theta_R), 0, sin(theta_R);
            sin(psi_R), 1, 0;
            cos(psi_R)*sin(theta_R), 0, -cos(theta_R)];        
        T_R_theta=[ -cos(psi_R)*sin(theta_R), 0, cos(theta_R);
            0, 0, 0;
            cos(psi_R)*cos(theta_R), 0, sin(theta_R)];
        T_R_psi=[ -cos(theta_R)*sin(psi_R), 0, 0;
            cos(psi_R), 0, 0;
            -sin(psi_R)*sin(theta_R), 0, 0];
        W_R_dot=T_R*Euler_R_ddot + T_R_theta*Euler_R_dot*Euler_R_dot(2) + T_R_psi*Euler_R_dot*Euler_R_dot(3);
        
        T_L = [-cos(psi_L)*cos(theta_L), 0, -sin(theta_L);
            sin(psi_L), 1, 0;
            -cos(psi_L)*sin(theta_L), 0, cos(theta_L)];
        T_L_theta = [cos(psi_L)*sin(theta_L), 0, -cos(theta_L);
            0, 0, 0;
            -cos(psi_L)*cos(theta_L), 0, -sin(theta_L)];
        T_L_psi = [sin(psi_L)*cos(theta_L), 0, 0;
            cos(psi_L), 0, 0;
            sin(psi_L)*sin(theta_L), 0, 0];
        W_L_dot=T_L*Euler_L_ddot + T_L_theta*Euler_L_dot*Euler_L_dot(2) + T_L_psi*Euler_L_dot*Euler_L_dot(3);
        varargout{3}=W_R_dot;
        varargout{4}=W_L_dot;        
    end
        
end
