clear all;
close all;
e1=[1 0 0]';
e2=[0 1 0]';
e3=[0 0 0]';

load VICON_data

k=1;
N=length(t);
R = zeros(3,3,N);
R_R = zeros(3,3,N);
R_L = zeros(3,3,N);
R_A = zeros(3,3,N);
Q_R = zeros(3,3,N);
Q_L = zeros(3,3,N);
Q_A = zeros(3,3,N);

for k=1:N
    %% right wing
    % fit a plane passing thorough T1 to RW1, RW2, RW3
    M_RR = (RW1(:,k)-T1(:,k))*(RW1(:,k)-T1(:,k))'...
        + (RW2(:,k)-T1(:,k))*(RW2(:,k)-T1(:,k))'...
        + (RW3(:,k)-T1(:,k))*(RW3(:,k)-T1(:,k))';
    [V D]=eig(M_RR);
    [~, I_min]=min(diag(D));
    n=V(:,I_min);
    tmp3 = cross(RW1(:,k)-T1(:,k), RW3(:,k)-T1(:,k));
    tmp3 = tmp3/norm(tmp3);
    if dot(tmp3, n) < 0
        n= -n;
    end
    
    R_R(:,3,k) = n;
    % project wing tip to the plane
    tmp2 = (eye(3)-n*n')*(RW2(:,k)-T1(:,k));
    R_R(:,2,k) = tmp2/norm(tmp2);
    R_R(:,1,k) = cross(R_R(:,2,k),R_R(:,3,k));
    
    
    %% left wing
    % fit a plane passing thorough T1 to LW1, LW2, LW3
    M_LL = (LW1(:,k)-T1(:,k))*(LW1(:,k)-T1(:,k))'...
        + (LW2(:,k)-T1(:,k))*(LW2(:,k)-T1(:,k))'...
        + (LW3(:,k)-T1(:,k))*(LW3(:,k)-T1(:,k))';
    [V D]=eig(M_LL);
    [~, I_min]=min(diag(D));
    n=V(:,I_min);
    tmp3 = cross(LW3(:,k)-T1(:,k), LW1(:,k)-T1(:,k));
    tmp3 = tmp3/norm(tmp3);
    if dot(tmp3, n) < 0
        n= -n;
    end
    
    R_L(:,3,k) = n;
    % project wing tip to the plane
    tmp2 = (eye(3)-n*n')*(T1(:,k)-LW2(:,k));
    R_L(:,2,k) = tmp2/norm(tmp2);
    R_L(:,1,k) = cross(R_L(:,2,k),R_L(:,3,k));
    
    %% thorax attitude
    R(:,1,k) = (T1(:,k) - T2(:,k))/norm(T1(:,k) - T2(:,k));
    R(:,2,k) = e2; % assume there is no body roll
    R(:,3,k) = cross( R(:,1,k), R(:,2,k));
    
    %% abdomen attitude
    tmpA = (A1(:,k) + A2(:,k))/2;
    R_A(:,1,k) = (T2(:,k) - tmpA)/norm(T2(:,k) - tmpA);
    R_A(:,2,k) = e2; % assume there is no body roll
    R_A(:,3,k) = cross( R_A(:,1,k), R_A(:,2,k));    
    
    %% relative attitude
    Q_R(:,:,k)= R(:,:,k)'*R_R(:,:,k);
    Q_L(:,:,k)= R(:,:,k)'*R_L(:,:,k);
    Q_A(:,:,k)= R(:,:,k)'*R_A(:,:,k);
    
    theta_B(k)=atan2(-R(3,1,k),R(1,1,k));
    theta_A(k)=atan2(-Q_A(3,1,k),Q_A(1,1,k));
end

save fit_VICON_data