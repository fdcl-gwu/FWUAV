function ani_sim_QS
%generate animation from sim_QS output
clear all;
close all;
global e1 e2 e3
e1=[1 0 0]';
e2=[0 1 0]';
e3=[0 0 1]';
II=eye(3);

load STLRead/fv_monarch;

load('morp_MONARCH');
bool_video=true;
file_to_save = 'sim_QS_xR_hover.avi';

load('sim_QS_xR_hover.mat','x','R','Q_R','Q_L','Q_A', 'N','F_R','F_L','des');

%% generate the initial object when k=1
k=1;

h_fig=figure('color','w','Unit','normalized','Position',[0.2 0.2 0.6 0.6]);
[h_body, h_wr, h_wl, h_ab]=patch_monarch(fv_body, fv_wr, fv_wl, fv_abdomen, x(:,k), R(:,:,k), Q_R(:,:,k), Q_L(:,:,k), Q_A(:,:,k), [0 0 0]);
axis([min(x(1,:)), max(x(1,:)), min(x(2,:)), max(x(2,:)), min(x(3,:)), max(x(3,:))] + ...
    0.05*[-1, 1, -1, 1, -1, 1]);
D_R=zeros(3,1);
D_L=zeros(3,1);
M_R=zeros(3,1);
M_L=zeros(3,1);
F_rot_R=zeros(3,1); 
F_rot_L=zeros(3,1);
M_rot_R=zeros(3,1);
M_rot_L=zeros(3,1);

[h_F_R, h_F_L]=patch_force(x(:,k), R(:,:,k), Q_R(:,:,k), Q_L(:,:,k),F_R(:,k),F_L(:,k),D_R,D_L,M_R,M_L,F_rot_R,F_rot_L,M_rot_R,M_rot_L,[0 0 0 0 1]);
axis([min(x(1,:)), max(x(1,:)), min(x(2,:)), max(x(2,:)), min(x(3,:)), max(x(3,:))] + ...
    0.05*[-1, 1, -1, 1, -1, 1]);
h_x=line(x(1,1:k),x(2,1:k),x(3,1:k),'color','k','LineWidth',2);
% plot3(des.x_fit_t(1,:),des.x_fit_t(2,:),des.x_fit_t(3,:),'b','LineWidth',2);

%% animation

if bool_video
    vidObj = VideoWriter(file_to_save);
    vidObj.FrameRate = 30; % To decrease the speed of video
    open(vidObj);
end

for k=floor(linspace(1, N, round(N/2)))
    disp(k);
    update_monarch(h_fig,[h_body, h_wr, h_wl, h_ab], fv_body, fv_wr, fv_wl, fv_abdomen, x(:,k), R(:,:,k), Q_R(:,:,k), Q_L(:,:,k), Q_A(:,:,k));
    update_force(h_F_R, h_F_L, x(:,k), R(:,:,k), Q_R(:,:,k), Q_L(:,:,k), F_R(:,k), F_L(:,k),D_R,D_L,M_R,M_L,F_rot_R,F_rot_L,M_rot_R,M_rot_L,[0 0 0 0 1]);
    drawnow;
    set(h_x,'XData',x(1,1:k),'YData', x(2,1:k), 'ZData', x(3,1:k));
    
    axis([min(x(1,:)), max(x(1,:)), min(x(2,:)), max(x(2,:)), min(x(3,:)), max(x(3,:))] + ...
    0.05*[-1, 1, -1, 1, -1, 1]);
    if bool_video
        ax = gca;
        ax.Units = 'pixels';
        pos = ax.Position;
        rect = [pos(3)*0.375, pos(4)*0.125, pos(3)*0.6, pos(4)]; % Region of the frame to capture
        writeVideo(vidObj,getframe(gcf,rect));
    end
end

if bool_video
    close(vidObj);
end

%% save
filename='ani_sim_QS';
varData = whos;
saveIndex = cellfun(@isempty, regexp({varData.class}, 'matlab.(graphics|ui)'));
saveVars = {varData(saveIndex).name};
save(filename,saveVars{:});
evalin('base',['load ' filename]);

end



function [h_F_R h_F_L]=patch_force(x,R,Q_R,Q_L, L_R,L_L,D_R,D_L,M_R,M_L,F_rot_R,F_rot_L,M_rot_R,M_rot_L,varargin)
global scale_Force
e2=[0 1 0]';
lwidth=1;
alength=0.01;
awidth=alength*tand(10);
scale_Force=0.8e1;
wingspan=0.06;

if nargin < 15
    bool=ones(5,1);
else
    bool = varargin{1};
end

if bool(1)
    acolor=[1 0 0]';
    x_cp = x + R*Q_R*0.6*wingspan*e2;
    h_L_R=patch_arrow(x_cp, x_cp + scale_Force*R*Q_R*L_R, acolor, lwidth, alength, awidth);
    x_cp = x - R*Q_L*0.6*wingspan*e2;
    h_L_L=patch_arrow(x_cp, x_cp + scale_Force*R*Q_L*L_L, acolor, lwidth, alength, awidth);
else
    h_L_R=[0 0]';
    h_L_L=[0 0]';
end


if bool(2)
    acolor=[0 0 1]';
    x_cp = x + R*Q_R*0.6*wingspan*e2;
    h_D_R=patch_arrow(x_cp, x_cp + scale_Force*R*Q_R*D_R, acolor, lwidth, alength, awidth);
    x_cp = x - R*Q_L*0.6*wingspan*e2;
    h_D_L=patch_arrow(x_cp, x_cp + scale_Force*R*Q_L*D_L, acolor, lwidth, alength, awidth);
else
    h_D_R=[0 0]';
    h_D_L=[0 0]';
end


if bool(3)
    acolor=[0.72,0.27,1.00]';
    x_cp = x + R*Q_R*0.6*wingspan*e2;
    h_LD_R=patch_arrow(x_cp, x_cp + scale_Force*R*Q_R*(D_R+L_R), acolor, lwidth, alength, awidth);
    x_cp = x - R*Q_L*0.6*wingspan*e2;
    h_LD_L=patch_arrow(x_cp, x_cp + scale_Force*R*Q_L*(D_L+L_L), acolor, lwidth, alength, awidth);
else
    h_LD_R=[0 0]';
    h_LD_L=[0 0]';
end


if bool(4)
    acolor=[0 1 0]';
    x_cp = x + R*Q_R*0.6*wingspan*e2;
    h_F_rot_R=patch_arrow(x_cp, x_cp + scale_Force*R*Q_R*(F_rot_R), acolor, lwidth, alength, awidth);
    x_cp = x - R*Q_L*0.6*wingspan*e2;
    h_F_rot_L=patch_arrow(x_cp, x_cp + scale_Force*R*Q_L*(F_rot_L), acolor, lwidth, alength, awidth);
else
    h_F_rot_R=[0 0]';
    h_F_rot_L=[0 0]';
end


if bool(5)
    acolor=[0 0 0]';
    x_cp = x + R*Q_R*0.6*wingspan*e2;
    h_F_tot_R=patch_arrow(x_cp, x_cp + scale_Force*R*Q_R*(D_R+L_R+F_rot_R), acolor, lwidth, alength, awidth);
    x_cp = x - R*Q_L*0.6*wingspan*e2;
    h_F_tot_L=patch_arrow(x_cp, x_cp + scale_Force*R*Q_R*(D_L+L_L+F_rot_L), acolor, lwidth, alength, awidth);
else
    h_F_tot_R=[0 0]';
    h_F_tot_L=[0 0]';
end

h_F_R=[h_L_R, h_D_R, h_LD_R, h_F_rot_R, h_F_tot_R];
h_F_L=[h_L_L, h_D_L, h_LD_L, h_F_rot_L, h_F_tot_L];
end

function update_force(h_F_R,h_F_L, x,R,Q_R,Q_L, L_R,L_L,D_R,D_L,M_R,M_L,F_rot_R,F_rot_L,M_rot_R,M_rot_L,varargin)
global scale_Force
e1=[1 0 0]';
e2=[0 1 0]';
alength=0.01;
awidth=alength*tand(10);
wingspan=0.06;

if nargin < 17
    bool=ones(5,1);
else
    bool = varargin{1};
end

if bool(1)
    x_cp = x + R*Q_R*0.6*wingspan*e2;
    update_arrow(h_F_R(:,1), x_cp, x_cp + scale_Force*R*Q_R*L_R, alength, awidth);
    x_cp = x - R*Q_L*0.6*wingspan*e2;
    update_arrow(h_F_L(:,1), x_cp, x_cp + scale_Force*R*Q_L*L_L, alength, awidth);
end
if bool(2)
    x_cp = x + R*Q_R*0.6*wingspan*e2;
    update_arrow(h_F_R(:,2), x_cp, x_cp + scale_Force*R*Q_R*D_R, alength, awidth);
    x_cp = x - R*Q_L*0.6*wingspan*e2;
    update_arrow(h_F_L(:,2), x_cp, x_cp + scale_Force*R*Q_L*D_L, alength, awidth);
end
if bool(3)
    x_cp = x + R*Q_R*0.6*wingspan*e2;
    update_arrow(h_F_R(:,3), x_cp, x_cp + scale_Force*R*Q_R*(L_R+D_R), alength, awidth);
    x_cp = x - R*Q_L*0.6*wingspan*e2;
    update_arrow(h_F_L(:,3), x_cp, x_cp + scale_Force*R*Q_L*(L_L+D_L), alength, awidth);
end
if bool(4)
    x_cp = x + R*Q_R*0.6*wingspan*e2;
    update_arrow(h_F_R(:,4), x_cp, x_cp + scale_Force*R*Q_R*(F_rot_R), alength, awidth);
    x_cp = x - R*Q_L*0.6*wingspan*e2;
    update_arrow(h_F_L(:,4), x_cp, x_cp + scale_Force*R*Q_L*(F_rot_L), alength, awidth);
end
if bool(5)
    x_cp = x + R*Q_R*wingspan*(0.6*e2+0.1*e1);
    update_arrow(h_F_R(:,5), x_cp, x_cp + scale_Force*R*Q_R*(D_R+L_R+F_rot_R), alength, awidth);
    x_cp = x + R*Q_L*wingspan*(-0.6*e2+0.1*e1);
    update_arrow(h_F_L(:,5), x_cp, x_cp + scale_Force*R*Q_L*(D_L+L_L+F_rot_L), alength, awidth);
end

end


function [v_body, v_wr, v_wl, v_ab]=compute_vertices(fv_body, fv_wr, fv_wl, fv_ab, x, R, Q_R, Q_L, Q_A)
v_body=zeros(size(fv_body.vertices));
v_wr=zeros(size(fv_wr.vertices));
v_wl=zeros(size(fv_wl.vertices));
v_ab=zeros(size(fv_ab.vertices));

for i=1:length(v_body)
    v_body(i,:)=(R*fv_body.vertices(i,:)'+x)';
end
for i=1:length(v_wr)
    v_wr(i,:)=(R*Q_R*(fv_wr.vertices(i,:)')+x)';
end
for i=1:length(v_wl)
    v_wl(i,:)=(R*Q_L*(fv_wl.vertices(i,:)')+x)';
end
ab_x_shift=-0.01179;
e1=[1 0 0]';
for i=1:length(v_ab)
    v_wl(i,:)=(R*(Q_A*(fv_ab.vertices(i,:)'-ab_x_shift*e1)+ab_x_shift*e1)+x)';
end

end

function [h_body, h_wr, h_wl, h_ab, h_FB, h_FR, f_FL]=patch_monarch(fv_body, fv_wr, fv_wl, fv_ab, x, R, Q_R, Q_L, Q_A, varargin)
if nargin < 10
    bool_FB=true; % show the body-fixed frame
    bool_FR=false; % show the right wing frame
    bool_FL=false % show the left wing frame
else
    bool_FB=varargin{1}(1);
    bool_FR=varargin{1}(2);
    bool_FL=varargin{1}(3);
end

II=eye(3);

[v_body, v_wr, v_wl, v_ab]=compute_vertices(fv_body, fv_wr, fv_wl, fv_ab, x, R, Q_R, Q_L, Q_A);

my_face_color=[0.8 0.8 1.0];
h_body=patch('faces', fv_body.faces, 'vertices', v_body);
hold on;
h_wr=patch('faces', fv_wr.faces, 'vertices', v_wr);
h_wl=patch('faces', fv_wl.faces, 'vertices', v_wl);
h_ab=patch('faces', fv_ab.faces, 'vertices', v_ab);
set([h_body, h_wr, h_wl, h_ab], ...
    'FaceColor', my_face_color, ...
    'EdgeColor',       'none',        ...
    'FaceLighting',    'gouraud',     ...
    'AmbientStrength', 0.15);
alpha([h_body, h_wr, h_wl],0.8);
axis('image');

view(180+10,30);
%view(-90,80);
% view(250,30);

set(gca,'Zdir','reverse','YDir','reverse');
axis([min(x(1,:)), max(x(1,:)), min(x(2,:)), max(x(2,:)), min(x(3,:)), max(x(3,:))] + ...
    0.05*[-1, 1, -1, 1, -1, 1]);
%axis([0 1.0 -0.1 0.1 -0.3 0.1]);
axis equal;
camlight headlight;
material dull;
set(gca,'visible','off');

lwidth=1;
alength=12;
awidth=alength*tand(10);
if bool_FB
    acolor=[0 0 1];
    for i=1:3
        h_FB(:,i)=patch_arrow(x', x'+140*(R*II(:,i))', acolor, lwidth, alength, awidth);
    end
end

acolor=[1 0 0];
mu_R=[12 8 -3]';
if bool_FR
    for i=1:3
        h_FR(:,i)=patch_arrow(R*mu_R, R*mu_R+ 120*R*Q_R*II(:,i), acolor, lwidth, alength, awidth);
    end
end

acolor=[1 0 0];
mu_L=[12 -8 -3]';
if bool_FL
    for i=1:3
        h_FL(:,i)=patch_arrow(R*mu_L, R*mu_L+ 120*R*Q_R*II(:,i), acolor, lwidth, alength, awidth);
    end
end

% xlabel('$x$','interpreter','latex');
% ylabel('$y$','interpreter','latex');
% zlabel('$z$','interpreter','latex');

end

function update_monarch(h_fig, h_objects, fv_body, fv_wr, fv_wl, fv_ab, x, R, Q_R, Q_L, Q_A)
h_body=h_objects(1);
h_wr=h_objects(2);
h_wl=h_objects(3);
h_ab=h_objects(4);
[v_body, v_wr, v_wl, v_ab]=compute_vertices(fv_body, fv_wr, fv_wl, fv_ab, x, R, Q_R, Q_L, Q_A);

set(h_body,'faces',fv_body.faces,'vertices',v_body);
set(h_wr,'faces',fv_wr.faces,'vertices',v_wr);
set(h_wl,'faces',fv_wl.faces,'vertices',v_wl);
set(h_ab,'faces',fv_ab.faces,'vertices',v_ab);
figure(h_fig);

drawnow;
end

function fig_note(fv_body, fv_wr, fv_wl,bool_print)
global e1 e2 e3
II=eye(3);

%% body-fixed frame
h_FB=figure('color','w');
x=100*[-1.5 1 -1]';
R=expmhat(-pi/12*e3)*expmhat(pi/6*e2)*expmhat(-pi/12*e1);
Q_R=eye(3);
Q_L=eye(3);
[h_body, h_wr, h_wl]=patch_monarch(fv_body, fv_wr, fv_wl, x, R, Q_R, Q_L,[1,0,0]);
alpha([h_body, h_wr, h_wl],0.25);

% show the inertial frame
acolor=[0 0 0];
lwidth=1;
alength=12;
awidth=alength*tand(10);
for i=1:3
    patch_arrow([0 0 0], 50*(II(:,i))', acolor, lwidth, alength, awidth);
end
axis auto;

%% right wing frame
h_FR=figure('color','w');
x=[0 0 0]';
R=eye(3);
Q_R=eye(3);
Q_L=eye(3);
[h_body, h_wr, h_wl]=patch_monarch(fv_body, fv_wr, fv_wl, x, R, Q_R, Q_L,[1 1 0]);
alpha([h_body, h_wr, h_wl],0.25);
view(-90,90);

axis auto;

%% left wing frame
h_FL=figure('color','w');
x=[0 0 0]';
R=eye(3);
Q_R=eye(3);
Q_L=eye(3);
[h_body, h_wr, h_wl]=patch_monarch(fv_body, fv_wr, fv_wl, x, R, Q_R, Q_L,[1 0 1]);
alpha([h_body, h_wr, h_wl],0.25);
view(-90,90);

axis auto;

%% stroke plane
beta=15*pi/180;
h_FS=figure('color','w');
mu_R=[12 8 -3]';
x=[0 0 0]';
R=expmhat(pi/6*e2);
patch_circle(R*mu_R, R*expmhat(beta*e2)*e1, 170);

[Q_R Q_L]=wing_attitude(beta,[-pi/6,0,0]);
[h_body, h_wr, h_wl]=patch_monarch(fv_body, fv_wr, fv_wl, x, R, Q_R, Q_L, [1 0 0]);
alpha([h_body, h_wr, h_wl],0.8);
hold on;
view(370,-0);
camlight;

x0=R*[mu_R(1) 0 0]'; % origin of the unit vector normal to the stroke plane
for i=1:3
    patch_arrow(x0, x0 + 100*R*expmhat(beta*e2)*II(:,i), [0 1 0], lwidth, alength, awidth);
end
plot_arc(x0, R*expmhat(beta*e2)*e1, R*e1, 70);

%% flapping angle
h_FR_phi=figure('color','w');
x=[0 0 0]';
R=expmhat(pi/6*e2);
beta=15*pi/180;

[Q_R Q_L]=wing_attitude(beta,[pi/3,0,0]);
[h_body, h_wr, h_wl]=patch_monarch(fv_body, fv_wr, fv_wl, x, R, Q_R, Q_L, [0 0 0]);
alpha([h_body, h_wr, h_wl],0.9);
hold on;
view(270,-5);
camlight;

axis auto;

x0=R*mu_R;
for i=2
    patch_arrow(x0, x0 + 100*R*Q_R*II(:,i), [1 0 0], lwidth, alength, awidth);
    patch_arrow(x0, x0 + 100*R*expmhat(beta*e2)*II(:,i), [0 1 0], lwidth, alength, awidth);
    plot_arc(x0, R*Q_R*II(:,i), R*expmhat(beta*e2)*II(:,i), 40);
end

%% pitch angle
h_FR_theta=figure('color','w');
x=[0 0 0]';
R=expmhat(pi/6*e2);
beta=15*pi/180;

[Q_R Q_L]=wing_attitude(beta,[0,pi/6,0]);
[h_body, h_wr, h_wl]=patch_monarch(fv_body, fv_wr, fv_wl, x, R, Q_R, Q_L, [0 0 0]);
alpha([h_body, h_wr, h_wl],0.9);
hold on;
view(360.54,2);
camlight;

axis auto;

x0=R*[0 mu_R(2) mu_R(3)]';
for i=1
    patch_arrow(x0, x0 + 100*R*Q_R*II(:,i), [1 0 0], lwidth, alength, awidth);
    patch_arrow(x0, x0 + 100*R*expmhat(beta*e2)*II(:,i), [0 1 0], lwidth, alength, awidth);
    plot_arc(x0, R*Q_R*II(:,i), R*expmhat(beta*e2)*II(:,i), 40);
end


%% devitation angle
h_FR_psi=figure('color','w');
x=[0 0 0]';
R=expmhat(pi/6*e2);
beta=15*pi/180;

[Q_R Q_L]=wing_attitude(beta,[0,0,pi/9]);
[h_body, h_wr, h_wl]=patch_monarch(fv_body, fv_wr, fv_wl, x, R, Q_R, Q_L, [0 0 0]);
alpha([h_body, h_wr, h_wl],0.9);
hold on;
view(270,35);
H = camlight(45, 80);
axis auto;

x0=R*[mu_R(1) mu_R(2) -20]';
for i=2
    patch_arrow(x0, x0 + 100*R*Q_R*II(:,i), [1 0 0], lwidth, alength, awidth);
    patch_arrow(x0, x0 + 100*R*expmhat(beta*e2)*II(:,i), [0 1 0], lwidth, alength, awidth);
    plot_arc(x0, R*Q_R*II(:,i), R*expmhat(beta*e2)*II(:,i), 40);
end

if bool_print
    figure(h_FB);
    print('monarch_FB','-depsc');
    figure(h_FR);
    print('monarch_FR','-depsc');
    figure(h_FL);
    print('monarch_FL','-depsc');
    figure(h_FS);
    print('monarch_FS','-depsc');
    figure(h_FR_phi);
    print('monarch_FR_phi','-depsc');
    figure(h_FR_theta);
    print('monarch_FR_theta','-depsc');
    figure(h_FR_psi);
    print('monarch_FR_psi','-depsc');

    !mv *.eps ../doc/Figs
end
end

function patch_circle(x_orig, x_normal, rad)
% generate a circle with
%   x_orig : the origin
%   x_normal : the vector normal to the circle
%   rad : radius

if size(x_normal,2)>1
    x_normal=x_normal(:);
end
if size(x_orig,2)>1
    x_org=x_orig(:);
end

N=501;
theta=linspace(0,2*pi,N);

bx=x_normal/norm(x_normal);
R=[bx null(bx')];
if det(R) < 0 
    R(:,3)=-R(:,3);
end

fv.vertices=zeros(N+1,3);
fv.faces=zeros(N,3);

fv.vertices(1,:)=x_orig';
for k=1:N
    fv.vertices(k+1,:) = (x_orig + rad*R*[0; cos(theta(k)); sin(theta(k))])';
end

fv.faces(:,1)=ones(N,1);
fv.faces(:,2)=(2:N+1)';
fv.faces(:,3)=[(3:N+1)'; 2];

face_color=[0.7 0.9 0.7];
h_circle=patch(fv);
set([h_circle], ...
    'FaceColor', face_color, ...
    'EdgeColor',       'none',        ...
    'FaceLighting',    'gouraud',     ...
    'AmbientStrength', 0.15);
alpha([h_circle],0.4);

end

function plot_arc(x0,u1,u2,rad)
% plot an arc with
%   x0 : origin
%   u1 : unit vector along the first direction
%   u2 : unit vector along the second direction
%   rad : radius of the arc
N=501;
theta=linspace(0,acos(dot(u1,u2)),N);
r1=u1;
r3=cross(u1,u2);
r3=r3/norm(r3);
r2=cross(r3,r1);
x_arc=zeros(3,N);
for k=1:N
    x_arc(:,k)=x0+rad*(r1*cos(theta(k))+r2*sin(theta(k)));
end
line(x_arc(1,:),x_arc(2,:),x_arc(3,:),'linewidth',1,'color',[0 0 0]);
end

function h=patch_arrow(p1,p2,acolor,lwidth,alength,awidth)
% acolor=[0 0 1];
% alength=1;
% awidth=1;
% lwidth=1;

if size(p1,1)>1
    p1=p1(:)';
end
if size(p2,1)>1
    p2=p2(:)';
end

dir=(p2-p1)/norm(p2-p1);
origin=p2-alength*dir;

tmp=rand(3,1);
base1=cross(dir,tmp)/norm(cross(dir,tmp));
base2=cross(dir,base1);

N=20;
theta=linspace(0,2*pi,N);
for a=1:N
    vertice(a,:)=origin+awidth*cos(theta(a))*base1+awidth*sin(theta(a))*base2;
end
vertice=[vertice;origin;origin+alength*dir];
faces=[1:N-1 1:N-1; 2:N 2:N; N+1*ones(1,N-1) N+2*ones(1,N-1)]';
ha=patch('Vertices',vertice,'Faces',faces,'FaceColor',acolor,'LineStyle','none');
hl=line([p1(1) origin(1)],[p1(2) origin(2)],[p1(3) origin(3)],'Color',acolor,'LineWidth',lwidth);

h=[ha; hl];
end

function update_arrow(h, p1,p2 ,alength,awidth)

if size(p1,1)>1
    p1=p1(:)';
end
if size(p2,1)>1
    p2=p2(:)';
end

dir=(p2-p1)/norm(p2-p1);
origin=p2-alength*dir;

tmp=rand(3,1);
base1=cross(dir,tmp)/norm(cross(dir,tmp));
base2=cross(dir,base1);

N=20;
theta=linspace(0,2*pi,N);
for a=1:N
    vertice(a,:)=origin+awidth*cos(theta(a))*base1+awidth*sin(theta(a))*base2;
end
vertice=[vertice;origin;origin+alength*dir];
faces=[1:N-1 1:N-1; 2:N 2:N; N+1*ones(1,N-1) N+2*ones(1,N-1)]';

set(h(1),'Vertices',vertice,'Faces',faces);
set(h(2),'XData',[p1(1) origin(1)], 'YData', [p1(2) origin(2)], 'ZData', [p1(3) origin(3)]);
end


