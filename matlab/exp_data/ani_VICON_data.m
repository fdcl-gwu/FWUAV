function ani_VICON_data
clear all;
close all;
global FV_sph;

load VICON_data;
load fit_VICON_data;

FV_sph = FV_sphere(0.003);
bool_moving_camera=true;
bool_video=true;

k=1;
if ~bool_moving_camera
    h_RW=patch_chain(T1(:,k), RW1(:,k), RW2(:,k), RW3(:,k), T2(:,k), [0 0 1]);
    h_LW=patch_chain(T1(:,k), LW1(:,k), LW2(:,k), LW3(:,k), T2(:,k), [0 1 0]);
    h_thorax=patch_chain(T1(:,k),T2(:,k),A1(:,k),A2(:,k), [1 0 0]);
    axis equal;
    axis([-0.1 0.7 -0.1 0.1 -0.2 0.1]);
else
    h_RW=patch_chain(T1(:,k)-T1(:,k), RW1(:,k)-T1(:,k), RW2(:,k)-T1(:,k), RW3(:,k)-T1(:,k), T2(:,k)-T1(:,k), [0 0 1]);
    h_LW=patch_chain(T1(:,k)-T1(:,k), LW1(:,k)-T1(:,k), LW2(:,k)-T1(:,k), LW3(:,k)-T1(:,k), T2(:,k)-T1(:,k), [0 1 0]);
    h_thorax=patch_chain(T1(:,k)-T1(:,k), T2(:,k)-T1(:,k), A1(:,k)-T1(:,k), A2(:,k)-T1(:,k), [1 0 0]);
    %h_QR=patch_frame(T1(:,k)-T1(:,k),0.05*R(:,:,k)*Q_R(:,:,k),[0 0 1]);
    %h_QL=patch_frame(T1(:,k)-T1(:,k),0.05*R(:,:,k)*Q_L(:,:,k),[0 1 0]);
    %h_R=patch_frame(T1(:,k)-T1(:,k),0.05*R(:,:,k),[1 0 0]);
    %h_QA=patch_frame(T1(:,k)-T1(:,k),0.05*R(:,:,k)*Q_A(:,:,k),[1 0 0]);
    axis equal;
    axis([-0.1 0.1 -0.1 0.1 -0.1 0.1]);
end
set(gca,'ZDir','reverse','YDir','reverse');
view(3);
camlight;
material dull;
grid on;

if bool_video
    vidObj = VideoWriter('ani_VICON_data.avi');
    open(vidObj);
end


if ~bool_moving_camera
    for k=floor(linspace(1,length(t),101))      
        update_chain(h_thorax,T1(:,k),T2(:,k),A1(:,k),A2(:,k));
        update_chain(h_RW, T1(:,k), RW1(:,k), RW2(:,k), RW3(:,k), T2(:,k));
        update_chain(h_LW, T1(:,k), LW1(:,k), LW2(:,k), LW3(:,k), T2(:,k));
        axis([-0.1 0.7 -0.1 0.1 -0.2 0.1]);
        drawnow;
        if bool_video
            writeVideo(vidObj,getframe(gcf));
        end
    end
else
    for k=floor(linspace(1,length(t),101))        
        update_chain(h_RW, T1(:,k)-T1(:,k), RW1(:,k)-T1(:,k), RW2(:,k)-T1(:,k), RW3(:,k)-T1(:,k), T2(:,k)-T1(:,k));
        update_chain(h_LW, T1(:,k)-T1(:,k), LW1(:,k)-T1(:,k), LW2(:,k)-T1(:,k), LW3(:,k)-T1(:,k), T2(:,k)-T1(:,k));
        update_chain(h_thorax, T1(:,k)-T1(:,k), T2(:,k)-T1(:,k), A1(:,k)-T1(:,k), A2(:,k)-T1(:,k));
        %update_frame(h_QR,T1(:,k)-T1(:,k),0.05*R(:,:,k)*Q_R(:,:,k));
        %update_frame(h_QL,T1(:,k)-T1(:,k),0.05*R(:,:,k)*Q_L(:,:,k));
        %update_frame(h_R,T1(:,k)-T1(:,k),0.05*R(:,:,k));
        %update_frame(h_QA,T1(:,k)-T1(:,k),0.05*R(:,:,k)*Q_A(:,:,k));
        axis(0.7*[-0.1 0.1 -0.1 0.1 -0.1 0.1]);
        view(20,20);
        drawnow;
        
        if bool_video
            writeVideo(vidObj,getframe(gcf));
        end
    end
end


if bool_video
    close(vidObj);
end

end

function FV = FV_sphere(rad)
% generate a patch object for a sphere with a given radius

[sph.x sph.y sph.z]=sphere(51);
sph.x=sph.x*rad;
sph.y=sph.y*rad;
sph.z=sph.z*rad;
FV = surf2patch(sph.x, sph.y, sph.z);
end

function h=patch_chain(varargin)
N=nargin-1;
color=varargin{N+1};

for i=1:N
    h.sph(i) = patch_sphere(varargin{i},color);
end

for i=1:N-1
    tmp=zeros(3,2);
    for j=1:3
        tmp(j,:)=[varargin{i}(j) varargin{i+1}(j)];
    end
    h.line(i) = line(tmp(1,:), tmp(2,:), tmp(3,:), 'color', color, 'LineWidth', 2);
end
end

function update_chain(h, varargin)
N=nargin-1;

for i=1:N
    update_sphere(h.sph(i),varargin{i});
end

for i=1:N-1
    tmp=zeros(3,2);
    for j=1:3
        tmp(j,:)=[varargin{i}(j) varargin{i+1}(j)];
    end
    set(h.line(i),'XData',tmp(1,:),'YData',tmp(2,:), 'ZData', tmp(3,:));
end
end

function h=patch_sphere(x,color)
global FV_sph;
if ~iscolumn(x)
    x=x';
end

N_vertices=length(FV_sph.vertices);
vertices=FV_sph.vertices + ones(N_vertices,1)*x';

h=patch('Faces',FV_sph.faces, 'Vertices', vertices,'FaceColor', color, 'LineStyle','none');
end

function update_sphere(h,x)
global FV_sph;
if ~iscolumn(x)
    x=x';
end

N_vertices=length(FV_sph.vertices);
vertices=FV_sph.vertices + ones(N_vertices,1)*x';

set(h,'Vertices', vertices);

end

function h = patch_frame(x,R,acolor)

lwidth=2;
alength=0.01;
awidth=alength*tand(12);

for i=1:3    
    h(:,i)=patch_arrow(x,x+R(:,i),acolor,lwidth,alength,awidth);
end

end

function update_frame(h,x,R)

lwidth=2;
alength=0.01;
awidth=alength*tand(12);

for i=1:3    
    update_arrow(h(:,i), x,x+R(:,i) ,alength,awidth);
end

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




