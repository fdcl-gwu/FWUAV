function tmp_arrow
close all;

x_cp=zeros(3,1);
lwidth=1;
alength=1.2;
awidth=alength*tand(10);
acolor=[1 0 0]';
scale_Force=1e6;
[h_La, h_Ll]=patch_arrow(x_cp, x_cp +[ 0 0 0]', acolor, lwidth, alength, awidth);
pause;
update_arrow(h_La, h_Ll, x_cp, x_cp +[ 10 10 0]', alength, awidth);
axis(10*[-1 1 -1 1 -1 1]);
axis equal;

end


function [ha hl]=patch_arrow(p1,p2,acolor,lwidth,alength,awidth)
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

end

function update_arrow(ha, hl, p1,p2 ,alength,awidth)

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

set(ha,'Vertices',vertice,'Faces',faces);
set(hl,'XData',[p1(1) origin(1)], 'YData', [p1(2) origin(2)], 'ZData', [p1(3) origin(3)]);

end


