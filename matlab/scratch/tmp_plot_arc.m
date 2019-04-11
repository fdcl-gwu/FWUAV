clear all;
close all;

x0=[1 1 1]';
u1=[1 0 0]';
u2=[0 1 1]';
rad=1;

u1=u1/norm(u1);
u2=u2/norm(u2);

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
line(x_arc(1,:),x_arc(2,:),x_arc(3,:),'linewidth',2);
figure;
plot3(x_arc(1,:),x_arc(2,:),x_arc(3,:));




