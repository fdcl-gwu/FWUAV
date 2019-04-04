clear all;
close all;

x_normal=rand(3,1);
x_orig=rand(3,1);
rad=20;

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
    fv.vertices(k+1,:) = (x_orig + R*[0; cos(theta(k)); sin(theta(k))])';
end

fv.faces(:,1)=ones(N,1);
fv.faces(:,2)=(2:N+1)';
fv.faces(:,3)=[(3:N+1)'; 2];

face_color=[0.4 0.4 0.4];
h_circle=patch(fv);
set([h_circle], ...
    'FaceColor', face_color, ...
    'EdgeColor',       'none',        ...
    'FaceLighting',    'gouraud',     ...
    'AmbientStrength', 0.15);
alpha([h_circle],0.8);
view(3);
