function fv_monarch
%% generate faces and vertices for Monarch
close all;
clear all;
e1=[1 0 0]';
e2=[0 1 0]';
e3=[0 0 1]';

%% import an STL 
fv = stlread('./Monarch_reduced_5.stl');
N=length(fv.faces);

% Ns=[1:1000:N N];
% for i=1:length(Ns)-1
%     patch_selected_faces(fv,Ns(i),Ns(i+1)-1);
% end

%% split faces for each part by trial and error
N1=4390;
N2=7890;
N3=11450;
N4=17148;
N5=22600;
N_abdomen=[11971, 13016];
%patch_selected_faces(fv,1:N1);
% patch_selected_faces(fv,N1+1:N2);
% patch_selected_faces(fv,N2+1:N3);
%patch_selected_faces(fv,N3+1:15000);
%patch_selected_faces(fv,11971:13016);
% patch_selected_faces(fv,N4+1:N5);
% patch_selected_faces(fv,N5+1:N);



N_body=[1:N1 N3+1:11970 13017:N4];
N_wing=[N1+1:N2 N4+1:N5];
N_abdomen=11971:13016;

% patch_selected_faces(fv,[N_body N_wing]);
% xlabel('x');
% ylabel('y');
% zlabel('z');

scale_factor=0.06/146;

%% extract vertices/faces for body
tmp=fv.faces(N_body,:);
vert_index=unique(tmp(:));
fv_body.vertices=fv.vertices(vert_index,:);
fv_body.faces=zeros(length(N_body),3);
for i=1:length(fv_body.faces)
    for j=1:3
        fv_body.faces(i,j)=find(vert_index==fv.faces(N_body(i),j));
    end
end
R=[-e2 -e1 -e3];
x_shift=60;
for i=1:length(fv_body.vertices)    
    fv_body.vertices(i,:)=scale_factor*((R*fv_body.vertices(i,:)')'+[20 0 0]);
end

figure;
patch_selected_faces(fv_body,1:length(fv_body.faces));
hold on;
xlabel('x');
ylabel('y');
zlabel('z');
set(gca,'Ydir','reverse');
grid on;



%% extract vertices/faces for abdomen
tmp=fv.faces(N_abdomen,:);
vert_index=unique(tmp(:));
fv_abdomen.vertices=fv.vertices(vert_index,:);
fv_abdomen.faces=zeros(length(N_abdomen),3);
for i=1:length(fv_abdomen.faces)
    for j=1:3
        fv_abdomen.faces(i,j)=find(vert_index==fv.faces(N_abdomen(i),j));
    end
end
R=[-e2 -e1 -e3];
x_shift=40;
for i=1:length(fv_abdomen.vertices)    
    fv_abdomen.vertices(i,:)=scale_factor*((R*fv_abdomen.vertices(i,:)')'+[20 0 0]);
end

fv_abdomen.x_shift=-0.01179;

patch_selected_faces(fv_abdomen,1:length(fv_abdomen.faces));
hold on;
xlabel('x');
ylabel('y');
zlabel('z');
set(gca,'Ydir','reverse');

%% extract vertices/faces for right wing
tmp=fv.faces(N_wing,:);
vert_index=unique(tmp(:));
fv_wr.vertices=fv.vertices(vert_index,:);
fv_wr.faces=zeros(length(N_wing),3);
for i=1:length(fv_wr.faces)
    for j=1:3
        fv_wr.faces(i,j)=find(vert_index==fv.faces(N_wing(i),j));
    end 
end
R=[-e2 -e1 -e3];
x_shift=20;
Q=expmso3(5*pi/180*e1);
for i=1:length(fv_wr.vertices)    
    fv_wr.vertices(i,:)=scale_factor*((Q*R*fv_wr.vertices(i,:)')'+[20 0 0]);
end

patch_selected_faces(fv_wr,1:length(fv_wr.faces));

%% create left wing
fv_wl=fv_wr;
O_rl=[e1 -e2 e3];
for i=1:length(fv_wl.vertices)    
    fv_wl.vertices(i,:)=(O_rl*fv_wl.vertices(i,:)')';
end

patch_selected_faces(fv_wl,1:length(fv_wl.faces));


view(-90,90);
camlight(0,45);

filename='fv_monarch';
save(filename);
evalin('base',['load ' filename]);
evalin('base',['!cp ' filename '.mat ..']);

end


function patch_selected_faces(fv,Ns)

fv_tmp=fv;
fv_tmp.faces=fv_tmp.faces(Ns,:);

patch(fv_tmp,'FaceColor',       [0.8 0.8 1.0], ...
         'EdgeColor',       'none',        ...
         'FaceLighting',    'gouraud',     ...
         'AmbientStrength', 0.15);
title(['Faces ' num2str(min(Ns)) ' to ' num2str(max(Ns))]);

% Add a camera light, and tone down the specular highlighting

material('dull');

% Fix the axes scaling, and set a nice view angle
axis('image');
view(3);
%camlight(200,-45);
%view([-135 35]);
set(gca,'Zdir','reverse','YDir','reverse');
end