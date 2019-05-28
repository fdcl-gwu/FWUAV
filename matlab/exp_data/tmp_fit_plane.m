clear all;
close all;

load VICON_data.mat

k=1;
M_RR = (RW1(:,k)-T1(:,k))*(RW1(:,k)-T1(:,k))'...
    + (RW2(:,k)-T1(:,k))*(RW2(:,k)-T1(:,k))'...
    + (RW3(:,k)-T1(:,k))*(RW3(:,k)-T1(:,k))';

tmp=cross( RW1(:,k)-T1(:,k), RW2(:,k)-T1(:,k) );
tmp/norm(tmp)

tmp=cross( RW1(:,k)-T1(:,k), RW3(:,k)-T1(:,k) );
tmp/norm(tmp)

tmp=cross( RW2(:,k)-T1(:,k), RW3(:,k)-T1(:,k) );
tmp/norm(tmp)

[V D]=eig(M_RR)

n=V(:,1)
abs(n'*(RW1(:,k)-T1(:,k))) + abs(n'*(RW2(:,k)-T1(:,k))) + abs(n'*(RW3(:,k)-T1(:,k)))

n=V(:,2)
abs(n'*(RW1(:,k)-T1(:,k))) + abs(n'*(RW2(:,k)-T1(:,k))) + abs(n'*(RW3(:,k)-T1(:,k)))

n=V(:,3)
abs(n'*(RW1(:,k)-T1(:,k))) + abs(n'*(RW2(:,k)-T1(:,k))) + abs(n'*(RW3(:,k)-T1(:,k)))

