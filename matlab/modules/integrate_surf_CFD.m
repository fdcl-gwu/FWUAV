function [F_I,F_B,F_R,M_I,M_B,M_R] = integrate_surf_CFD(path,timeStep,Q_R,R,xState)
fid = fopen([path,'/output_aero_surf',sprintf('%05d',timeStep),'.dat'],'r');
datacell = textscan(fid,'%f%f%f%f%f%f%f','HeaderLines',1,'Collect',1);
fclose(fid);
surf.data = datacell{1};
BN = surf.data(:,1)';
xI = surf.data(:,2)';
yI = surf.data(:,3)';
zI = surf.data(:,4)';
FxI = surf.data(:,5)';
FyI = surf.data(:,6)';
FzI = surf.data(:,7)';

xCG = xState(1); yCG = xState(2); zCG = xState(3);

for i = 1:length(BN)
    FI_C(:,i) = [FxI(i); FyI(i); FzI(i)];
    FB_C(:,i) = R'*FI_C(:,i);
    FR_C(:,i) = Q_R'*R'*FI_C(:,i);
    
    xB_C(:,i) = R'*[xI(i)-xCG; yI(i)-yCG; zI(i)-zCG];
    xR_C(:,i) = Q_R'*R'*[xI(i)-xCG; yI(i)-yCG; zI(i)-zCG];
    
    MI_C(:,i) = [FzI(i)*yI(i)-FyI(i)*zI(i); ...
                 FxI(i)*zI(i)-FzI(i)*xI(i); ...
                 FyI(i)*xI(i)-FxI(i)*yI(i)];
    
    MB_C(:,i) = [FB_C(3,i)*xB_C(2,i)-FB_C(2,i)*xB_C(3,i); ...
                 FB_C(1,i)*xB_C(3,i)-FB_C(3,i)*xB_C(1,i); ...
                 FB_C(2,i)*xB_C(1,i)-FB_C(1,i)*xB_C(2,i)];
    
    MR_C(:,i) = [FR_C(3,i)*xR_C(2,i)-FR_C(2,i)*xR_C(3,i); ...
                 FR_C(1,i)*xR_C(3,i)-FR_C(3,i)*xR_C(1,i); ...
                 FR_C(2,i)*xR_C(1,i)-FR_C(1,i)*xR_C(2,i)];
    
end

F_I = sum([FxI;FyI;FzI],2);
F_B = sum(FB_C,2);
F_R = sum(FR_C,2);

M_I = sum(MI_C,2);
M_B = sum(MB_C,2);
M_R = sum(MR_C,2);

end

