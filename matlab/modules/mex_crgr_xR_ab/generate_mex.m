% coder.extrinsic(@fcn);
% coder.varsize('X', [5001 length(X0)], [1 0] );

load('../../sim_QS_xR_hover_control_opt_mc.mat', 'INSECT', 'WK_R', 'WK_L', 't', 'des');
X0 = des.X0;
cfg = coder.config('mex');
% cfg.EnableAutoParallelization = true; % maybe not
cfg.InlineBetweenUserAndMathWorksFunctions = 'Always';
% cfg.DynamicMemoryAllocation = 'Off';
% codegen crgr_xR -args {INSECT, WK_R, WK_L, t, X0} -profile -launchreport -config cfg
N_max = 5001; % maximum number of integration points
codegen crgr_xR -args {INSECT, WK_R, WK_L, coder.typeof(t,[1,N_max],[0,1]), X0} -config cfg -launchreport

!cp crgr_xR_mex* ../