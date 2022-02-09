% coder.extrinsic(@fcn);
% coder.varsize('X', [5001 length(X0)], [1 0] );

% load('../../sim_QS_xR_pitch.mat', 'INSECT', 'WK', 't', 'X0');
load('crgr_sample_input.mat', 'INSECT', 'WK', 't', 'X0');
cfg = coder.config('mex');
% cfg.EnableAutoParallelization = true; % maybe not
cfg.InlineBetweenUserAndMathWorksFunctions = 'Always';
% cfg.DynamicMemoryAllocation = 'Off';
N_max = 5001; % maximum number of integration points
WK = orderfields(WK);
codegen crgr_xR_pitch -args {INSECT, WK, WK, coder.typeof(t,[1,N_max],[0,1]), X0} -config cfg -launchreport
copyfile crgr_xR_pitch_mex* ../
