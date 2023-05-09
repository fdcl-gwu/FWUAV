%% Simple NN Model with algorithms : Multilayer Perceptron (2 layer linear activation)
evalin('base','clear all');
close all;
addpath('./modules', './sim_data', './plotting');
filename = 'iterative_learning_coil_noisyinp_test';
dnnTrain = dnnTrainFuns;

% load('sim_QS_xR_hover_control_opt_bigxR_T2_N30000', 't', 'X_ref0', ...
%     'WK_R', 'WK_L', 'INSECT', 'N_single', 'N_iters', 'N_per_iter', ...
%     'get_args', 'param_type', 'm', 'N_dang', 'Weights', 'N_sims', ...
%     'X0_pert', 'opt_param', 'des', 'opt_complete', 'cost', 'X_ref');

% X0_pert_flat = squeeze(X0_pert(:, 1, :));
% opt_param_flat = reshape(permute(opt_param, [1, 3, 2]), 3*N_sims, 6*N_iters);
% err = zeros(3*N_sims, 12);
% for i=1:(3*N_sims)
%     err(i, :) = get_error(X0_pert_flat(i, :)', des.X0)';
% end
% inputs = (1 ./ Weights.PerturbVariables)' .* err'; targets = opt_param_flat';
% opt_complete = opt_complete & (cost(:, end) < cost(:, 1));
% % inputs = inputs(:, opt_complete); targets = targets(:, opt_complete);
% % idx = 1:1000;
% idx = 1:450;
% inputs = inputs(:, opt_complete([idx, N_sims+idx, 2*N_sims+idx]));
% targets = targets(:, opt_complete([idx, N_sims+idx, 2*N_sims+idx]));

% load('python_training_data.mat');%, 'input_coil', 'target_coil', 'input_zeros');
idx = 1:5000;
load('python_multiple_input.mat');
% idx = 1:1500;
inputs = input_coil(idx, :)'; targets = target_coil(idx, :)';
inputs_multiple = {squeeze(input_coil(idx, 1, :))'; squeeze(input_coil(idx, 2, :))'};

% N_zero = 1000; % 25, 100
N_data = size(inputs, 2);
N_zero = round(N_data/15);
% N_zero = 300;
N_zero = 900;

inputs(:, (N_data+1):(N_data+N_zero)) = zeros(24, N_zero);
for i=1:length(inputs_multiple)
    inputs_multiple{i}(:, (N_data+1):(N_data+N_zero)) = zeros(12, N_zero);
end
targets(:, (N_data+1):(N_data+N_zero)) = zeros(60, N_zero);
% N_zero = size(input_zeros', 2);
% inputs = [inputs, input_zeros'];
% targets = [targets, zeros(60, N_zero)];

N_features = size(inputs, 1)/2;
[N_outputs, N_data] = size(targets);
N_iters = 5;

% inputs_orig = inputs; targets_orig = targets;
rand_idx = randperm(N_data);
inputs = inputs(:, rand_idx); targets = targets(:, rand_idx);
% inputs_orig = inputs_multiple;
inputs_mean = {}; inputs_scale = {};
for i=1:length(inputs_multiple)
    inputs_multiple{i} = inputs_multiple{i}(:, rand_idx);
    [inputs_multiple{i}, inputs_mean{i}, inputs_scale{i}] = normalize(inputs_multiple{i}, 2, 'range', [-1, 1]);
end
[targets, targets_mean, targets_scale] = normalize(targets, 2, 'range', [-1, 1]);
large_inp_idx = find(vecnorm(inputs_multiple{1}, 2, 1) > 1);

load('python_data.mat', 'Mean', 'Cov');
rng shuffle;
% GMModel = gmdistribution(bgm_means, permute(bgm_covars, [2,3,1]), bgm_weights);
% inputs = inputs_orig + random(GMModel, size(inputs, 2))';
% inputs = inputs_orig + mvnrnd(Mean, Cov, size(inputs, 2))';

%% Deep Neural Network
load('iterative_learning_dnn_tuning.mat', 'BayesObject');
optVars = BayesObject.bestPoint;
optVars.InitialLearnRate = 1e-3;

drop_prob = 0.0; % 0.25;
layers = [featureInputLayer(12, 'Normalization', 'none')
    
    fullyConnectedLayer(36) % 'Name','fc1'
    leakyReluLayer
    dropoutLayer(drop_prob)
    
%     fullyConnectedLayer(96)
%     leakyReluLayer
%     dropoutLayer(drop_prob)

%     concatenationLayer(1,2)
%     fullyConnectedLayer(60)

%     regressionLayer
    ];
lGraph = layerGraph(layers);
% lGraph = connectLayers(lGraph, 'input', 'concat/in2');

fc_idx = [];
for i = 1:length(lGraph.Layers)
    if contains(lGraph.Layers(i).Name, 'fc')
        fc_idx = [fc_idx, i];
    end
end

fc_idx = fc_idx(end);
grad_idx = [5,6]; % from dlnet.Learnables
% grad_idx = [1,2,3,4];

% rand_idx = randperm(size(inputs, 2));
% train_split = 0.8;
% train_idx = rand_idx(1:round(train_split*length(rand_idx)));
% val_idx = rand_idx(round(train_split*length(rand_idx)):end);
% inputs_train = inputs(:, train_idx); targets_train = targets(:, train_idx);
% inputs_val = inputs(:, val_idx); targets_val = targets(:, val_idx);

% miniBatchSize = 6000;
% maxEpochs = 1000;
% validationFrequency = 1;
% options = trainingOptions('adam', ...
%     'InitialLearnRate',optVars.InitialLearnRate, ...
%     'MaxEpochs',maxEpochs, ...
%     'MiniBatchSize',miniBatchSize, ...
%     'L2Regularization',optVars.L2Regularization, ...
%     'Shuffle','every-epoch', ...
%     'Verbose',false, ...
%     'Plots','training-progress', ...
%     'ValidationFrequency',validationFrequency, ...
%     'ValidationPatience',5, ...
%     'OutputNetwork','best-validation-loss', ...
%     'ExecutionEnvironment','gpu');
% [XTrain, YTrain, XValidation, YValidation] = train_val_split(inputs', targets', 0.8);
% inputs = XTrain'; targets = YTrain';

%% Setup
% rng default;
% alpha = 0.5;
% alphas = linspace(0.1, 0.5, N_iters);
alphas = 0.75 * ones(1, N_iters);
N_inputs = N_features;
solver = @fmincon;
problem.solver = func2str(solver);
problem.options = optimoptions(solver);
problem.options.Algorithm = 'interior-point'; % interior-point, sqp
% problem.options.SubproblemAlgorithm = 'cg';
problem.options.Display = 'iter';
problem.options.UseParallel = true;
problem.options.MaxFunctionEvaluations = 1e6;
problem.options.MaxIterations = 10; % 10, 25
% problem.options.ConstraintTolerance = 1e-10;
% problem.options.StepTolerance = 1e-10;
problem.options.ObjectiveLimit = 0;

%% Algorithm
y_star = targets;
perf = zeros(N_iters, 1); perf_z = zeros(N_iters, 1);
perf_yz = zeros(N_iters, 1);
cons = zeros(N_iters, 1); cons_z = zeros(N_iters, 1);
cs = cell(N_iters, 1); es = cell(N_iters, 1); os = cell(N_iters, 1);

ttotal = tic;
% diary filename
% PRETRAINING
disp('Starting initial training');
tinit = tic;

%%
% Training options
numEpochs = 100; % 50, 150;
learningRate = 1e-2;
executionEnvironment = "cpu";
validationFrequency = 5;
% validationCheck = true;
regularizationCoeff = 1e-2;

miniBatchSize = 500; % numEpochs = 100
miniBatchSize = N_data; % numEpochs = 250 (more needed with a larger batch)
% mbq = minibatchqueue(dsTrain,...
%     MiniBatchSize=miniBatchSize,...
%     MiniBatchFormat={'BC','BC','BC'},...
%     PartialMiniBatch="discard", ...
%     OutputEnvironment="cpu");
[mbq, mbq_val, mbq_full] = dnnTrain.train_val_ds(inputs_multiple, targets, 1, miniBatchSize, "cpu");
validationCheck = false;
[dlnet, netParams] = dnnTrain.initialize_nn_param(layers);
[dlnet, netParams, loss] = dnnTrain.train_custom_network(dlnet, netParams, mbq, ...
    numEpochs, learningRate, executionEnvironment, ...
    mbq_val, validationFrequency, validationCheck,...
    regularizationCoeff);

% shuffle(mbq_full);
[inp1, inp2, out] = next(mbq_full);
Y = dnnTrain.predictnet(dlnet, netParams, inp1, inp2);
mse(Y, out)
regression(extractdata(out)', extractdata(Y)', 'one')

% shuffle(mbq_val);
% [inp1, inp2, out] = next(mbq_val);
% Y = predictnet(dlnet, netParams, inp1, inp2);
% mse(Y, out)
% regression(extractdata(out)', extractdata(Y)', 'one')
% histogram(extractdata(log(rmse(Y, out, 1))))

% zero_inp1 = (dlarray(zeros(12,1), "CB") - inputs_mean{1}) ./ inputs_scale{1};
% zero_inp2 = (dlarray(zeros(12,1), "CB") - inputs_mean{2}) ./ inputs_scale{2};
% zero_out = (dlarray(zeros(60,1), "CB") - targets_mean) ./ targets_scale;

[zero_inp1, zero_inp2, zero_out] = dnnTrain.norm_data(dlarray(zeros(12,1), "CB"), ...
    dlarray(zeros(12,1), "CB"), dlarray(zeros(60,1), "CB"), ...
    inputs_mean, inputs_scale, targets_mean, targets_scale);

zero_pred = dnnTrain.predictnet(dlnet, netParams, zero_inp1, zero_inp2);
mse(zero_pred, zero_out)

% net = trainNetwork(inputs', targets', lGraph, options);
time_init = toc(tinit);
% net_init = net;
disp('Finished initial training');
% reset(mbq_full);
% [train_inp1, train_inp2, train_out] = next(mbq_full);
y = dnnTrain.predictnet(dlnet, netParams, dlarray(inputs_multiple{1}, "CB"), dlarray(inputs_multiple{2}, "CB"));
% w_star = get_fc_weights(net, fc_idx);
% w_y = w_star;

for i=1:N_iters
    %% CONSTRAINTS / MASTER STEP
	disp('Starting adjustment training for iteration ' + string(i));
%     options.ValidationData = {XValidation,((1-alphas(i))*YValidation + alphas(i)*net.predict(XValidation))};
    [mbq, mbq_val, mbq_full] = dnnTrain.train_val_ds(inputs_multiple, extractdata((1-alphas(i))*y_star + alphas(i)*y), ...
        1, miniBatchSize, "cpu");
    [dlnet, netParams] = dnnTrain.initialize_nn_param(layers);
    [dlnet, netParams, loss] = dnnTrain.train_custom_network(dlnet, netParams, mbq, ...
        numEpochs, learningRate, executionEnvironment, ...
        mbq_val, validationFrequency, false,...
        regularizationCoeff);
%     net = trainNetwork(inputs', ((1-alphas(i))*y_star + alphas(i)*y)', lGraph, options);
	disp('Finished adjustment training for iteration ' + string(i));
%         w0 = (1-alpha)*w_star + alpha*w_y;

%%
    w0 = dnnTrain.get_dl_weights(dlnet, netParams, fc_idx);
    N_w = length(w0);
    del = 0.2;
    problem.x0 = w0;
    % problem.lb = w0 - del * abs(w0);
    % problem.ub = w0 + del * abs(w0);

	corr_mat = eye(60) / N_outputs;
	F = zeros(length(w0));
%     grad_log_pi = zeros(length(w0), N_outputs);

%     dlGraph = net.layerGraph;
%     dlGraph = removeLayers(dlGraph, 'regressionoutput');
%     dlnet = dlnetwork(dlGraph);

	parfor k = 1:length(large_inp_idx)
		j = large_inp_idx(k);
%         inp = dlarray(inputs_orig(:, j), "CB");
        inp1 = dlarray(inputs_multiple{1}(:, j), "CB"); inp2 = dlarray(inputs_multiple{2}(:, j), "CB");
        [prediction, gradval] = dlfeval(dnnTrain.jac_out_weights, dlnet, netParams, inp1, inp2, length(w0), grad_idx);
        grad_log_pi = - gradval * corr_mat * (extractdata(prediction) - targets(:, j));
		F = F + grad_log_pi * grad_log_pi';
	end
	F = double(F) / length(large_inp_idx);

	problem.objective = @(w) fisher_norm(w, w0, F);
    % problem.objective = @(w) weight_fun(w, w0);
%     zero_inp = dlarray(zeros(N_inputs, 1), "CB");
%     zero_out = dlarray(zeros(N_outputs, 1), "CB");
    problem.nonlcon = @(w) zero_cons_dnn(w, dlnet, netParams, fc_idx, ...
        zero_inp1, zero_inp2, zero_out, length(w0), grad_idx, dnnTrain);
    problem.options.HessianApproximation = 'lbfgs';
    problem.options.SpecifyConstraintGradient = true;
    problem.options.ConstraintTolerance = 1e-10;
    problem.options.MaxIterations = 1000; % 25, 50
	problem.options.SpecifyObjectiveGradient = true;

	disp('Starting optimization for iteration ' + string(i));
	topt = tic;
    [w, prediction, exitflag, output] = solver(problem);
	es{i} = exitflag; os{i} = output;
	time_opt = toc(topt)
	disp('Finished optimization for iteration ' + string(i));
    [dlnet, netParams] = dnnTrain.set_dl_weights(dlnet, netParams, fc_idx, w);
%     net_z = net;

%% TRAINING / LEARNING
%     z = net.predict(inputs')';
    z = dnnTrain.predictnet(dlnet, netParams, dlarray(inputs_multiple{1}, "CB"), dlarray(inputs_multiple{2}, "CB"));
%     cons_z(i) = norm(net.predict(zeros(1,N_features)));
    zero_pred = dnnTrain.predictnet(dlnet, netParams, zero_inp1, zero_inp2);
    cons_z(i) = mse(zero_pred, zero_out);
    perf_z(i) = mse(z, y_star);

	disp('Starting unconstrained training for iteration ' + string(i));
% 	% inputs = inputs_orig + mvnrnd(Mean, Cov, size(inputs, 2))';
%     net = trainNetwork(inputs', z', lGraph, options);
    [mbq, mbq_val, mbq_full] = dnnTrain.train_val_ds(inputs_multiple, extractdata(z), 1, miniBatchSize, "cpu");
    [dlnet, netParams] = dnnTrain.initialize_nn_param(layers);
    [dlnet, netParams, loss] = dnnTrain.train_custom_network(dlnet, netParams, mbq, ...
        numEpochs, learningRate, executionEnvironment, ...
        mbq_val, validationFrequency, false,...
        regularizationCoeff);
	disp('Finished unconstrained training for iteration ' + string(i));
%     y = net.predict(inputs')';
    y = dnnTrain.predictnet(dlnet, netParams, dlarray(inputs_multiple{1}, "CB"), dlarray(inputs_multiple{2}, "CB"));

%     cons(i) = norm(net.predict(zeros(1,N_features)));
    zero_pred = dnnTrain.predictnet(dlnet, netParams, zero_inp1, zero_inp2);
    cons(i) = mse(zero_pred, zero_out);
    perf(i) = mse(y, y_star);
    perf_yz(i) = mse(y, z);
end

time_taken = toc(ttotal);
disp('Number of weight params near their bounds are,');
N_bad_param = sum(abs(abs((w-w0)./w0) - del) <= del * 1e-2);
disp(N_bad_param);
% diary off;

%% Performance and results
% plotregression(targets, control_net(inputs));
% plotregression(targets(:, tr.testInd), control_net(inputs(:, tr.testInd)));

allvars = whos;
tosave = cellfun(@isempty, regexp({allvars.class}, '^matlab\.(ui|graphics)\.'));
save(filename, allvars(tosave).name)
evalin('base',['load ' filename]);

%%

function err = get_error(X, Xd)
err = zeros(12, 1);
R = reshape(X(4:12), 3, 3); Rd = reshape(Xd(4:12), 3, 3);
err(1:3) = X(1:3) - Xd(1:3);
err(4:6) = 0.5*vee(Rd'*R - R'*Rd);
% [err(4), err(6), err(5)] = dcm2angle(R'*Rd, 'xzy');
err(7:9) = X(13:15) - Xd(13:15);
err(10:12) = X(16:18) - (R'*Rd*Xd(16:18));
end

function [J, dJ] = fisher_norm(w, w0, F)
J = 0.5 * (w-w0)' * F * (w-w0);

if nargout > 1
    dJ = F * (w-w0);
end

end

%%
function [XTrain, YTrain, XValidation, YValidation] = train_val_split(X, Y, split_ratio)
    N_data = size(X, 1);
    rand_idx = randperm(N_data);
    train_idx = rand_idx(1:round(split_ratio*N_data));
    val_idx = rand_idx((round(split_ratio*N_data)+1):end);
    XTrain = X(train_idx, :); YTrain = Y(train_idx, :);
    XValidation = X(val_idx, :); YValidation = Y(val_idx, :);
end

% function [c, ceq] = zero_cons_dnn(w, net, fc_idx, zero_inp, zero_out, num_wts, grad_idx)
function [c, ceq, gradc, gradceq] = zero_cons_dnn(w, dlnet, netParams, fc_idx, ...
    zero_inp1, zero_inp2, zero_out, num_wts, grad_idx, dnnTrain)
c = [];

% net = set_fc_weights(net, fc_idx, w);
% ceq = double(norm(net.predict(zeros(1, N_inp))'));

[dlnet, netParams] = dnnTrain.set_dl_weights(dlnet, netParams, fc_idx, w);
% ceq = double(norm(extractdata(predict(dlnet, zero_inp))));

if nargout > 2
	gradc = [];
    [ceq, gradceq] = dlfeval(dnnTrain.jac_err_weights, dlnet, netParams, zero_inp1, zero_inp2, zero_out, num_wts, grad_idx);
    ceq = double(extractdata(ceq));
end

end

function get_bayesian_params(beta, alpha, loss, gradientsSubnet, gradientsParams, N_data)
%             [beta, alpha] = get_bayesian_params(beta, alpha, loss, gradientsSubnet, gradientsParams, out.size(2));
% Does not work with this simple method because gamma value will always be 1.
tnet.Learnables = gradientsSubnet;
gradLoss = get_dl_weights(tnet, gradientsParams, []);
N_param = length(gradLoss);
lambda_max = beta * gradLoss'*gradLoss + alpha;
tr_inv = (N_param - 1) * 1/alpha + 1 / lambda_max;
gamma = N_param - alpha * tr_inv;
% J = rand(N_data, N_param);
% S = svd(J);
% eig_vals = beta * S.^2 + alpha * ones(size(S));
% tr_inv = sum(1./eig_vals) + abs(diff(size(J))) * 1 / alpha
end
