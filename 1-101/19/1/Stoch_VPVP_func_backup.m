function [] = Stoch_VPVP_func(results_folder,beginn, endd)
tic
if(isdeployed==false)
    addpath(genpath(pwd));
end

% results_folder = [pwd '/Results/New/calcium_multi1/2/'];
FileNameAndLocation=[mfilename('fullpath')];
currentfile=strcat(FileNameAndLocation, '.m');
newbackup = [results_folder mfilename '_backup.m'];
copyfile(currentfile,newbackup);


%% Initialization of parameters
randomness = 1;
plotting = 1;
load_previous_data = 0;

%% Stochastic Search parameters
iterations = 600000;%100000
n_alphas = 3;
random_seed = 3;
alphas_dispersion = 1;
err_weight = 0.1;
W_start = 45000; %1000
random_factor1 = 0.1;
random_factor2 = 0.1;
w_random_factor1 = 0.01;
w_random_factor2 = 0.1;
high_jump_period_W = 30;
high_jump_period_a = 31;

%% VPVP parameters
optim_freq = 3000;
lambda = 3e8;
downsample = 10;
MaxIter_X = 200;
MaxIter = 50;
optTol_X = 1e-5/downsample;
progTol_X = 1e-9/downsample;
optTol_aW = 1e-3/downsample;
progTol_aW = 1e-6/downsample;

%% Load or create data to fit

% origin_of_data = 'calcium';
origin_of_data = 'huntington';
% origin_of_data = 'simulation';

switch origin_of_data
    case 'calcium';
        load svd_masked;
%         beginn = 365;
%         endd = 465;
        ncomp = 6;
        X_true = 10*U(beginn:endd-1,2:ncomp+1);

    case 'huntington';
        load subj_000_samp_0;
        beginn = 1;
        endd = 101;
        components = [1 2 3 4 5 6 7 8 9 14 15 16 17 18 20];
        ncomp = 15;
        X_true = fmri_compnent_tc(beginn:endd-1,components);
        A = mean(std(X_true));
        X_true = X_true*(1/A);
        x_start = fmri_compnent_tc(beginn,components)'*(1/A);

    case 'simulation';
        load('workspace_diagW0_7.mat','alphas','W','x_start','y_start','param')
        alphas = struct2cell(alphas);
        [X_true, Y] = calcium_dynamic(x_start, y_start, alphas, W, param);
        X_true = X_true';
        ncomp = 6;
end

if randomness
    rng('shuffle')
else
    rng('default')
    rng(random_seed);
end

% Integration settings
x_start = X_true(1,:)';
y_start = 1*randn(ncomp,1);
W = zeros(ncomp);
Wp = W;
alphas{1} = 0.578*ones(ncomp,1);
alphas{2} = 1*ones(ncomp,1);
alphas{3} = 11.38*ones(ncomp,1);
N = size(X_true,1);
param.N = N;
param.dt = 0.1;
param.ncomp = ncomp;
param.n_alphas = n_alphas;


% intialization of some varaibles
mean_corr = [];
mean_error = [];
iter_improve = [];
max_mean_corr = 0;
ce_min = -inf;
iter_ini = 1;
if randomness
    optim_start = floor(rand*1000)
else
    optim_start = 1;
end

if load_previous_data
load('/Users/ger/Documents/LearningSmoother/Results Oficina/calcium_best_fittings/params_VARIOUS_Stoch_VPVP_Search.mat')
if ~isa(alphas,'cell')
    alphas = struct2cell(alphas);
end
    if exist('y_start_new','var') == 1
        y_start = y_start_new;
        W = Wans;
    end
%     param.n_alphas = param.amount_of_alphas;
    [X, Y] = calcium_dynamic(x_start, y_start, alphas, W, param);
    c = measure_corr(X',X_true);
    error = SMAPE(X', X_true);
    cmin = min(c);
    corr_error = c - err_weight.*error;
    ce_min = min(corr_error);
    figure(1);clf;
    plott(X_true', X, beginn, endd, W, c, ncomp)
    W_start = 0;
    iterations = 200000;
        if exist('iter','var') == 1
            iter_ini = iter+1;
        end
%     optim_start = 5000;%5000
end
improvement_count = size(mean_error,1);
optim_steps = [];
%% stoch optimization
for iter=iter_ini:iterations
    if mod(iter,100) == 0,
        disp(['>>>>> iter = ', num2str(iter)]);
    end

    if mod(iter,high_jump_period_W) == 0,%600
        w_random_factor = w_random_factor2;
    else
        w_random_factor = w_random_factor1;
    end

    % start changing W only later in the process, when the indivudual time
    % series are (hopefully) modeleled well enough by now; alternatively,
    % maybe another condition here - like current model fit' - e.g. minimal
    % correlation?
    if iter > W_start
        stoch_step_W = w_random_factor*randn(size(W));
        stoch_step_W(eye(size(stoch_step_W))~=0)=0;
        Wp = W + stoch_step_W;
    end;
    if mod(iter,high_jump_period_a) == 0,%500
        random_factor = random_factor2;
        full_stochastic_step
    else
        random_factor = random_factor1;
        stochastic_step
    end

    %%  generate time series (X) given the current parameters
    [X, Y] = calcium_dynamic(x_start, y_startp, alphasp, Wp, param);
    if size(X,2) == N % check if integration was successfully completed
        % measure how close is simulated data respect to the true data
        c = measure_corr(X',X_true);
        error = SMAPE(X', X_true);
%         c_last = measure_corr(X(:,1:end-8)',X_true(1:end-8,:));
%         corr_error = c_last/2 + c - err_weight.*error; % NORMALIZAR EL ERROR
        corr_error = c - err_weight.*error; % NORMALIZAR EL ERROR
        if min(corr_error) > ce_min & ~sum(sum(isnan(X)))
            improvement_count = improvement_count + 1
            ce_min = min(corr_error);
            cmin = min(c);
            disp(['Stochastic improvement >>> min_corr = ',num2str(cmin)]);

            y_start = y_startp;
            alphas = alphasp;
            W = Wp;
            mean_c = mean(c);
            mean_corr = [mean_corr; mean_c];
            mean_error = [mean_error; mean(error)];
            iter_improve = [iter_improve; iter];

            if mean_c > max_mean_corr
                best_alphas = alphas;
                best_W = W;
                best_y_start = y_start;
                best_c = c;
                max_mean_corr = mean_c;
                best_iter = iter;
            end

            if plotting == 3
                figure(1);clf;
                plott(X_true', X, beginn, endd, Wp, c, ncomp)
                figure(2)
                plot(mean_corr,'-.O'); grid on; hold on;
                plot(mean_error,'-.O'); grid on;
                legend('mean correlation', 'mean SMAPE error','Location','NorthWest')
                drawnow; hold off;
            end
        end

        %% Optimization
%         if 0
        if or(mod(iter,optim_freq)==0, iter == optim_start),
            disp('>>>>> VPVP')
            VPVP_params.lambda = lambda;
            VPVP_params.downsample = downsample;
            VPVP_params.MaxIter_X = MaxIter_X;
            VPVP_params.MaxIter = MaxIter;
            VPVP_params.optTol_X = optTol_X;
            VPVP_params.progTol_X = progTol_X;
            VPVP_params.optTol_aW = optTol_aW;
            VPVP_params.progTol_aW = progTol_aW;
            if or(mod(iter,100*optim_freq)==0, iter == optim_start)
                save([results_folder 'VPVP_calcium_iter_preOptim' num2str(iter)],'y_start','alphas','W','iter','c','mean_corr','mean_error','iter_improve','optim_start','optim_steps','VPVP_params')
            end
            [y_startp, alphasp, Wp] = call_VPVP(X_true, y_start, alphas, W, VPVP_params, param);

            [X, Y] = calcium_dynamic(x_start, y_startp, alphasp, Wp, param);
            if size(X,2) == N
                disp('VPVP accepted')
                improvement_count = improvement_count + 1
                optim_steps = [optim_steps improvement_count]
                c = measure_corr(X',X_true);
                error = SMAPE(X', X_true);
%                 c_last = measure_corr(X(:,1:end-8)',X_true(1:end-8,:));
                corr_error = c - err_weight.*error;
%                 corr_error = c_last/2 + c - err_weight.*error; % NORMALIZAR EL ERROR
                y_start = y_startp;
                alphas = alphasp;
                W = Wp;
                if min(c) > cmin & ~sum(sum(isnan(X)))
                    disp(['VPVP improvement >>> corr_min = ',num2str(min(c))]);
                end
                if min(corr_error) > ce_min & ~sum(sum(isnan(X)))
                    disp(['VPVP improvement >>> ce_min = ',num2str(min(corr_error))]);
                end
                % change ce_min even if there was no improvement
                cmin = min(c);
                ce_min = min(corr_error);
                mean_corr = [mean_corr; mean(c)];
                mean_error = [mean_error; mean(error)];
                iter_improve = [iter_improve; iter];

                if mean_c > max_mean_corr
                    best_alphas = alphas;
                    best_W = W;
                    best_y_start = y_start;
                    best_c = c;
                    max_mean_corr = mean_c;
                    best_iter = iter;
                end

                if plotting >= 2
                    figure(1);clf;
                    plott(X_true', X, beginn, endd, W, c, ncomp)
                    figure(2)
                    plot(mean_corr,'-.O'); grid on; hold on;
                    plot(mean_error,'-.O'); grid on;
                    legend('mean correlation', 'mean SMAPE error','Location','NorthWest')
                    drawnow; hold off;
                end
                if or(mod(iter,100*optim_freq)==0, iter == optim_start)
                    save([results_folder 'VPVP_calcium_iter_postOptim' num2str(iter)],'y_start','alphas','W','iter','c','mean_corr','mean_error','iter_improve','optim_start','optim_steps','VPVP_params')
                end
            end
        end

%         if mod(iter,1000)==0,
%             disp(['>>>>> iter ',num2str(iter)]);
%         end
    end
end

if plotting >= 1
    [X, Y] = calcium_dynamic(x_start, best_y_start, best_alphas, best_W, param);
    c = measure_corr(X',X_true);
    figure(3)
    plott(X_true', X, beginn, endd, W, c, ncomp)
    title(['Best Fit - iter ' num2str(best_iter)])
    figure(2)
    plot(mean_corr,'-.O'); grid on; hold on;
    plot(mean_error,'-.O'); grid on;
    legend('mean correlation', 'mean SMAPE error','Location','NorthWest')
    drawnow; hold off;
end

best_c
max_mean_corr
save([results_folder 'best_params'],'best_alphas','best_W','best_y_start','best_c','max_mean_corr','best_iter')
toc
end
