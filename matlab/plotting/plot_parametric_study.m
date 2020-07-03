%% File to plot the data
addpath('../modules', '../sim_data', '../');
set(0,'DefaultAxesFontName','times');
set(0,'DefaultAxesFontSize',18);
type = 'xR'; % types = {'x', 'xR', 'full'}

%%
if strcmp(type, 'x')
    load('parametric_study.mat');
    h_f_a_long = figure;
    h_f_a_long.PaperUnits = 'inches';
    h_f_a_long.PaperPosition = [0 0 13 9];
    nr = 2;
    nc = 3;
    %
    ic = 1;

    ir = 1;
    subplot(nr,nc,ic);
    plot(eps, squeeze(f_a_m(ic,ir,:,1)));
    ylabel('$\bar f_{a_1}$','interpreter','latex');
    hold on;
    yyaxis right;
    plot(eps, squeeze(f_a_m(ic,ir,:,2)), 'r');
    title(sprintf('$m_p$ = %0.3e \n $m_n$ = %0.3e\n', params.df_a_1_by_dphi_m), 'interpreter','latex');

    ir = 3;
    subplot(nr,nc,ic+nc);
    plot(eps, squeeze(f_a_m(ic,ir,:,2)), 'r');
    xlabel(sprintf('$\\Delta\\phi_{m_s}$'),'interpreter','latex');
    ylabel('$\bar f_{a_3}$','interpreter','latex');
    title(sprintf('$m_n$ = %0.3e\n', params.df_a_3_by_dphi_m),'interpreter','latex');
    %
    ic = 2;

    ir = 1;
    subplot(nr,nc,ic);
    plot(eps, squeeze(f_a_m(ic,ir,:,1)));
    hold on;
    yyaxis right;
    plot(eps, squeeze(f_a_m(ic,ir,:,2)), 'r');
    title(sprintf('$m_p$ = %0.3e \n $m_n$ = %0.3e\n', params.df_a_1_by_dtheta_0),'interpreter','latex');

    ir = 3;
    subplot(nr,nc,ic+nc);
    plot(eps, squeeze(f_a_m(ic,ir,:,2)), 'r');
    xlabel(sprintf('$\\Delta\\theta_0$'),'interpreter','latex');
    title(sprintf('$m_n$ = %0.3e\n', params.df_a_3_by_dtheta_0),'interpreter','latex');
    %
    ic = 4;

    ir = 1;
    subplot(nr,nc,ic-1);
    plot(eps, squeeze(f_a_m(ic,ir,:,1)));
    hold on;
    yyaxis right;
    plot(eps, squeeze(f_a_m(ic,ir,:,2)), 'r');
    title(sprintf('$m_p$ = %0.3e \n $m_n$ = %0.3e\n', params.df_a_1_by_dtheta_A_m),'interpreter','latex');

    ir = 3;
    subplot(nr,nc,ic+nc-1);
    plot(eps, squeeze(f_a_m(ic,ir,:,1)));
    hold on;
    yyaxis right;
    plot(eps, squeeze(f_a_m(ic,ir,:,2)), 'r');
    xlabel('$\Delta\theta_{A_m}$','interpreter','latex');
    title(sprintf('$m_p$ = %0.3e \n $m_n$ = %0.3e\n', params.df_a_3_by_dtheta_A_m),'interpreter','latex');
    %
    print(h_f_a_long, 'hover_param_study_long', '-depsc', '-r0');

    set(0,'DefaultAxesFontSize',24);
    h_f_a_lat = figure;
    plot(eps(eps<0), squeeze(f_a_m(3,2,eps<0,1)));
    ylabel('$\bar f_{a_2}$','interpreter','latex');
    hold on;
    plot(eps(eps>0), squeeze(f_a_m(3,2,eps>0,2)), 'r');
    xlabel(sprintf('$\\Delta\\phi_{m_k}$'),'interpreter','latex');
    title(sprintf('$m$ = %0.3e\n', params.df_a_2_by_dphi_m(2)),'interpreter','latex');
    print(h_f_a_lat, 'hover_param_study_lat', '-depsc', '-r0');
end

%%
if strcmp(type, 'xR')
    load('parametric_study_xR.mat');
    xlabels = {'$\Delta\phi_{m_s}$', '$\Delta\theta_{0_s}$', '$\Delta\phi_{m_k}$',...
        '$\Delta\phi_{0_s}$', '$\Delta\theta_{0_k}$', '$\Delta\psi_{0_k}$'};
    ylabels_f = {'$\bar f_{a_1}$', '$\bar f_{a_2}$' '$\bar f_{a_3}$'};
    ylabels_M = {'$\bar M_{a_1}$', '$\bar M_{a_2}$' '$\bar M_{a_3}$'};
    nr = 3;
    nc = 6;
    %
    h_f_a = figure;
    h_f_a.PaperUnits = 'inches';
    h_f_a.PaperPosition = [0 0 16 7];
    for ir=1:3
        for ic=1:6
            subplot(nr,nc,(ir-1)*nc+ic);
            plot(eps, squeeze(f_a_m(ic,ir,:,1)));
            if ic==1
                ylabel(ylabels_f{ir},'interpreter','latex');
            end
%             hold on;
%             yyaxis right;
%             plot(eps, squeeze(f_a_m(ic,ir,:,2)), 'r');
            if ir==3
                xlabel(xlabels{ic},'interpreter','latex');
            end
        end
    end
    print(h_f_a, 'hover_param_study_xRf', '-depsc', '-r0');
    %
    h_f_a = figure;
    h_f_a.PaperUnits = 'inches';
    h_f_a.PaperPosition = [0 0 16 7];
    for ir=1:3
        for ic=1:6
            subplot(nr,nc,(ir-1)*nc+ic);
            plot(eps, squeeze(M_a_m(ic,ir,:,1)));
            if ic==1
                ylabel(ylabels_M{ir},'interpreter','latex');
            end
%             hold on;
%             yyaxis right;
%             plot(eps, squeeze(M_a_m(ic,ir,:,2)), 'r');
            if ir==3
                xlabel(xlabels{ic},'interpreter','latex');
            end
        end
    end
    print(h_f_a, 'hover_param_study_xRM', '-depsc', '-r0');

%     title(sprintf('$m_p$ = %0.3e \n $m_n$ = %0.3e\n', params.df_a_1_by_dphi_m), 'interpreter','latex');
end
