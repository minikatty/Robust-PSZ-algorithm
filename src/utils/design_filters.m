function [filters_w, algorithm_result] = design_filters( ...
        algorithm_name, design_data, para)
% DESIGN_FILTERS - A wrapper to call different filter design algorithms.

    algorithm_result = [];

    % Extract the necessary ATFs for filter design (usually control points)
    ATF_BZ_ctrl = design_data.ATF_BZ.ctrl;
    ATF_DZ_ctrl = design_data.ATF_DZ.ctrl;
    ATF_desired = design_data.ATF_desired;

    % virtual_src_idx = 13;  % virtual source index
    % ATF_desired = squeeze(ATF_BZ_ctrl(:, virtual_src_idx, :)); % M*L
    % load('./data/ATF_desired_plane.mat');
    % ATF_desired = ATF_desired_plane;
    % algorithms_to_test = {'ACC','PM','ACC_PM', 'wcRACC','RPM',...
    % 'RACC_PM_GLS','RACC_PM_Subpro', 'POTDC_RACC'};
    switch algorithm_name
        case 'ACC'
            fprintf(['Designing filters with unregularized ACC ', ...
                '(beta=0; conditioning warnings suppressed locally)...\n']);
            filters_w = ACC_Unregularized(ATF_BZ_ctrl, ATF_DZ_ctrl);

        case 'ACC_Reg'
            fprintf('Designing filters with ACC-Reg (beta=1e-6)...\n');
            filters_w = ACC(ATF_BZ_ctrl, ATF_DZ_ctrl);

        case 'M_ACC'
            fprintf(['Designing filters with M-ACC ', ...
                '(bright-zone complex-pressure uniformity)...\n']);
            [filters_w, algorithm_result] = M_ACC( ...
                ATF_BZ_ctrl, ATF_DZ_ctrl, para);

        case 'PM'
            % --- PRESSURE MATCHING ALGORITHM GOES HERE ---
            fprintf('Designing filters with PM...\n');
            filters_w = PM(ATF_BZ_ctrl, ATF_DZ_ctrl, ATF_desired);

        case 'ACC_PM'

            fprintf('Designing filters with ACC_PM...\n');
            filters_w = ACC_PM(ATF_BZ_ctrl, ATF_DZ_ctrl, ATF_desired, para);

        case 'wcRACC'
            fprintf('Designing filters with worst case ACC...\n');  
            wopt = wcACC(ATF_BZ_ctrl, ATF_DZ_ctrl, para);
            algorithm_result = wopt;
            if isfield(wopt,'scale')
               filters_w = wopt.w;
            else
               filters_w = wopt;
            end

        case 'RACC_Matched'
            fprintf('Designing parameter-matched conventional RACC filters...\n');
            if isfield(para, 'fre_idces')
                wopt = ConventionalRACCMatched( ...
                    ATF_BZ_ctrl, ATF_DZ_ctrl, para, para.fre_idces);
            else
                wopt = ConventionalRACCMatched(ATF_BZ_ctrl, ATF_DZ_ctrl, para);
            end
            filters_w = wopt.w;

        case 'NoCT_WCRACC'
            fprintf('Designing NoCT-WCRACC filters...\n');
            if isfield(para, 'fre_idces')
                wopt = NoCT_WCRACC( ...
                    ATF_BZ_ctrl, ATF_DZ_ctrl, para, para.fre_idces);
            else
                wopt = NoCT_WCRACC(ATF_BZ_ctrl, ATF_DZ_ctrl, para);
            end
            algorithm_result = wopt;
            filters_w = wopt.w;

        case 'Full_WCRACC'
            fprintf('Designing Full-WCRACC filters...\n');
            if isfield(para, 'fre_idces')
                wopt = Full_WCRACC( ...
                    ATF_BZ_ctrl, ATF_DZ_ctrl, para, para.fre_idces);
            else
                wopt = Full_WCRACC(ATF_BZ_ctrl, ATF_DZ_ctrl, para);
            end
            algorithm_result = wopt;
            filters_w = wopt.w;

        case 'RPM'
            % This is a placeholder. Replace with your actual function call.
            fprintf('Designing filters with MyRPM...\n');
            wopt = RPM(ATF_BZ_ctrl, ATF_DZ_ctrl, ATF_desired, para);
            algorithm_result = wopt;
            if isfield(wopt,'scale')
               filters_w = wopt.w;
            else
               filters_w = wopt;
            end
            
         case 'RACC_PM_GLS'
             if isfield(para,'fre_idces')
                 % --- RACC_PM Global Large Scale ---
                fprintf('Designing filters with global large scale RACC_PM @ few frequency bins...\n');
                fre_idces = para.fre_idces;
                wopt = RACC_PM_GLS_V2(ATF_BZ_ctrl, ATF_DZ_ctrl, ATF_desired, para, fre_idces);
                algorithm_result = wopt;
                if isfield(wopt,'scale')
                   filters_w = wopt.w;
                else
                   filters_w = wopt;
                end
             else
                % --- RACC_PM Global Large Scale ---
                fprintf('Designing filters with global large scale RACC_PM...\n');
                wopt = RACC_PM_GLS(ATF_BZ_ctrl, ATF_DZ_ctrl, ATF_desired, para);
                algorithm_result = wopt;
                if isfield(wopt,'scale')
                   filters_w = wopt.w;
                else
                   filters_w = wopt;
                end
             end

         case 'RACC_PM_Subpro'
             if isfield(para,'fre_idces')
                fprintf(['Designing filters with the maintained ', ...
                    'RACC_PM_Sub implementation @ selected bins...\n']);
                fre_idces = para.fre_idces;
                wopt = RACC_PM_Sub(ATF_BZ_ctrl, ATF_DZ_ctrl, ...
                    ATF_desired, para, fre_idces);
                algorithm_result = wopt;
                if isfield(wopt,'scale')
                   filters_w = wopt.w;
                else
                   filters_w = wopt;
                end
             else
                % --- Decompose the large scale problem into sub-contraints ---
                fprintf('Designing filters with Subproblems RACC_PM ...\n');
                wopt = RACC_PM_Sub(ATF_BZ_ctrl, ATF_DZ_ctrl, ATF_desired, para);
                algorithm_result = wopt;
                if isfield(wopt,'scale')
                   filters_w = wopt.w;
                else
                   filters_w = wopt;
                end
             end
         case 'POTDC_RACC'
            % --- Decompose the large scale problem into sub-contraints ---
            fprintf('Designing filters with POTDC_RACC ...\n');
            wopt = POTDC_RACC(ATF_BZ_ctrl, ATF_DZ_ctrl, para);
            algorithm_result = wopt;
            if isfield(wopt,'scale')
               filters_w = wopt.w;
            else
               filters_w = wopt;
            end            
        otherwise
            error('Unknown algorithm name: %s', algorithm_name);
    end

end
