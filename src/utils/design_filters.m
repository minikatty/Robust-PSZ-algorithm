function filters_w = design_filters(algorithm_name, design_data, para)
% DESIGN_FILTERS - A wrapper to call different filter design algorithms.

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
            % --- ACC ALGORITHM GOES HERE ---
            fprintf('Designing filters with ACC...\n');
            filters_w = ACC(ATF_BZ_ctrl, ATF_DZ_ctrl);

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
            if isfield(wopt,'scale')
               filters_w = wopt.w;
            else
               filters_w = wopt;
            end

        case 'POTDC-RACC'

            fprintf('Designing filters with POTDC_RACC...\n');
            wopt = POTDC_RACC(ATF_BZ_ctrl, ATF_DZ_ctrl, para);
            if isfield(wopt,'scale')
               filters_w = wopt.w;
            else
               filters_w = wopt;
            end

        case 'RPM'
            % This is a placeholder. Replace with your actual function call.
            fprintf('Designing filters with MyRPM...\n');
            wopt = RPM(ATF_BZ_ctrl, ATF_DZ_ctrl, ATF_desired, para);
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
                if isfield(wopt,'scale')
                   filters_w = wopt.w;
                else
                   filters_w = wopt;
                end
             else
                % --- RACC_PM Global Large Scale ---
                fprintf('Designing filters with global large scale RACC_PM...\n');
                wopt = RACC_PM_GLS(ATF_BZ_ctrl, ATF_DZ_ctrl, ATF_desired, para);
                if isfield(wopt,'scale')
                   filters_w = wopt.w;
                else
                   filters_w = wopt;
                end
             end

         case 'RACC_PM_Subpro'
             if isfield(para,'fre_idces')
                 % --- RACC_PM Global Large Scale ---
                fprintf('Designing filters with global Subproblem RACC_PM @ few frequency bins...\n');
                fre_idces = para.fre_idces;
                wopt = RACC_PM_Sub_V2(ATF_BZ_ctrl, ATF_DZ_ctrl, ATF_desired, para, fre_idces);
                if isfield(wopt,'scale')
                   filters_w = wopt.w;
                else
                   filters_w = wopt;
                end
             else
                % --- Decompose the large scale problem into sub-contraints ---
                fprintf('Designing filters with Subproblems RACC_PM ...\n');
                wopt = RACC_PM_Sub(ATF_BZ_ctrl, ATF_DZ_ctrl, ATF_desired, para);
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
            if isfield(wopt,'scale')
               filters_w = wopt.w;
            else
               filters_w = wopt;
            end            
        otherwise
            error('Unknown algorithm name: %s', algorithm_name);
    end

end