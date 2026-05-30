function [gamma,epsilon,freq_params] = get_real_measurement_bound(folder)
% ref:       reference matrix
% file_num:  total number of files
% folder:    folder path, e.g. 'data/Cabin_Measurements'
    files = dir(fullfile(folder, 'ATF*'));
    file_num = numel(files);
    ref = load(fullfile(folder,'Nominal_ATF.mat'));
    freq_params = ref.freq_params;freq_params.file_num = file_num;
    ATF_BZ_ref = ref.ATF_BZ; ATF_DZ_ref = ref.ATF_DZ;
    max_norm_deltaHB = zeros(size(ATF_BZ_ref,3),1);
    max_norm_deltaHD = zeros(size(ATF_DZ_ref,3),1);
    max_norm_deltaRB = zeros(size(ATF_BZ_ref,3),1);
    max_norm_deltaRD = zeros(size(ATF_DZ_ref,3),1);

    for i = 1:file_num
        filename = fullfile(folder, sprintf('ATF_%d.mat', i));
        current_data = load(filename);
        deltaHB = current_data.ATF_BZ - ATF_BZ_ref; 
        deltaHD = current_data.ATF_DZ - ATF_DZ_ref;
        deltaRB = pagemtimes(current_data.ATF_BZ,'ctranspose',current_data.ATF_BZ,'none') - ...
                    pagemtimes(ATF_BZ_ref,'ctranspose',ATF_BZ_ref,'none');
        deltaRD = pagemtimes(current_data.ATF_DZ,'ctranspose',current_data.ATF_DZ,'none') - ...
                    pagemtimes(ATF_DZ_ref,'ctranspose',ATF_DZ_ref,'none');
        norm_deltaHB = squeeze(pagenorm(deltaHB, 'fro'));    
        norm_deltaHD = squeeze(pagenorm(deltaHD, 'fro'));
        norm_deltaRB = squeeze(pagenorm(deltaRB, 'fro'));  
        norm_deltaRD = squeeze(pagenorm(deltaRD, 'fro'));
        max_norm_deltaHB = max(max_norm_deltaHB, norm_deltaHB);
        max_norm_deltaHD = max(max_norm_deltaHD, norm_deltaHD);
        max_norm_deltaRB = max(max_norm_deltaRB, norm_deltaRB);
        max_norm_deltaRD = max(max_norm_deltaRD, norm_deltaRD);
    end    
    epsilon.B = max_norm_deltaHB;    
    epsilon.D = max_norm_deltaHD;
    gamma.B   = max_norm_deltaRB;
    gamma.D   = max_norm_deltaRD;
end