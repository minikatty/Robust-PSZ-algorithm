function [gamma,epsilon] = get_bound_paras(mode, design_filename)

    switch mode
        case 'snr', keyword = 'SNR';
        case 'temperature', keyword = 'T';
        case 'position', keyword = 'Pos';
    end

    files = dir(fullfile(['./data/SimulateRIR/' mode], ['*_' keyword '-*.mat']));
    if isempty(files), error('未找到文件'); end

    extract_num = @(name) str2double(regexp(name, ['Data_' keyword '-(-?\d+(?:\.\d+)?)'], 'tokens', 'once'));
    vals_cell = arrayfun(@(x) extract_num(x.name), files,'UniformOutput', false);
    vals = cell2mat(vals_cell); 
    current_val = extract_num(design_filename);
    [~, max_idx] = max(abs(vals - current_val));   
    bound_file = files(max_idx).name;
    bound_data = load(bound_file);
    current_data = load(design_filename);
    deltaHB = current_data.ATF_BZ.ctrl - bound_data.ATF_BZ.ctrl; 
    deltaHD = current_data.ATF_DZ.ctrl - bound_data.ATF_DZ.ctrl;  
    deltaRB = pagemtimes(current_data.ATF_BZ.ctrl,'ctranspose',current_data.ATF_BZ.ctrl,'none') - ...
        pagemtimes(bound_data.ATF_BZ.ctrl,'ctranspose',bound_data.ATF_BZ.ctrl,'none');
    deltaRD = pagemtimes(current_data.ATF_DZ.ctrl,'ctranspose',current_data.ATF_DZ.ctrl,'none') - ...
    pagemtimes(bound_data.ATF_DZ.ctrl,'ctranspose',bound_data.ATF_DZ.ctrl,'none');

    epsilon.B = squeeze(pagenorm(deltaHB, 'fro'));    
    epsilon.D = squeeze(pagenorm(deltaHD, 'fro'));
    gammaB = squeeze(pagenorm(deltaRB, 'fro'));
    gamma.B = gammaB;
    gammaD = squeeze(pagenorm(deltaRD, 'fro'));
    gamma.D = gammaD;
end