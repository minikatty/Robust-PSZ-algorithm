function [gamma, epsilon, sorted_vals] = get_bound_paras_all(mode, design_filename)
% 获取所有环境下的不确定性边界矩阵
% 输出:
%   gamma, epsilon: 结构体，包含 .B 和 .D 字段，尺寸均为[频点数 x 文件数]
%   sorted_vals:    排好序的环境变量值（如 51 个温度点），方便画图时直接作为 X 轴
    switch mode
        case 'snr', keyword = 'SNR';
        case 'temperature', keyword = 'T';
        case 'position', keyword = 'Pos';
    end

    folder_path = fullfile('./data/SimulateRIR/', mode);
    files = dir(fullfile(folder_path, ['*_' keyword '-*.mat']));
    if isempty(files), error('未找到文件'); end

    % 1. 提取所有文件对应的环境变量数值
    extract_num = @(name) str2double(regexp(name, ['Data_' keyword '-(-?\d+(?:\.\d+)?)'], 'tokens', 'once'));
    vals_cell = arrayfun(@(x) extract_num(x.name), files, 'UniformOutput', false);
    vals = cell2mat(vals_cell); 
    
    % 保证输出矩阵的列顺序与 X 轴完美对应
    [sorted_vals, sort_idx] = sort(vals);
    files = files(sort_idx);
    file_num = length(files);

    current_data = load(design_filename);
    num_freqs = size(current_data.ATF_BZ.ctrl, 3);
    
    RB_nominal = pagemtimes(current_data.ATF_BZ.ctrl, 'ctranspose', current_data.ATF_BZ.ctrl, 'none');
    RD_nominal = pagemtimes(current_data.ATF_DZ.ctrl, 'ctranspose', current_data.ATF_DZ.ctrl, 'none');

    % allocate
    epsilon.B = zeros(num_freqs, file_num);
    epsilon.D = zeros(num_freqs, file_num);
    gamma.B   = zeros(num_freqs, file_num);
    gamma.D   = zeros(num_freqs, file_num);

    for i = 1:file_num
        bound_file = fullfile(files(i).folder, files(i).name);
        bound_data = load(bound_file);

        deltaHB = current_data.ATF_BZ.ctrl - bound_data.ATF_BZ.ctrl; 
        deltaHD = current_data.ATF_DZ.ctrl - bound_data.ATF_DZ.ctrl;  

        deltaRB = RB_nominal - pagemtimes(bound_data.ATF_BZ.ctrl, 'ctranspose', bound_data.ATF_BZ.ctrl, 'none');
        deltaRD = RD_nominal - pagemtimes(bound_data.ATF_DZ.ctrl, 'ctranspose', bound_data.ATF_DZ.ctrl, 'none');

        epsilon.B(:, i) = squeeze(pagenorm(deltaHB, 'fro'));    
        epsilon.D(:, i) = squeeze(pagenorm(deltaHD, 'fro'));
        gamma.B(:, i)   = squeeze(pagenorm(deltaRB, 'fro'));
        gamma.D(:, i)   = squeeze(pagenorm(deltaRD, 'fro'));
    end
end