function filename = get_data_filename(data_dir, mode, parameter)
    switch mode
        case 'snr',       filename = fullfile(data_dir, sprintf('Data_SNR-%.2f.mat', parameter));
        case 'temperature', filename = fullfile(data_dir, sprintf('Data_T-%.2f.mat', parameter));
        case 'position',    filename = fullfile(data_dir, sprintf('Data_Pos-%.4f.mat', parameter));
    end
end