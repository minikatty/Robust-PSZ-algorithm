function export_public_listening_demos(renderer_file, output_dir, speech_file)
%EXPORT_PUBLIC_LISTENING_DEMOS Export the public 10-algorithm demo set.
%   The package contains three program types (speech, synthetic music, and
%   a 1-kHz tone).  For each item, the desired BZ reference, the desired
%   silent DZ reference, and the BZ/DZ outputs of ten algorithms are
%   exported at the horizontal 18-cm proxy pairs.
%
%   The synthetic music is generated in this file and can be redistributed.
%   No EBU SQAM audio or derivative is used by this exporter.

project_root = fileparts(fileparts(mfilename('fullpath')));
if nargin < 1 || isempty(renderer_file)
    renderer_file = fullfile(project_root, 'RQ_response', ...
        'private_inputs', 'broadband_renderers.mat');
end
if nargin < 2 || isempty(output_dir)
    output_dir = fullfile(project_root, 'RQ_response', 'listening_demo');
end
if nargin < 3 || isempty(speech_file)
    speech_file = fullfile(project_root, 'RQ_response', ...
        'private_inputs', '61-70968-0035.flac');
end
if ~exist(renderer_file, 'file')
    error('PublicDemo:MissingRenderer', ...
        'Renderer file does not exist: %s', renderer_file);
end
if ~exist(output_dir, 'dir')
    mkdir(output_dir);
end

loaded = load(renderer_file, 'renderers', 'reference', ...
    'point_definition');
renderers = loaded.renderers;
reference = loaded.reference;
point_definition = loaded.point_definition;
fs = reference.fs;

horizontal_indices = point_definition.horizontal_local_indices;
bz_indices = point_definition.bz_indices(horizontal_indices);
dz_indices = point_definition.dz_indices(horizontal_indices);
if numel(bz_indices) ~= 2 || numel(dz_indices) ~= 2
    error('PublicDemo:InvalidPointDefinition', ...
        'Expected two horizontal proxy points in each zone.');
end

algorithm_fields = { ...
    'ACC', 'ACC_Reg', 'PM', 'ACC_PM', 'wcRACC', ...
    'NoCT_WCRACC', 'Full_WCRACC', 'POTDC_RACC', ...
    'RPM', 'RACC_PM_Subpro'};
algorithm_labels = { ...
    'ACC', 'ACC-Reg', 'PM', 'ACC-PM', 'WCRACC', ...
    'NoCT-WCRACC', 'Full-WCRACC', 'POTDC-RACC', ...
    'WCRPM', 'RACC-PM'};
for i = 1:numel(algorithm_fields)
    if ~isfield(renderers, algorithm_fields{i})
        error('PublicDemo:MissingAlgorithm', ...
            'Renderer is missing algorithm field %s.', algorithm_fields{i});
    end
end

speech = prepare_speech(speech_file, fs);
music = synthesize_music(fs, 8.0);
tone = synthesize_tone(fs, 6.0, 1000);

items = struct( ...
    'id', {'S01', 'M00', 'T01'}, ...
    'category', {'speech', 'synthetic_music', 'tone'}, ...
    'source', {speech, music, tone});

rows = cell(0, 8);
for item_index = 1:numel(items)
    item = items(item_index);
    source = normalize_active_rms(item.source, fs);
    output_length = numel(source) + ...
        size(renderers.(algorithm_fields{1}).effective_irs, 2) - 1;
    ideal_bz = reference.target_pressure_pa .* render_fractional_delays( ...
        source, reference.total_delay_samples(bz_indices), output_length);
    ideal_dz = zeros(size(ideal_bz));

    rendered_bz = cell(numel(algorithm_fields), 1);
    rendered_dz = cell(numel(algorithm_fields), 1);
    global_peak = max(abs(ideal_bz(:)));
    for algorithm_index = 1:numel(algorithm_fields)
        effective_irs = renderers.(algorithm_fields{algorithm_index}) ...
            .effective_irs;
        rendered_bz{algorithm_index} = render_effective_irs( ...
            source, effective_irs(bz_indices, :));
        rendered_dz{algorithm_index} = render_effective_irs( ...
            source, effective_irs(dz_indices, :));
        global_peak = max(global_peak, ...
            max(abs(rendered_bz{algorithm_index}(:))));
        global_peak = max(global_peak, ...
            max(abs(rendered_dz{algorithm_index}(:))));
    end
    common_gain = min(1, 0.95 / max(global_peak, eps));

    reference_name = sprintf('%s_reference.wav', item.id);
    audiowrite(fullfile(output_dir, reference_name), ...
        ideal_bz .* common_gain, fs, 'BitsPerSample', 24);
    rows(end + 1, :) = {item.id, item.category, 'reference', ...
        'BZ', 'horizontal', reference_name, fs, common_gain}; %#ok<AGROW>

    dz_reference_name = sprintf('%s_DZ_reference.wav', item.id);
    audiowrite(fullfile(output_dir, dz_reference_name), ...
        ideal_dz, fs, 'BitsPerSample', 24);
    rows(end + 1, :) = {item.id, item.category, 'reference', ...
        'DZ', 'horizontal', dz_reference_name, fs, common_gain}; %#ok<AGROW>

    for algorithm_index = 1:numel(algorithm_fields)
        label = algorithm_labels{algorithm_index};
        safe_label = regexprep(label, '[^A-Za-z0-9-]', '-');
        filename = sprintf('%s_%s.wav', item.id, safe_label);
        audiowrite(fullfile(output_dir, filename), ...
            rendered_bz{algorithm_index} .* common_gain, fs, ...
            'BitsPerSample', 24);
        rows(end + 1, :) = {item.id, item.category, label, ...
            'BZ', 'horizontal', filename, fs, common_gain}; %#ok<AGROW>

        dz_filename = sprintf('%s_DZ_%s.wav', item.id, safe_label);
        audiowrite(fullfile(output_dir, dz_filename), ...
            rendered_dz{algorithm_index} .* common_gain, fs, ...
            'BitsPerSample', 24);
        rows(end + 1, :) = {item.id, item.category, label, ...
            'DZ', 'horizontal', dz_filename, fs, common_gain}; %#ok<AGROW>
    end
end

manifest = cell2table(rows, 'VariableNames', { ...
    'ItemID', 'Category', 'Algorithm', 'Zone', 'Orientation', ...
    'Filename', 'SampleRate_Hz', 'CommonDigitalGain'});
writetable(manifest, fullfile(output_dir, 'manifest.csv'));

metadata_file = fullfile(output_dir, 'generation_metadata.txt');
fid = fopen(metadata_file, 'w');
cleanup = onCleanup(@() fclose(fid));
fprintf(fid, 'Generated: %s\n', ...
    char(datetime('now', 'Format', 'yyyy-MM-dd HH:mm:ss')));
fprintf(fid, 'Renderer source file: %s\n', ...
    char(string(get_file_name(renderer_file))));
fprintf(fid, 'Sample rate: %d Hz\n', fs);
fprintf(fid, ['Channels: horizontal 18-cm BZ and DZ ', ...
    'pressure-proxy pairs\n']);
fprintf(fid, 'Algorithms: %s\n', strjoin(algorithm_labels, ', '));
fprintf(fid, ['Level handling: one common digital safety gain per ', ...
    'program item, shared by both zones, the references, and all ', ...
    'algorithms.\n']);
fprintf(fid, ['M00 is generated procedurally by ', ...
    'export_public_listening_demos.m; no SQAM material is used.\n']);
clear cleanup;

fprintf('Exported %d WAV files to %s.\n', height(manifest), output_dir);
end

function signal = prepare_speech(filename, target_fs)
if ~exist(filename, 'file')
    error('PublicDemo:MissingSpeech', ...
        'LibriSpeech source does not exist: %s', filename);
end
[signal, source_fs] = audioread(filename);
signal = mean(signal, 2);
if source_fs ~= target_fs
    [p, q] = rat(target_fs / source_fs, 1e-12);
    signal = resample(signal, p, q);
end
signal = double(signal(:));
signal = signal - mean(signal);
[b, a] = butter(6, 100 / (target_fs / 2), 'high');
signal = filtfilt(b, a, signal);
end

function signal = synthesize_music(fs, duration_s)
number_of_samples = round(duration_s * fs);
signal = zeros(number_of_samples, 1);
note_duration = 0.5;
melody_midi = [60, 64, 67, 72, 69, 67, 64, 62, ...
    60, 64, 67, 76, 72, 69, 67, 64];
chord_roots = [48, 53, 55, 48];
for note_index = 1:numel(melody_midi)
    first = round((note_index - 1) * note_duration * fs) + 1;
    last = min(number_of_samples, ...
        round(note_index * note_duration * fs));
    n = (0:(last - first)).';
    t = n / fs;
    frequency = 440 * 2^((melody_midi(note_index) - 69) / 12);
    envelope = note_envelope(numel(n), fs, 0.025, 0.12);
    note = (sin(2*pi*frequency*t) + ...
        0.32*sin(2*pi*2*frequency*t) + ...
        0.12*sin(2*pi*3*frequency*t)) .* envelope;
    signal(first:last) = signal(first:last) + 0.55 * note;
end
for chord_index = 1:numel(chord_roots)
    first = round((chord_index - 1) * 2 * fs) + 1;
    last = min(number_of_samples, round(chord_index * 2 * fs));
    n = (0:(last - first)).';
    t = n / fs;
    envelope = note_envelope(numel(n), fs, 0.08, 0.35);
    chord = zeros(size(t));
    for interval = [0, 4, 7]
        frequency = 440 * 2^((chord_roots(chord_index) + ...
            interval - 69) / 12);
        chord = chord + sin(2*pi*frequency*t) + ...
            0.18*sin(2*pi*2*frequency*t);
    end
    signal(first:last) = signal(first:last) + 0.12 * chord .* envelope;
end
signal = signal - mean(signal);
end

function signal = synthesize_tone(fs, duration_s, frequency_hz)
n = (0:(round(duration_s * fs) - 1)).';
signal = sin(2*pi*frequency_hz*n/fs);
fade_samples = round(0.05 * fs);
fade = 0.5 - 0.5*cos(pi*(0:(fade_samples - 1)).' / fade_samples);
signal(1:fade_samples) = signal(1:fade_samples) .* fade;
signal((end-fade_samples+1):end) = ...
    signal((end-fade_samples+1):end) .* flipud(fade);
end

function envelope = note_envelope(number_of_samples, fs, attack_s, release_s)
envelope = ones(number_of_samples, 1);
attack_samples = min(number_of_samples, max(1, round(attack_s * fs)));
release_samples = min(number_of_samples, max(1, round(release_s * fs)));
envelope(1:attack_samples) = linspace(0, 1, attack_samples).';
envelope((end-release_samples+1):end) = ...
    envelope((end-release_samples+1):end) .* ...
    linspace(1, 0, release_samples).';
end

function signal = normalize_active_rms(signal, fs)
signal = double(signal(:));
frame_length = round(0.02 * fs);
hop_length = round(0.01 * fs);
if numel(signal) < frame_length
    active_rms = sqrt(mean(signal.^2));
else
    number_of_frames = 1 + ...
        floor((numel(signal) - frame_length) / hop_length);
    frame_rms = zeros(number_of_frames, 1);
    for frame_index = 1:number_of_frames
        first = (frame_index - 1) * hop_length + 1;
        frame = signal(first:(first + frame_length - 1));
        frame_rms(frame_index) = sqrt(mean(frame.^2));
    end
    active = frame_rms >= max(frame_rms) * 10^(-40/20);
    active_rms = sqrt(mean(frame_rms(active).^2));
end
signal = signal ./ max(active_rms, eps);
end

function output = render_fractional_delays(source, delays, output_length)
source = double(source(:));
fft_length = 2^nextpow2(output_length + 1024);
source_spectrum = fft(source, fft_length);
signed_bins = [0:(fft_length/2), (-fft_length/2+1):-1].';
output = zeros(output_length, numel(delays));
for delay_index = 1:numel(delays)
    response = exp(-1i * 2*pi .* signed_bins .* ...
        delays(delay_index) / fft_length);
    response(fft_length/2 + 1) = cos(pi * delays(delay_index));
    current = real(ifft(source_spectrum .* response));
    output(:, delay_index) = current(1:output_length);
end
end

function output = render_effective_irs(source, effective_irs)
output_length = numel(source) + size(effective_irs, 2) - 1;
fft_length = 2^nextpow2(output_length);
source_spectrum = fft(source, fft_length);
response_spectrum = fft(effective_irs, fft_length, 2);
output = real(ifft(response_spectrum .* source_spectrum.', [], 2)).';
output = output(1:output_length, :);
end

function name = get_file_name(path_value)
[~, stem, extension] = fileparts(path_value);
name = [stem, extension];
end
