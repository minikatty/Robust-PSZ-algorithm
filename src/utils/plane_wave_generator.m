function [ATF_desired_plane] = plane_wave_generator(pos, freq_params)

src_pos = pos.src;
pos_center =pos.pos_center;
mics_pos = pos.mics_pos;
virtual_src_idx = pos.virtual_src_idx;
pos_spk_13 = src_pos(virtual_src_idx, :);   % 13号喇叭的坐标 (起点)

%% 2. 计算传播方向向量 (Direction Vector)
vec_prop = pos_center - pos_spk_13; 
vec_prop_norm = vec_prop / norm(vec_prop); 

%% 3. 计算波矢量 (Wave Vector k)
f_target = freq_params.target_freqs;            % 目标频率
c = 343;                    % 声速
k_mag = (2 * pi .* f_target) / c;

% k 向量 = 波数 * 单位方向向量
vec_k = k_mag .* vec_prop_norm(:)'; 

%% 4. 生成平面波目标向量 (Desired Response)
% 公式：d = exp(-j * k_vec * r)
phase_delay = mics_pos * vec_k'; 
% 得到目标向量 (M x num_target_freqs)
ATF_desired_plane = exp(-1j * phase_delay); 
% ATF_desired_plane = compute_atf(d_plane_wave, freq_params);

end