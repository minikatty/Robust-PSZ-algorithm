clear; clc; close all;
addpath(genpath('src'));
addpath(genpath('data'));

geometry = load(fullfile('data', 'arrayGeometry', 'array_layout.mat'));
pos.src = geometry.roomArray.s;
pos.pos_center = [2.6, 2.0, 1.6];
pos.mics_pos = geometry.roomArray.BZ_ctrl;
pos.virtual_src_idx = 13;
freq_params = configure_freq_parameters();
ATF_desired_plane = plane_wave_generator(pos, freq_params);
save(fullfile('data', 'ATF_desired_plane.mat'), 'ATF_desired_plane');
