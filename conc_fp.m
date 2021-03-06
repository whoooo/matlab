% combine multiple fingerprints into single file
clear all

% number of fingerprints to combine
n_fp = 20;

% has trouble with s2?
folder = 'Z:\jtobin\gunshots\fingerprintLib\f_domain\';
fp_loc(1,:) = strcat(folder, 'A_33_s1_2048_4096_48k.coe');
fp_loc(2,:) = strcat(folder, 'W_27_s1_2048_4096_48k.coe');
fp_loc(3,:) = strcat(folder, 'B_23_s1_2048_4096_48k.coe');
fp_loc(4,:) = strcat(folder, 'B_23_s1_2048_4096_48k.coe');
fp_loc(5,:) = strcat(folder, 'E_21_s1_2048_4096_48k.coe');
fp_loc(6,:) = strcat(folder, 'E_21_s1_2048_4096_48k.coe');
fp_loc(7,:) = strcat(folder, 'F_46_s1_2048_4096_48k.coe');
fp_loc(8,:) = strcat(folder, 'F_46_s1_2048_4096_48k.coe');
fp_loc(9,:) = strcat(folder, 'J_0._s1_2048_4096_48k.coe');
fp_loc(10,:) = strcat(folder, 'J_0._s1_2048_4096_48k.coe');
fp_loc(11,:) = strcat(folder, 'L_21_s1_2048_4096_48k.coe');
fp_loc(12,:) = strcat(folder, 'L_21_s1_2048_4096_48k.coe');
fp_loc(13,:) = strcat(folder, 'Q_21_s1_2048_4096_48k.coe');
fp_loc(14,:) = strcat(folder, 'Q_21_s1_2048_4096_48k.coe');
fp_loc(15,:) = strcat(folder, 'R_27_s1_2048_4096_48k.coe');
fp_loc(16,:) = strcat(folder, 'R_27_s1_2048_4096_48k.coe');
fp_loc(17,:) = strcat(folder, 'S_16_s1_2048_4096_48k.coe');
fp_loc(18,:) = strcat(folder, 'S_16_s1_2048_4096_48k.coe');
fp_loc(19,:) = strcat(folder, 'V_25_s1_2048_4096_48k.coe');
fp_loc(20,:) = strcat(folder, 'V_25_s1_2048_4096_48k.coe');
% initialize array
fp_data(:, 1) = zeros;
% shot_name(n_fp,:) = zeros;
% fp_temp(n_fp,:) = zeros;

% store shot name, open file, read in binary fingerprint data, concatenate with other fingerprints. repeat
for i = 1 : n_fp    
    temp_fp_loc = fp_loc(i,:);
    shot_name_ind = strfind(temp_fp_loc, '_');
    shot_name(i,:) = temp_fp_loc(shot_name_ind(2)-1 : shot_name_ind(2) +2);
    fileID = fopen(temp_fp_loc);  
    fp_temp(i,:) = fread(fileID, '*char');
    fp_data = strcat(fp_data, fp_temp(i, 64:end - 2));
%     fp_data(i,:) = fp_temp(i, 64:end);
    fclose(fileID);
end

% radix and initialization vector for vivado mem generator
header = fp_temp(1,1:64);

% combined fingerprint data to save to file
fp_data = strcat(header, fp_data);

% get n_samples, nfft, and sampling rate for filename
fp_desc_ind = strfind(temp_fp_loc, '_s');
fp_desc = temp_fp_loc(fp_desc_ind(1) + 1 : end);

% concatenate shot names (R_27, A_33, etc)
shot_names = zeros;
for i = 1 : n_fp
    shot_names = strcat(shot_names, shot_name(i,:), '_');
end

% file name to be saved
save_loc = 'Z:\jtobin\gunshots\fingerprintLib\f_domain\combinations\';
fp_name = strcat(save_loc, shot_names, fp_desc);

% write data to file 
fp_write = fopen(fp_name, 'wt');
fprintf(fp_write, '%s ', fp_data);
fprintf(fp_write, '%c', ';');
fclose(fp_write);
% fprintf(fp_write, '%s\n', 'memory_initialization_radix=2;', 'memory_initialization_vector= ');
    
