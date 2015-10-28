% combine multiple fingerprints into single file
clear all

% number of fingerprints to combine
n_fp = 4;

folder = 'Z:\jtobin\gunshots\fingerprintLib\f_domain\';
fp_loc(1,:) = strcat(folder, 'A_33_s1_2048_4096_48k.coe');
fp_loc(2,:) = strcat(folder, 'J_0._s1_2048_4096_48k.coe');
fp_loc(3,:) = strcat(folder, 'Q_21_s1_2048_4096_48k.coe');
fp_loc(4,:) = strcat(folder, 'R_27_s1_2048_4096_48k.coe');

% initialize array
fp_data(:, 1) = zeros;

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
    
