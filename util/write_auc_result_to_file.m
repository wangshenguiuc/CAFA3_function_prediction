function result = write_auc_result_to_file(file_name,avg_mic_auroc,avg_mac_auroc,std_mac_auroc,specie)
%WRITE_AUC_RESULT_TO_FILE Summary of this function goes here
%   Detailed explanation goes here
result = [avg_mic_auroc,std_mac_auroc,avg_mac_auroc];
dlmwrite(['../result/PSB/',specie,'/',file_name,'.txt'],result,'delimiter','\t');


end

