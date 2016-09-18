function write_result_to_file(model_name,mic_auroc,mac_auroc,result_file)

fid = fopen(result_file, 'a+');
fprintf(fid,'%s\t',model_name);

for i=1:length(mic_auroc)
    fprintf(fid,'%f\t',mic_auroc(i));
end

for i=1:length(mac_auroc)
    fprintf(fid,'%f\t',mac_auroc(i));
end

fprintf(fid,'\n');
fclose(fid);
end

