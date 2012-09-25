function subfnPrintModelResults(data,fid)

fprintf(fid,'Model\n')
fprintf(fid,'%8s\t%8s\t%8s\t%8s\t%8s\n',' ','coeff','se','t','p')

fieldN = fieldnames(data);
for i = 1:length(fieldN)
    if ~strcmp(fieldN{i}, 'Model') && ~strcmp(fieldN{i}, 'Outcome')
        fprintf(fid,'%8s\t',fieldN{i})
        Str = sprintf('data.%s.beta',fieldN{i});
        fprintf('%8.4f\t',eval(Str));
        Str = sprintf('data.%s.se',fieldN{i});
        fprintf('%8.4f\t',eval(Str));
        Str = sprintf('data.%s.t',fieldN{i});
        fprintf('%8.4f\t',eval(Str));
        Str = sprintf('data.%s.p',fieldN{i});
        fprintf('%8.4f\n',eval(Str));
    end
end
        
