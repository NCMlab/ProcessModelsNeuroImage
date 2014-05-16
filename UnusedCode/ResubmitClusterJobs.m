OutPath = pwd



List = [1:130]
F = dir('Results*.mat');
for i = 1:length(F)
    name = F(i).name;
    number = str2num(name(9:12));
    List(number) = 0;
end
List = List(find(List));



for i = 1:length(List)
    
    fprintf(1,['qsub -e ' OutPath ' -o ' OutPath ' ' sprintf('job_%04d.sh\n',List(i))]);
end
