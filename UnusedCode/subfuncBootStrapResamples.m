function Resamples = subfuncBootStrapResamples(Nboot,NSub,Stratification)

Resamples = uint16(zeros(NSub,Nboot));

if ~isempty(Stratification)
    Gr1 = find(Stratification == 0);
    NGr1 = length(Gr1);
    Gr2 = find(Stratification == 1);
    NGr2 = length(Gr2);
else
    NGr1 = [];
    NGr2 = [];
end

for i = 1:Nboot
    if isempty(Stratification)
        Samp =  floor(NSub*rand(NSub,1))+1;
    else
        Samp1 = floor(NGr1*rand(NGr1,1))+1;
        Samp2 = floor(NGr2*rand(NGr2,1))+1+NGr1;
        Samp = [Samp1; Samp2];
    end
    Resamples(:,i) = Samp;
end