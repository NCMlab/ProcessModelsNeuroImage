function [BCaci PERci p Z] = subfnFindConfidenceIntervals(data,bstat,pointEst,Thresholds)
Sbstat = sort(bstat);
[nboot Nmed NParameters] = size(bstat);
% Find the number of non-zero bstat values
if length(find(squeeze(bstat(1,1,:)))) == 1
    NParameters = 1;
end
BCaci = cell(Nmed,NParameters);
PERci = cell(Nmed,NParameters);

for i = 1:length(data.Thresholds)
    % The p and Z addition does not need to be recalculated at each
    % threshold tested but it is for now.
    [Alpha1 Alpha2 p Z] = subfnFindBCaLimits(bstat,pointEst.values,data.Thresholds(i),data);
    temp = num2str(data.Thresholds(i));

    if ~isnan(Alpha1) | ~isnan(Alpha2)
        for k = 1:NParameters
            for j = 1:Nmed
                % [Sbstat(ceil(Alpha1(j,k)*nboot),j,k) Sbstat(ceil(Alpha2(j,k)*nboot),j,k)]
                BCaci{j,k} = setfield(BCaci{j,k},['alpha' temp(3:end)],[Sbstat(ceil(Alpha1(j,k)*nboot),j,k) Sbstat(ceil(Alpha2(j,k)*nboot),j,k)]);
                PERci{j,k} = setfield(PERci{j,k},['alpha' temp(3:end)],[Sbstat(ceil(data.Thresholds(i)/2*nboot),j,k) Sbstat(ceil((1-data.Thresholds(i)/2)*nboot),j,k)]);

            end
        end
    else
        for k = 1:NParameters
            for j = 1:Nmed
                % [Sbstat(ceil(Alpha1(j,k)*nboot),j,k) Sbstat(ceil(Alpha2(j,k)*nboot),j,k)]
                BCaci{j,k} = setfield(BCaci{j,k},['alpha' temp(3:end)],[-99 -99]);
                PERci{j,k} = setfield(PERci{j,k},['alpha' temp(3:end)],[-99 -99]);
            end
        end
    end
end
