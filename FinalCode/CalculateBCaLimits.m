function [Alpha1, Alpha2, Z, p] = CalculateBCaLimits(JackKnifeData,PointEstimate, BootStrapData,alpha)    

    [m n] = size(PointEstimate);
    N = size(JackKnifeData,3);
    Nboot = size(BootStrapData,3);
    zA = norminv(alpha/2);
    z1mA = norminv(1 - alpha/2);
    Alpha1 = zeros(m,n);
    Alpha2 = zeros(m,n);
    Z = zeros(m,n);
    p = zeros(m,n);
    testValue = 0;
    for i = 1:m
        for j = 1:n
            if PointEstimate(i,j) ~= 0
                % For diagnostic purposes
                % fprintf(1,'i = %d,j = %d\n',i,j);
                zh0 = norminv(length(find(BootStrapData(i,j,:) < PointEstimate(i,j)))/Nboot);
                ThetaDiff = (sum(JackKnifeData(i,j,:))/N) - squeeze(JackKnifeData(i,j,:));
                acc = (sum(ThetaDiff.^3))/(6*(sum(ThetaDiff.^2))^(3/2));
                Alpha1(i,j) = normcdf(zh0 + (zh0+zA)/(1 - acc*(zh0 + zA)));
                Alpha2(i,j) = normcdf(zh0 + (zh0+z1mA)/(1 - acc*(zh0 + z1mA)));
                % Find percentile of the distribution below the null value
                PCTlower = sum(BootStrapData(i,j,:) < testValue)./Nboot;
                PCTupper = sum(BootStrapData(i,j,:) > testValue)./Nboot;
                IsLower = PCTlower < PCTupper;
                
                
                % Check to make sure the tails are not zero
                PCTlower = max(PCTlower,1./Nboot);
                PCTupper = min(PCTupper, 1-1./Nboot);
                
                PCTupper = max(PCTupper,1/Nboot);
                PCTlower = min(PCTlower,1-1/Nboot);
                
                % Z-score
                ZPCTlower = norminv(PCTlower) - zh0;
                ZPCTupper = norminv(PCTupper) - zh0;
                % adjust for acceleration
                Zlower = (ZPCTlower.*(1-acc.*zh0)-zh0)./(1+acc.*ZPCTlower);
                Zupper = (ZPCTupper.*(1-acc.*zh0)-zh0)./(1+acc.*ZPCTupper);
                % which tail are we in?
                
                ZTemp = Zupper;
                ZTemp(IsLower) = -Zlower(IsLower);
                Z(i,j) = ZTemp;
                
                pTemp = normcdf(ZTemp);
                pTemp = [pTemp; 1 - pTemp];
                p(i,j) = 2.*min(pTemp);
            end
        end
        
    end