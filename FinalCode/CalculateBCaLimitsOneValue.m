function [Alpha1, Alpha2, Z, p,  BCaCI] = CalculateBCaLimitsOneValue(JackKnifeData,PointEstimate, BootStrapData,alpha)


N = size(JackKnifeData,1); % NUMBER OF SUBJECTS
Nboot = size(BootStrapData,1);

zA = norminv(alpha/2);
z1mA = norminv(1 - alpha/2);

testValue = 0;

if PointEstimate ~= 0%{1}(i,j) ~= 0
    % For diagnostic purposes
    % fprintf(1,'i = %d,j = %d\n',i,j);
    zh0 = norminv(length(find(BootStrapData < PointEstimate))/Nboot);
    ThetaDiff = (sum(JackKnifeData)/N) - squeeze(JackKnifeData);
    acc = (sum(ThetaDiff.^3))/(6*(sum(ThetaDiff.^2))^(3/2));
    Alpha1= normcdf(zh0 + (zh0+zA)/(1 - acc*(zh0 + zA)));
    Alpha2 = normcdf(zh0 + (zh0+z1mA)/(1 - acc*(zh0 + z1mA)));
    PCTlower = prctile(BootStrapData,Alpha1*100,1);
    PCTupper = prctile(BootStrapData,Alpha2*100,1);
    BCaCI = [PCTlower PCTupper];
    
    % Find percentile of the distribution below the null value
    PCTlower = sum(BootStrapData < testValue)./Nboot;
    PCTupper = sum(BootStrapData > testValue)./Nboot;


        
        % Check to make sure the tails are not zero
        PCTlower = max(PCTlower,1./Nboot);
        PCTupper = min(PCTupper, 1-1./Nboot);
        
        PCTupper = max(PCTupper,1/Nboot);
        PCTlower = min(PCTlower,1-1/Nboot);
        
        % Z-score
        ZPCTlower = norminv(PCTlower) - zh0;
        ZPCTupper = norminv(PCTupper) - zh0;
        % adjust for acceleration
        Z_1 = (ZPCTlower.*(1-acc.*zh0)-zh0)./(1+acc.*ZPCTlower);
        Z_2 = (ZPCTupper.*(1-acc.*zh0)-zh0)./(1+acc.*ZPCTupper);

        % which tail are we in?
        % If there is more of the distribution to the left of zero, then we
        % are in the lower tail
        if PCTlower < PCTupper
            Z = Z_2;
        else
            % otherwise we are in the other tail
            Z = -Z_1;
        end
        
        pTemp = normcdf(Z);
        pTemp = [pTemp; 1 - pTemp];
        p = 2.*min(pTemp);

end


    