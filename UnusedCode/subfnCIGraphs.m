function subfnCIGraphs(data,Parameters)
ModelNum = data.ModelNum;
switch ModelNum
    case '14'
        A = data.X;
        B = data.M;
        V = data.V;
        a = Parameters.Model1{1}.AgeGroup.beta;
        
        b = Parameters.Model2.FreeSurfer.beta;
        v = Parameters.Model2.FreeSurfer_x_CogRes.beta;
        PE = a.*(b+v.*V);
        COV = cov([A B V]);
        VARa = COV(1,1);
        VARb = COV(2,2);
        VARv = COV(3,3);
        COVbv = COV(2,3);
        SOVAR = VARa.*(b + v.*V).^2 + (a^2+VARa).*(VARb+2*sqrt(COVbv).*V+VARv.*V.*V);
        
        figure(1)
        clf
        hold on
        plot(V,PE)
        plot(V,PE+1.96.*sqrt(SOVAR))
        plot(V,PE-1.96.*sqrt(SOVAR))
        line([-3 3],[ 0 0])