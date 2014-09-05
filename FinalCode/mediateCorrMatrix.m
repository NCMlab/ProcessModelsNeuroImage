function [Unstand, Stand, T] = mediateCorrMatrix(CorrMatrix,StandardDevations,N)
% Calculate the mediation model effects from a correlation matrix
Nvar = size(CorrMatrix,1);
if Nvar ~= length(StandardDevations)
    error('length of standard deviation vector does not match the size of the correlation matrix');
else
    % convert the correlation matrix to a covariance matrix
    CovMatrix = diag(StandardDevations)*CorrMatrix*diag(StandardDevations);
    % calculate the parameters
    a = CovMatrix(1,2)/CovMatrix(1,1);
    Model2beta = CovMatrix([1 2],3)'/CovMatrix(1:2,1:2);
    b = Model2beta(2);
    cP = Model2beta(1);
    c = CovMatrix(1,3)/CovMatrix(1,1);
    Unstand = {};
    Unstand.a = a;
    Unstand.b = b;
    Unstand.c = c;
    Unstand.cP = cP;
    % calculate the Standardized parameters
    Sa = CorrMatrix(1,2)/CorrMatrix(1,1);
    Model2beta = CorrMatrix([1 2],3)'/CorrMatrix(1:2,1:2);
    Sb = Model2beta(2);
    ScP = Model2beta(1);
    Sc = CorrMatrix(1,3)/CorrMatrix(1,1);
    Stand = {};
    Stand.a = Sa;
    Stand.b = Sb;
    Stand.c = Sc;
    Stand.cP = ScP;

    % Model 1: X -> M
    R21 = Stand.a*CorrMatrix(1,2);
    % Model 2: X, M -> Y
    R22 = Stand.b*CorrMatrix(2,3) + Stand.cP*CorrMatrix(1,3); 
    % This has to be M predicted by X
    R22b = R21;
    % This has to be Y predicted by X
    R22cP = Stand.c*CorrMatrix(1,3);

    Unstand.SEa = (StandardDevations(2)/StandardDevations(1))*sqrt((1-CorrMatrix(1,2)^2)/(N-2));
    Unstand.SEb = (StandardDevations(3)/StandardDevations(2))*sqrt(1/(1-R21))*sqrt((1-R22)/(N-2-1));
    Unstand.SEc = (StandardDevations(3)/StandardDevations(1))*sqrt((1-CorrMatrix(1,3)^2)/(N-2));
    Unstand.SEcP =(StandardDevations(3)/StandardDevations(1))*sqrt(1/(1-R21))*sqrt((1-R22)/(N-2-1));
    
    
    T = {};
    T.a = Unstand.a/Unstand.SEa;
    T.b = Unstand.b/Unstand.SEb;
    T.c = Unstand.c/Unstand.SEc;
    T.cP = Unstand.cP/Unstand.SEcP;
end


