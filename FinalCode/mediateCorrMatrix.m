function [a, b, c, cP, Sa, Sb, Sc, ScP] = mediateCorrMatrix(CorrMatrix,StandardDevations)
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
    
    % calculate the Standardized parameters
    Sa = CorrMatrix(1,2)/CorrMatrix(1,1);
    Model2beta = CorrMatrix([1 2],3)'/CorrMatrix(1:2,1:2);
    Sb = Model2beta(2);
    ScP = Model2beta(1);
    Sc = CorrMatrix(1,3)/CorrMatrix(1,1);
end


