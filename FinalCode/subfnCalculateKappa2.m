function k2 = subfnCalculateKappa2(X, M, Y, a, b)
% The following calculation of Kappa2 follows the description given in
% "Effect Size Measures for Mediation Models: Quantitative Strategies for
% Communicating Indirect Effects"
% by Preacher and Kelley, 2011, Psychological Methods
%
% The code is a copy of the R code that accompanies the paper and found at
% (9/20/2012) http://cran.r-project.org/web/packages/MBESS/index.html
%
% Find the covariance matrix of the variables
try
    COV = cov([X M Y]);
    s2 = {};
    s2.X = COV(1,1);
    s2.M = COV(2,2);
    s2.Y = COV(3,3);
    s = {};
    
    s.YX = COV(1,3);
    s.XM = COV(2,1);
    s.YM = COV(3,2);
    % Calculate the permissible values of a
    perm_a =[(s.YM * s.YX + sqrt(s2.M * s2.Y - s.YM^2) * sqrt(s2.X * s2.Y - s.YX^2))/(s2.X * s2.Y), ...
        (s.YM * s.YX - sqrt(s2.M * s2.Y - s.YM^2) * sqrt(s2.X * s2.Y - s.YX^2))/(s2.X * s2.Y)];
    % Calculate the permissible values of b
    perm_b = [sqrt(s2.X * s2.Y - s.YX^2)/sqrt(s2.X * s2.M - s.XM^2), ...
        -sqrt(s2.X * s2.Y - s.YX^2)/sqrt(s2.X * s2.M - s.XM^2)];
    % Check the a values and find the one in teh same direction as the point estimate of a
    if a > 0
        temp = perm_a(find(perm_a>0));
        max_a = max(temp);
    elseif a < 0
        temp = perm_a(find(perm_a < 0));
        max_a = min(temp);
    else
        max_a = 0;
    end
    % Check the b values and find the one in teh same direction as the point
    % estimate of b
    if b > 0
        temp = perm_b(find(perm_b>0));
        max_b = max(temp);
    elseif b < 0
        temp = perm_b(find(perm_b < 0));
        max_b = min(temp);
    else
        max_b = 0;
    end
    % calculate kappa
    %fprintf(1,'%d %d %d %d\n',a, b, max_a, max_b);
    perm_ab = max_a*max_b;
    k2 = (a*b)/perm_ab;
    
catch 
    k2 = 0;
    %fprintf(1,'%d %d %d %d\n',a, b, max_a, max_b);
end
