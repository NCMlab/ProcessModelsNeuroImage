function beta = subfnregress(y,design)
% Check to see if the y variable is dichotomous
if isequal(logical(y),y)
    beta = subfnLogisticRegress(y,design);
else
    design = x2fx(design);
    [Q,R] = qr(design,0);
    beta = R\(Q'*y);
end