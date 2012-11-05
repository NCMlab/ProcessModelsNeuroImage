function beta = subfnregress(y,design)
design = x2fx(design);
[Q,R] = qr(design,0);
beta = R\(Q'*y);