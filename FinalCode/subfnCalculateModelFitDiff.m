function diffS = subfnCalculateModelFitDiff(S,S1);
diffS = {};
diffS.rsquare = S.rsquare - S1.rsquare;
diffS.adjrsquare = S.adjrsquare - S1.adjrsquare;
diffS.fstat.dfr = S.fstat.dfr - S1.fstat.dfr;
diffS.fstat.dfe = S.fstat.dfe;
diffS.fstat.f = diffS.rsquare/(diffS.fstat.dfr)/((1 - S.rsquare)/diffS.fstat.dfe);
diffS.fstat.pval = fcdf(1/diffS.fstat.f, diffS.fstat.dfe, diffS.fstat.dfr);