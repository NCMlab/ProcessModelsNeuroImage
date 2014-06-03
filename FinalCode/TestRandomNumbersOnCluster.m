function TestRandomNumbersOnCluster

rng('shuffle','multFibonacci');
D = randperm(39);
fprintf(1,'%d ',D);




