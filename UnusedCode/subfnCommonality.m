function [C, V, T, r, R2, R2_13, R2_12, R2_23] = subfnCommonality(dependent, independent)
[M N] = size(independent);
C = 0;
V = 0;
T = 0;
R2_13 = [];
R2_12 = [];
R2_23 = [];
switch N
    case 2
        r = corr([independent dependent]);
        R2 = r.^2;
        R2_1 = corr(dependent,independent(:,1))^2;
        R2_2 = corr(dependent,independent(:,2))^2;
        
        COL = [1 2];
        [B BINT R] = regress(dependent,[independent(:,COL) ones(M,1)]);
        pred = [independent(:,COL) ones(M,1)]*B;
        R2_12 = corr(dependent, pred)^2;
        
        V = zeros(N,N);
        V(1,1) = R2_12 - R2_2;
        V(2,2) = R2_12 - R2_1;
        C = R2_1 + R2_2 - R2_12;
        T = R2_12;
        
    case 3
        %%
        r = corr([independent dependent]);
        R2 = r.^2;
        % Total PVAF between depedent and each independent
        R2_1 = corr(dependent,independent(:,1))^2;
        R2_2 = corr(dependent,independent(:,2))^2;
        R2_3 = corr(dependent,independent(:,3))^2;
        
        COL = [1 2 3];
        [B BINT R] = regress(dependent,[independent(:,COL) ones(M,1)]);
        pred = [independent(:,COL) ones(M,1)]*B;
        R2_123 = corr(dependent, pred)^2;
        
        COL = [1 2];
        [B] = regress(dependent,[independent(:,COL) ones(M,1)]);
        pred = [independent(:,COL)  ones(M,1)]*B;
        R2_12 = corr(dependent, pred)^2;
        
        COL = [1 3];
        [B] = regress(dependent,[independent(:,COL) ones(M,1)]);
        pred = [independent(:,COL)  ones(M,1)]*B;
        R2_13 = corr(dependent, pred)^2;
        
        COL = [2 3];
        [B] = regress(dependent,[independent(:,COL) ones(M,1)]);
        pred = [independent(:,COL)  ones(M,1)]*B;
        R2_23 = corr(dependent, pred)^2;
        
        V = zeros(1,N);
        V(1,1) = R2_123 - R2_23;
        V(2,2) = R2_123 - R2_13;
        V(3,3) = R2_123 - R2_12;
        V(1,2) = R2_13 + R2_23 - R2_3 - R2_123;
        V(1,3) = R2_12 + R2_23 - R2_2 - R2_123;
        V(2,3) = R2_12 + R2_13 - R2_1 - R2_123;
        V = V;
        C = R2_1 + R2_2 + R2_3 - R2_12 - R2_13 - R2_23 + R2_123;
        C = C;
        T = R2_123;
        
    case 5
        
        r = corr([independent dependent]);
        R2 = r.^2;
        R2_1 = corr(dependent,independent(:,1))^2;
        R2_2 = corr(dependent,independent(:,2))^2;
        R2_3 = corr(dependent,independent(:,3))^2;
        R2_4 = corr(dependent,independent(:,4))^2;
        R2_5 = corr(dependent,independent(:,5))^2;
        
        COL = [1 2 3 4 5];
        [B BINT R] = regress(dependent,[independent(:,COL) ones(M,1)]);
        pred = [independent(:,COL) ones(M,1)]*B;
        R2_12345 = corr(dependent, pred)^2;
        
        COL = [1 2 3 4];
        [B BINT R] = regress(dependent,[independent(:,COL) ones(M,1)]);
        pred = [independent(:,COL) ones(M,1)]*B;
        R2_1234 = corr(dependent, pred)^2;
        
        COL = [1 2 3 5];
        [B BINT R] = regress(dependent,[independent(:,COL) ones(M,1)]);
        pred = [independent(:,COL) ones(M,1)]*B;
        R2_1235 = corr(dependent, pred)^2;
        
        COL = [1 2 4 5];
        [B BINT R] = regress(dependent,[independent(:,COL) ones(M,1)]);
        pred = [independent(:,COL) ones(M,1)]*B;
        R2_1245 = corr(dependent, pred)^2;
        
        COL = [1 3 4 5];
        [B BINT R] = regress(dependent,[independent(:,COL) ones(M,1)]);
        pred = [independent(:,COL) ones(M,1)]*B;
        R2_1345 = corr(dependent, pred)^2;
        
        COL = [2 3 4 5];
        [B BINT R] = regress(dependent,[independent(:,COL) ones(M,1)]);
        pred = [independent(:,COL) ones(M,1)]*B;
        R2_2345 = corr(dependent, pred)^2;

        COL = [1 2 3];
        [B BINT R] = regress(dependent,[independent(:,COL) ones(M,1)]);
        pred = [independent(:,COL) ones(M,1)]*B;
        R2_123 = corr(dependent, pred)^2;
        
          COL = [3 4 5];
        [B BINT R] = regress(dependent,[independent(:,COL) ones(M,1)]);
        pred = [independent(:,COL) ones(M,1)]*B;
        R2_345 = corr(dependent, pred)^2;
        
        COL = [2 4 5];
        [B BINT R] = regress(dependent,[independent(:,COL) ones(M,1)]);
        pred = [independent(:,COL) ones(M,1)]*B;
        R2_245 = corr(dependent, pred)^2;
              
        COL = [2 3 5];
        [B BINT R] = regress(dependent,[independent(:,COL) ones(M,1)]);
        pred = [independent(:,COL) ones(M,1)]*B;
        R2_235 = corr(dependent, pred)^2; 
        
        COL = [2 3 4];
        [B BINT R] = regress(dependent,[independent(:,COL) ones(M,1)]);
        pred = [independent(:,COL) ones(M,1)]*B;
        R2_234 = corr(dependent, pred)^2; 
        
        COL = [1 4 5];
        [B BINT R] = regress(dependent,[independent(:,COL) ones(M,1)]);
        pred = [independent(:,COL) ones(M,1)]*B;
        R2_145 = corr(dependent, pred)^2; 
        
        COL = [1 3 5];
        [B BINT R] = regress(dependent,[independent(:,COL) ones(M,1)]);
        pred = [independent(:,COL) ones(M,1)]*B;
        R2_135 = corr(dependent, pred)^2; 
        
        COL = [1 3 4];
        [B BINT R] = regress(dependent,[independent(:,COL) ones(M,1)]);
        pred = [independent(:,COL) ones(M,1)]*B;
        R2_134 = corr(dependent, pred)^2; 
        
        COL = [1 2 5];
        [B BINT R] = regress(dependent,[independent(:,COL) ones(M,1)]);
        pred = [independent(:,COL) ones(M,1)]*B;
        R2_125 = corr(dependent, pred)^2; 
        
         COL = [1 2 4];
        [B BINT R] = regress(dependent,[independent(:,COL) ones(M,1)]);
        pred = [independent(:,COL) ones(M,1)]*B;
        R2_124 = corr(dependent, pred)^2; 
        
        COL = [1 2 3];
        [B BINT R] = regress(dependent,[independent(:,COL) ones(M,1)]);
        pred = [independent(:,COL) ones(M,1)]*B;
        R2_123 = corr(dependent, pred)^2; 
        
        V = zeros(N,N);
        V(1,1) = R2_12345 - R2_2345;
        V(2,2) = R2_12345 - R2_1345;
        V(3,3) = R2_12345 - R2_1245;
        V(4,4) = R2_12345 - R2_1235;
        V(5,5) = R2_12345 - R2_1234;
        V(1,2) = R2_1345 + R2_2345 - R2_345 - R2_12345;
        V(1,3) = R2_1245 + R2_2345 - R2_245 - R2_12345;
        V(1,4) = R2_1235 + R2_2345 - R2_235 - R2_12345;
        V(1,5) = R2_1234 + R2_2345 - R2_234 - R2_12345;
        V(2,3) = R2_1245 + R2_1345 - R2_145 - R2_12345;
        V(2,4) = R2_1235 + R2_1345 - R2_135 - R2_12345;
        V(2,5) = R2_1234 + R2_1345 - R2_134 - R2_12345;
        V(3,4) = R2_1235 + R2_1245 - R2_125 - R2_12345;
        V(3,5) = R2_1234 + R2_1245 - R2_124 - R2_12345;
        V(4,5) = R2_1234 + R2_1235 - R2_123 - R2_12345;
        %V = V*100;
        C = -99;
        T = R2_12345;
    
end
    
    