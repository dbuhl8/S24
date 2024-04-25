

A = [0.0, 10.0, -10.0; -100.0, -1.0, 0.0; 0.0, 10.0, -100.0];

[V, D] = eig(A);
for dt = 10.^-(1:3)
    Current_DT = dt
    for i = 1:3
        lambda_A = D(i, i);
        Current_Lambda = lambda_A
        C = [0.0, 1.0, 0.0; 0.0, 0.0, 1.0; (5.*dt*lambda_A)/12., (-16.0*dt*lambda_A)/12., 1 + (23.*dt*lambda_A)/12.];
        [C1, C2] = eig(C);
        for j = 1:3
            lambda_C(j) = abs(C2(j, j));
        end 
        if (max(lambda_C) <= 1.)
            Current_DT_Lambda = true
        else 
            Current_DT_Lambda = false
        end
        %For_Current_Lambda_DT = lambda_C
    end 
end
