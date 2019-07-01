FUNCTIONS = """

function MM(km, Vmax, S)
        Vmax * S / (km + S)
    end

    function MMWithKcat(km, kcat, S, E)
        kcat * E * S / (km + S)
    end

    function NonCompetitiveInhibition(km, ki, Vmax, n, I, S)
        Vmax * S / ( (km + S) * (1 + (I / ki)^n ) )
    end

    function MA1(k, S)
        k * S
    end

    function MA2(k, S1, S2)
        k * S1 * S2
    end

    function MA1Mod(k, S, M)
        k * S * M
    end

    function MA2Mod(k, S1, S2, M)
        k * S1 * S2 * M
    end

    function CompetitiveInhibitionWithKcat(km, ki, kcat, E, I, S)
        (kcat * E * S) / (km + S + ((km * I )/ ki)  )
    end    

    function CompetitiveInhibition(Vmax, km, ki, I, S)
        Vmax * S / (km + S + ((km * I )/ ki)  )
    end

    function Hill(km, beta, n, X)
        beta * X^n / (km*n + X*n)
    end

    function HillWithKcat(km, kcat, n, X, E)
        kcat*E* X^n / (km*n + X*n)
    end
"""

TEST_MODEL1 = """
model test_model1
    compartment cell = 1
    A in cell
    B in cell 
    C in cell
    
    A = 0;
    B = 0;
    C = 0;
    S = 0;
    
    k1 = 0.1;
    k2 = 0.1;
    k3 = 0.1;
    k4 = 0.1;
    k5 = 0.1;
    k6 = 0.1;
    k7 = 0.1;
    
    R1: => A ; cell * k1;
    R2: A => ; cell * A * k2
    R3: A => B ; cell * k3*A*S;
    R4: B => A ; cell * k4*B;
    R5: B => ; cell * k5*B*C;
    R6: => C ; cell * k6*B;
    R7: C => ; cell * k7*C;
    
end
"""

TEST_MODEL2 = f"""
{FUNCTIONS}
model test_model2
    compartment cell = 1
    A in cell
    B in cell 
    C in cell

    A = 50;
    B = 35;
    C = 10;
    Ap = 0;
    Bp = 0;
    Cp = 0;
    I = 0;
    S = 1;

    r1km        = 60;
    r1kcat      = 0.5;
    r1n         = 2;
    r2km        = 30;
    r2kcat      = 4;
    r3vmax      = 50;
    r3km        = 90;
    r3ki        = 0.5;
    r3kcat      = 0.5;
    r4km        = 30;
    r4vmax      = 5;
    r5km        = 50;
    r5kcat      = 0.5;
    r6km        = 60;
    r6vmax      = 3.5;

    R1 : A => Ap ;  HillWithKcat(r1km, r1kcat, r1n, A, S); 
    R2 : Ap => A ;  MMWithKcat(r2km, r2kcat, Ap, Cp);
    R3 : B => Bp ;  CompetitiveInhibitionWithKcat(r3km, r3ki, r3kcat, Ap, I, B);
    R4 : Bp => B ;  MM(r4km, r4vmax, Bp);
    R5 : C => Cp ;  MMWithKcat(r5km, r5kcat, C, Bp)
    R6 : Cp => C ;  MM(r6km, r6vmax, Cp);

end
"""





















