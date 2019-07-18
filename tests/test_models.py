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



SIMPLE_AKT_MODEL_SBML = """
// Created by libAntimony v2.9.4
function NonCompetitiveInhibitionWithKcat(km, ki, kcat, E, n, I, S)
  kcat*E*S/((km + S)*(1 + (I/ki)^n));
end

function MMWithKcat_1(IRS1, Insulin, kIRS1Phos_kcat, kIRS1Phos_km)
  kIRS1Phos_kcat*Insulin*IRS1/(kIRS1Phos_km + IRS1);
end

function MMWithKcat_2(IRS1pS636_639, S6KpT389, kIRS1Dephos_kcat, kIRS1Dephos_km)
  kIRS1Dephos_kcat*S6KpT389*IRS1pS636_639/(kIRS1Dephos_km + IRS1pS636_639);
end

function MMWithKcat_3(Akt, IRS1pS636_639, _kAktPhos_kcat, _kAktPhos_km)
  _kAktPhos_kcat*IRS1pS636_639*Akt/(_kAktPhos_km + Akt);
end

function Function_for_R4(AktpT308, Erk_pT202_Y204, _kAktDephos)
  _kAktDephos*AktpT308*Erk_pT202_Y204;
end

function Function_for_R5(AktpT308, Erk_pT202_Y204, TSC2, _kTSC2Phos_kcat, _kTSC2Phos_ki, _kTSC2Phos_km)
  NonCompetitiveInhibitionWithKcat(_kTSC2Phos_km, _kTSC2Phos_ki, _kTSC2Phos_kcat, AktpT308, 1, Erk_pT202_Y204, TSC2);
end

function MM_1(TSC2pT1462, _kTSC2Dephos_km, _kTSC2Dephos_vmax)
  _kTSC2Dephos_vmax*TSC2pT1462/(_kTSC2Dephos_km + TSC2pT1462);
end

function MMWithKcat_4(AktpT308, PRAS40, kPras40PhosByAkt_kcat, kPras40PhosByAkt_km)
  kPras40PhosByAkt_kcat*AktpT308*PRAS40/(kPras40PhosByAkt_km + PRAS40);
end

function MMWithKcat_5(FourEBP1, TSC2, kFourEBP1Phos_kcat, kFourEBP1Phos_km)
  kFourEBP1Phos_kcat*TSC2*FourEBP1/(kFourEBP1Phos_km + FourEBP1);
end

function MMWithKcat_6(S6K, TSC2, kS6KPhos_kcat, kS6KPhos_km)
  kS6KPhos_kcat*TSC2*S6K/(kS6KPhos_km + S6K);
end

function MMWithKcat_7(Erk, Insulin, kErkPhos_kcat, kErkPhos_km)
  kErkPhos_kcat*Insulin*Erk/(kErkPhos_km + Erk);
end

function Function_for_R15(Erk_pT202_Y204, Feedback, kErkDephos)
  kErkDephos*Erk_pT202_Y204*Feedback;
end

function Function_for_R16(Erk_pT202_Y204, kFeedbackIn)
  kFeedbackIn*Erk_pT202_Y204;
end


model *SimpleAktModel()

  // Compartments and Species:
  compartment Cell;
  species $IRS1 in Cell, IRS1pS636_639 in Cell, $Akt in Cell, AktpT308 in Cell;
  species $TSC2 in Cell, TSC2pT1462 in Cell, $PRAS40 in Cell, PRAS40pT246 in Cell;
  species $S6K in Cell, S6KpT389 in Cell, $FourEBP1 in Cell, FourE_BP1pT37_46 in Cell;
  species $Erk in Cell, Erk_pT202_Y204 in Cell, Feedback in Cell;

  // Assignment Rules:
  IRS1 := IRS1_tot - IRS1pS636_639;
  IRS1_tot := 1.925974 + offset_amount;
  Akt := Akt_tot - AktpT308;
  Akt_tot := 1.241997 + offset_amount;
  TSC2 := TSC2_tot - TSC2pT1462;
  TSC2_tot := 1.136033 + offset_amount;
  PRAS40 := PRAS40_tot - PRAS40pT246;
  PRAS40_tot := 0.981968 + offset_amount;
  S6K := S6K_tot - S6KpT389;
  S6K_tot := 1.330735 + offset_amount;
  FourEBP1 := FourEBP1_tot - FourE_BP1pT37_46;
  FourEBP1_tot := 0.458272 + offset_amount;
  Erk := Erk_tot - Erk_pT202_Y204_obs;
  Erk_tot := 1.305048 + offset_amount;
  Erk_pT202_Y204_obs := Erk_pT202_Y204;
  IRS1pS636_639_obs := IRS1pS636_639;
  AktpT308_obs := AktpT308;
  TSC2pT1462_obs := TSC2pT1462;
  PRAS40pT246_obs := PRAS40pT246;
  S6KpT389_obs := S6KpT389;
  FourE_BP1pT37_46_obs := FourE_BP1pT37_46;

  // Reactions:
  R1: $IRS1 => IRS1pS636_639; Cell*MMWithKcat_1(IRS1, Insulin, kIRS1Phos_kcat, kIRS1Phos_km);
  R2: IRS1pS636_639 => $IRS1; Cell*MMWithKcat_2(IRS1pS636_639, S6KpT389, kIRS1Dephos_kcat, kIRS1Dephos_km);
  R3: $Akt => AktpT308; Cell*MMWithKcat_3(Akt, IRS1pS636_639, _kAktPhos_kcat, _kAktPhos_km);
  R4: AktpT308 => $Akt; Cell*Function_for_R4(AktpT308, Erk_pT202_Y204, _kAktDephos);
  R5: $TSC2 => TSC2pT1462; Cell*Function_for_R5(AktpT308, Erk_pT202_Y204, TSC2, _kTSC2Phos_kcat, _kTSC2Phos_ki, _kTSC2Phos_km);
  R6: TSC2pT1462 => $TSC2; Cell*MM_1(TSC2pT1462, _kTSC2Dephos_km, _kTSC2Dephos_vmax);
  R7: $PRAS40 => PRAS40pT246; Cell*MMWithKcat_4(AktpT308, PRAS40, kPras40PhosByAkt_kcat, kPras40PhosByAkt_km);
  R9: PRAS40pT246 => $PRAS40; Cell*kPras40Dephos*PRAS40pT246;
  R10: $FourEBP1 => FourE_BP1pT37_46; Cell*MMWithKcat_5(FourEBP1, TSC2, kFourEBP1Phos_kcat, kFourEBP1Phos_km);
  R11: FourE_BP1pT37_46 => $FourEBP1; Cell*kFourEBP1Dephos*FourE_BP1pT37_46;
  R12: $S6K => S6KpT389; Cell*MMWithKcat_6(S6K, TSC2, kS6KPhos_kcat, kS6KPhos_km);
  R13: S6KpT389 => $S6K; Cell*kS6KDephos*S6KpT389;
  R14: $Erk => Erk_pT202_Y204; Cell*MMWithKcat_7(Erk, Insulin, kErkPhos_kcat, kErkPhos_km);
  R15: Erk_pT202_Y204 => $Erk; Cell*Function_for_R15(Erk_pT202_Y204, Feedback, kErkDephos);
  R16:  => Feedback; Cell*Function_for_R16(Erk_pT202_Y204, kFeedbackIn);
  R17: Feedback => ; Cell*kFeedbackOut*Feedback;
  IRS1 has substance_per_volume;

  // Species initializations:
  IRS1pS636_639 = 0.861333;
  IRS1pS636_639 has substance_per_volume;
  Akt has substance_per_volume;
  AktpT308 = 0.486243;
  AktpT308 has substance_per_volume;
  TSC2 has substance_per_volume;
  TSC2pT1462 = 0.644957;
  TSC2pT1462 has substance_per_volume;
  PRAS40 has substance_per_volume;
  PRAS40pT246 = 0.38719;
  PRAS40pT246 has substance_per_volume;
  S6K has substance_per_volume;
  S6KpT389 = 0.395656;
  S6KpT389 has substance_per_volume;
  FourEBP1 has substance_per_volume;
  FourE_BP1pT37_46 = 0.488169;
  FourE_BP1pT37_46 has substance_per_volume;
  Erk has substance_per_volume;
  Erk_pT202_Y204 = 0.115661;
  Erk_pT202_Y204 has substance_per_volume;
  Feedback = 0;
  Feedback has substance_per_volume;

  // Compartment initializations:
  Cell = 1;
  Cell has volume;

  // Variable initializations:
  Insulin = 1;
  offset_amount = 1;
  kFeedbackOut = 0.01326991221;
  kFeedbackIn = 1092.520199;
  kFourEBP1Dephos = 0.04760716822;
  kS6KDephos = 0.02513594021;
  kErkDephos = 1.211544484e-05;
  kPras40Dephos = 0.1001194064;
  kErkPhos_km = 348.4143484;
  kFourEBP1Phos_km = 82.73582413;
  kS6KPhos_km = 9711.278764;
  kPras40PhosByAkt_km = 2684.754046;
  kS6KPhos_kcat = 350.115286;
  kErkPhos_kcat = 9.538348056;
  kFourEBP1Phos_kcat = 4.768585613;
  kPras40PhosByAkt_kcat = 286.5647513;
  kIRS1Phos_km = 4469.348507;
  kIRS1Dephos_km = 1732.050115;
  kIRS1Dephos_kcat = 31.97748246;
  kIRS1Phos_kcat = 91.25250293;
  _kAktDephos = 403.462;
  _kAktPhos_km = 23.32;
  _kAktPhos_kcat = 4887.04;
  _kTSC2Phos_kcat = 1285.13;
  _kTSC2Phos_km = 9527.75;
  _kTSC2Phos_ki = 1621.49;
  _kTSC2Dephos_vmax = 1136.83;
  _kTSC2Dephos_km = 7776.66;

  // Other declarations:
  var IRS1_tot, Akt_tot, TSC2_tot, PRAS40_tot, S6K_tot, FourEBP1_tot, Erk_tot;
  var Erk_pT202_Y204_obs, IRS1pS636_639_obs, AktpT308_obs, TSC2pT1462_obs;
  var PRAS40pT246_obs, S6KpT389_obs, FourE_BP1pT37_46_obs;
  const Cell, Insulin, offset_amount, kFeedbackOut, kFeedbackIn, kFourEBP1Dephos;
  const kS6KDephos, kErkDephos, kPras40Dephos, kErkPhos_km, kFourEBP1Phos_km;
  const kS6KPhos_km, kPras40PhosByAkt_km, kS6KPhos_kcat, kErkPhos_kcat, kFourEBP1Phos_kcat;
  const kPras40PhosByAkt_kcat, kIRS1Phos_km, kIRS1Dephos_km, kIRS1Dephos_kcat;
  const kIRS1Phos_kcat, _kAktDephos, _kAktPhos_km, _kAktPhos_kcat, _kTSC2Phos_kcat;
  const _kTSC2Phos_km, _kTSC2Phos_ki, _kTSC2Dephos_vmax, _kTSC2Dephos_km;

  // Unit definitions:
  unit length = metre;
  unit area = metre^2;
  unit volume = litre;
  unit time_unit = time_unit;
  unit substance = mole;
  unit extent = substance;
  unit substance_per_volume = mole / litre;

  // Display Names:
  time_unit is "time";
end
"""

















