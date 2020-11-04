#include "math.h"
#include "udf.h"
#define  prim_index 0  /*index of primary phase*/
#define  sec1_index 1  /*index of primary phase*/

#define R 8.314
#define pi 3.141592
double alphax = 0.1, MWP[20], MW[25], Mol_Weight, wH2 = -0.215993, wCO = 0.045, wCO2 = 0.223621, wH2O = 0.344, wP[20], w[24], long_term[24][25], fugacitycoeff[24], react3, prod3;
double ac[24], k[24], r[24], a, b, P, T, Tr[24], Pr[24], A, B, V_gas, V_liq, D[24], y[24],x[24], alpha[24], Zv, Zl, rho_gas1, rho_liquid, rhol, asum, bsum, change,dfdz, m[24][24], tot_gas_moles;
double  P_ref, T_ref, aCO, bCO, cCO, dCO, Hf298CO, Gf298CO, GammaCO, GfTCO, aH2, bH2, cH2, dH2, Hf298H2, Gf298H2, GammaH2, GfTH2, aH2O, bH2O, cH2O, dH2O, Hf298H2O, Gf298H2O, GammaH2O;
double GammaP[6], GftP[6], Hft298P[6], GftRxns[6], K_equilibriums[6], Hf298P[6], Gf298P[6], GfTP[6], GfTrxnP[7], GfTrxn1, GfTrxn2, prod[18], react[18], rate1, rate2, activity[24], rate[20];
double  GfTH2O, aCH4, bCH4, cCH4, dCH4, Hf298CH4, Gf298CH4, GammaCH4, GfTCH4, aC2H6, bC2H6, cC2H6, dC2H6, GfTCH2, kf[6], MWsyngas[4] = { 2.01588,28.0101,44.0095,18.0153 }, rate3;
double prod1, react1, prod2, react2, K_equilibriumCH2, K_equilibriumCH4, K_equilibriumC2H6, activityCH2, density_catalyst, porosity_catalyst, ACH2, BCH2, C, K_equilibriumP[6];
double X, omegaA, omegaB, alpha0[24], alpha1[24], alphaT[24], ai[24], bi[24], Q[12][12], q[24], I, lnphi[24], phi[24],a0,b0,cube;
double value, f[300], dfdz, change, Zv, V_gas, V_liq, p0, q0, r0, phi0, a, b, rho_liquid, rhol, Zl, y1, y2, y3, z[300],test_term;
double PcH2 = 1313000.0, PcCO = 3494000.0, PcCO2 = 7383000.0, PcH2O = 22064000.0, PcP[20], TcP[20], TcH2 = 33.19, TcCO = 132.85, TcCO2 = 304.21, TcH2O = 647.14, Tc[24], Pc[24];
double acube, rootacube, bbya, twointo, phinew;
int n, s, p, i, j, flag, n_iterations,pt;

DEFINE_ADJUST(MY_FUNC, md)
{
   cell_t c;
    Thread* t;
    Domain* d;
    sub_domain_loop(d, md, pt)
    {
        if (DOMAIN_ID(d) == 2)
            thread_loop_c(t, d)
        {
            begin_c_loop_all(c, t)
            {
                for (n = 0; n < 20; n++)
                {
                    MWP[n] = 2.01590803640286 + 14.026579298921 * (n + 1);
                }
                for (n = 0; n < 4; n++)
                {
                    MW[n] = MWsyngas[n];
                }
                for (n = 0; n < 20; n++)
                {
                    MW[n + 4] = MWP[n];
                }

                  tot_gas_moles = C_YI(c, t, 0) / MW[0] + C_YI(c, t, 1) / MW[1] + C_YI(c, t, 2) / MW[2] + C_YI(c, t, 3) / MW[3] + C_YI(c, t, 4) / MW[4] + C_YI(c, t, 5) / MW[5] + C_YI(c, t, 6) / MW[6] + C_YI(c, t, 7) / MW[7] + C_YI(c, t, 8) / MW[8] + C_YI(c, t, 9) / MW[9] + C_YI(c, t, 10) / MW[10] + C_YI(c, t, 11) / MW[11]+ C_YI(c, t, 12) / MW[12] + C_YI(c, t, 13) / MW[13] + C_YI(c, t, 14) / MW[14] + C_YI(c, t, 15) / MW[15] + C_YI(c, t, 16) / MW[16] + C_YI(c, t, 17) / MW[17] + C_YI(c, t, 18) / MW[18] + C_YI(c, t, 19) / MW[19] + C_YI(c, t, 20) / MW[20] + C_YI(c, t, 21) / MW[21] + C_YI(c, t, 22) / MW[22] + C_YI(c, t, 23) / MW[23];;
              //    printf("\nweight frac: C_YI H2,CO,CO2,H2O,C1-C8 = %f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f", C_YI(c, t, 0), C_YI(c, t, 1), C_YI(c, t, 2), C_YI(c, t, 3), C_YI(c, t, 4), C_YI(c, t, 5), C_YI(c, t, 6), C_YI(c, t, 7), C_YI(c, t, 8), C_YI(c, t, 9), C_YI(c, t, 10), C_YI(c, t, 11));

                y[0] = (C_YI(c, t, 0) / MW[0]) / tot_gas_moles;                            /*H2*/      
                y[1] = (C_YI(c, t, 1) / MW[1]) / tot_gas_moles;                            /*CO*/       
                y[2] = (C_YI(c, t, 2) / MW[2]) / tot_gas_moles;                            /*CO2*/      
                y[3] = (C_YI(c, t, 3) / MW[3]) / tot_gas_moles;                            /*H2O*/      
                y[4] = (C_YI(c, t, 4) / MW[4]) / tot_gas_moles;                            /*CH4*/       
                y[5] = (C_YI(c, t, 5) / MW[5]) / tot_gas_moles;                            /*C2*/          
                y[6] = (C_YI(c, t, 6) / MW[6]) / tot_gas_moles;                            /*C3*/        
                y[7] = (C_YI(c, t, 7) / MW[7]) / tot_gas_moles;                            /*C4*/            
                y[8] = (C_YI(c, t, 8) / MW[8]) / tot_gas_moles;                            /*C5*/           
                y[9] = (C_YI(c, t, 9) / MW[9]) / tot_gas_moles;                            /*C6*/        
                y[10] = (C_YI(c, t, 10) / MW[10]) / tot_gas_moles;                         /*C7*/          
                y[11] = (C_YI(c, t, 11) / MW[11]) / tot_gas_moles;                          /*C8*/
                y[12] = (C_YI(c, t, 12) / MW[12]) / tot_gas_moles;                            /*H2*/
                y[13] = (C_YI(c, t, 13) / MW[13]) / tot_gas_moles;                            /*CO*/
                y[14] = (C_YI(c, t, 14) / MW[14]) / tot_gas_moles;                            /*CO2*/
                y[15] = (C_YI(c, t, 15) / MW[15]) / tot_gas_moles;                            /*H2O*/
                y[16] = (C_YI(c, t, 16) / MW[16]) / tot_gas_moles;                            /*CH4*/
                y[17] = (C_YI(c, t, 17) / MW[17]) / tot_gas_moles;                            /*C2*/
                y[18] = (C_YI(c, t, 18) / MW[18]) / tot_gas_moles;                            /*C3*/
                y[19] = (C_YI(c, t, 19) / MW[19]) / tot_gas_moles;                            /*C4*/
                y[20] = (C_YI(c, t, 20) / MW[20]) / tot_gas_moles;                            /*C5*/
                y[21] = (C_YI(c, t, 21) / MW[21]) / tot_gas_moles;                            /*C6*/
                y[22] = (C_YI(c, t, 22) / MW[22]) / tot_gas_moles;                         /*C7*/
                y[23] = (C_YI(c, t, 23) / MW[23]) / tot_gas_moles;
             //  printf("\ntot_moles & mole frac: yH2,CO,CO2,H2O,C1-C8 = %f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f", tot_gas_moles,y[0], y[1], y[2], y[3], y[4], y[5], y[6], y[7], y[8], y[9], y[10], y[11]);
               C_UDMI(c, t, 0) = y[0];
               C_UDMI(c, t, 1) = y[1];
               C_UDMI(c, t, 2) = y[2];
               C_UDMI(c, t, 3) = y[3];
               C_UDMI(c, t, 4) = y[4];
               C_UDMI(c, t, 5) = y[5];
               C_UDMI(c, t, 6) = y[6];
               C_UDMI(c, t, 7) = y[7];
               C_UDMI(c, t, 8) = y[8];
               C_UDMI(c, t, 9) = y[9];
               C_UDMI(c, t, 10) = y[10];
               C_UDMI(c, t, 11) = y[11];

               PcP[0] = 4599000;
                for (n = 0; n < 19; n++)
                {
                    PcP[n + 1] = (0.418739207076228 - (0.418739207076228 - 77.8476928857587) * exp(-1.35147456886766 * pow((n + 2 + 4.97907115750929), 0.384275980761449))) * 1000000;
                }
                for (n = 0; n < 20; n++)
                {
                    TcP[n] = 956.881457325687 - (956.881457325687 - 181.630935936582) * exp(-0.169720243192845 * pow((n + 1 - 0.975912439965372), 0.720233680827972));
                }
                Tc[0] = TcH2;  Tc[1] = TcCO;  Tc[2] = TcCO2;  Tc[3] = TcH2O;
                Pc[0] = PcH2;  Pc[1] = PcCO;  Pc[2] = PcCO2;  Pc[3] = PcH2O;

                wP[0] = 0.011; wP[1] = 0.099; 
                for (n = 0; n < 18; n++)
                {
                    wP[n + 2] = -16.8218079682128 + 3.58424535476472 * pow((n + 3 + 115.312654179984), 0.325781056639631);
                }
                w[0] = wH2;  w[1] = wCO;  w[2] = wCO2;  w[3] = wH2O;

                for (n = 0; n < 20; n++)
                {
                    Pc[n + 4] = PcP[n];
                    Tc[n + 4] = TcP[n];
                    w[n + 4] = wP[n];
                }
               
                
                T = 493.15;
                P = 2000000;
                X = 0.253076587;
                omegaA = 8 * (5 * X + 1) / (49 - 37 * X);
                omegaB = X / (X + 3);

                for (n = 0; n < 24; n++)
                {
                    ac[n] = omegaA * pow((R * Tc[n]),2)/ Pc[n]; /* attraction term; pressure in Pa (Pa.m^6/mol^2)*/
                    Tr[n] = T/ Tc[n];
                    alpha0[n] = pow(Tr[n],(1.948150 * (0.911807 - 1)))*exp(0.125283 * (1 - pow(Tr[n],(1.948150 * 0.911807))));
                    alpha1[n] = pow(Tr[n],(2.812520 * (0.784054 - 1)))*exp(0.511614 * (1 - pow(Tr[n],(2.812520 * 0.784054))));
                    if (Tr[n] > 1)
                    {
                        alpha0[n] = pow(Tr[n],(-0.2 * (4.963070 - 1)))*exp(0.401219 * (1 - pow(Tr[n],(-0.2 * 4.963070))));
                        alpha1[n] = pow(Tr[n],(-8.0 * (1.248089 - 1)))*exp(0.024955 *(1 -  pow(Tr[n],(-8.0 * 1.248089))));
                    }
                    alphaT[n] = alpha0[n] + w[n]*(alpha1[n] - alpha0[n]);
                    ai[n] = alphaT[n]*ac[n]; /* (Pa.m^6/mol^2)*/
                    bi[n] = omegaB * (R * Tc[n])/ Pc[n]; /* repulsive term; (m^3/mol)*/
                    C_UDMI(c, t, 160 + n) = bi[n];
                    C_UDMI(c, t, 185 + n) = ai[n];
                }
            /*    printf("\n The ac[%d],Tr[%d],alpha0[%d],alpha1[%d], alphaT[%d], ai[%d] and bi[%d]: %f,%f,%f,%f,%f,%f,%f", 1, 1, 1, 1, 1, 1, 1, ac[1], Tr[1], alpha0[1], alpha1[1], alphaT[1], ai[1], bi[1]);*/

                asum=0.0; 
                bsum=0.0;
                for (n = 0; n < 24; n++)
                {
                    for (p = 0; p < 24; p++)
                    {
                        m[n][p] = pow((ai[n] * ai[p]), 0.5);

                        asum=asum+ y[n] * y[p] * m[n][p];
                    }
                    bsum = bsum+ y[n] * bi[n];
                }
                for (p = 0; p < 24; p++)
                {
                    C_UDMI(c, t, 91+p) = m[0][p];
                }
               
                a = asum; b = bsum;

             //   printf("\n the a = %f and b = %f", a,b);

                C_UDMI(c, t, 25) = a;
                C_UDMI(c, t, 26) = b;

                Mol_Weight = 0;
                for (n = 0; n < 24; n++)
                {
                    Mol_Weight = Mol_Weight + MW[n] * y[n];
                }
                A = a * P / (R * R * T * T);
                B = b * P / (R * T);
              //  printf("\n the Pressure = %f and Temperature = %f", P, T);


                C_UDMI(c, t, 64) = A;
                C_UDMI(c, t, 65) = B;
                C_UDMI(c, t, 271) = C_VOF(c, t);
             
              //  Zv= y1; 
           
               n_iterations = N_ITER;
               if (n_iterations < 800)
               {
                   C_UDMI(c, t, 159) = 1;

                   z[0] = 1.3;

                   for (n = 0; n < 300; n++)
                   {
                       f[n] = pow(z[n], 3) - (1 - B) * (pow(z[n], 2)) + (A - 3 * pow(B, 2) - 2 * B) * z[n] - (A * B - B * B - pow(B, 3));
                       dfdz = 3 * pow(z[n], 2) - 2 * (1 - B) * z[n] + (A - 3 * pow(B, 2) - 2 * B);
                       z[n + 1] = z[n] - alphax * (f[n] / dfdz);

                       change = (z[n + 1] - z[n]) / z[n + 1];
                       if (fabs(change) < 0.0001 && (0 < f[n]))
                       {
                           break;
                       }
                   }

                   Zv = z[n + 1];
                   I = 1.0 / (2.0 * sqrt(2.0)) * log((Zv + (1.0 + sqrt(2.0)) * B) / (Zv + (1.0 - sqrt(2.0)) * B));
                   C_UDMI(c, t, 63) = I;
                   C_UDMI(c, t, 116) = n + 1;
                   C_UDMI(c, t, 273) = z[n + 1];
                   C_UDMI(c, t, 274) = Zv;

                   V_gas = z[n + 1] * 8.314 * T / P;
                   rho_gas1 = Mol_Weight / (V_gas * 1000);
                   C_UDMI(c, t, 22) = rho_gas1;

                   z[0] = 0.15;

                   for (n = 0; n < 300; n++)
                   {
                       f[n] = pow(z[n], 3) - (1 - B) * (pow(z[n], 2)) + (A - 3 * pow(B, 2) - 2 * B) * z[n] - (A * B - B * B - pow(B, 3));
                       dfdz = 3 * pow(z[n], 2) - 2 * (1 - B) * z[n] + (A - 3 * pow(B, 2) - 2 * B);
                       z[n + 1] = z[n] - alphax * (f[n] / dfdz);
                       change = (z[n + 1] - z[n]) / z[n + 1];
                       if (fabs(change) < 0.0001)
                       {
                           break;
                       }
                   }
                   Zl = z[n + 1];
                   
                   V_liq = z[n + 1] * 8.314 * T / P;
                   rho_liquid = Mol_Weight / (V_liq * 1000);     rhol = 1 / V_liq;    /* rho = mol/volume (0.028010099kg/m3)*/
               }
               else if (n_iterations > 800)
               {
                   C_UDMI(c, t, 159) = 0;
                   Zv = 1.0;
               }

               C_UDMI(c, t, 23) = Zl;
                /*C_UDMI(c, t, 24) = rho_liquid;*/
               
                for (i = 0; i < 24; i++)
                {
                    long_term[i][0] = 0.0;
                }

                for (i = 0; i < 24; i++)
                {
                    for (j = 0; j < 24; j++)
                    {
                        long_term[i][j + 1] = long_term[i][j] + m[i][j] * y[j];
                    }
                    C_UDMI(c, t, i + 39) = long_term[i][12];
                }
               
               
              //  Zv = 0.97;
             
                for (i = 0; i < 24; i++)
                {
                    q[i] = (A / B) * (2 * long_term[i][12] / A - bi[i] / b);
                    lnphi[i] = (bi[i] / b) * (Zv - 1) - log(Zv - B) - q[i] * I;
                    phi[i] = exp(lnphi[i]);
                 //   C_UDMI(c, t, i + 51) = q[i];
                    C_UDMI(c, t, i + 27) = phi[i];
                }
               


               T_ref = 298.15;
                P_ref = 1.0e5;                                              /*Defining reference temperatures and pressures*/
                aCO = 30.08625492918;
                bCO = -0.0080311809989;
                cCO = 0.00001827934;
                dCO = -0.00000000676765;
                Hf298CO = -110530.0;                                                               /*Defining Ideal Gas Enthalpy of CO at 298K*/
                Gf298CO = -137150.0;                                                                /*Defining Ideal Gas Gibbs Free Energy of CO at 298K*/
                GammaCO = Hf298CO - 30.0862549291816950 * T_ref - 0.5 * -0.0080311809989949 * pow(T_ref, 2.0) - 0.333 * 0.0000182793417809 * pow(T_ref, 3.0) - 0.25 * (-0.0000000067676554) * pow(T_ref, 4.0);
                GfTCO = Gf298CO / T_ref - (-GammaCO * (1.0 / T - 1.0 / T_ref) + aCO * log(T / T_ref) + 0.5 * bCO * (T - T_ref) + 0.1666 * cCO * (pow(T, 2.0) - pow(T_ref, 2.0)) + 0.0833 * dCO * (pow(T, 3.0) - pow(T_ref, 3.0)));

                aH2 = 22.37072390439;
                bH2 = 0.0372441295;
                cH2 = -0.00006590339;
                dH2 = 0.0000000393527806;                                                       /*Defining Ideal Gas Enthalpy of H2 at 298K*/
                Hf298H2 = 0.0;                                                                  /*Defining Ideal Gas Gibbs Free Energy of H2 at 298K*/
                Gf298H2 = 0.0;
                GammaH2 = Hf298H2 - 22.3707239043913050 * T_ref - 0.5 * 0.0372441295211312 * pow(T_ref, 2.0) - 0.333 * -0.0000659033918954 * pow(T_ref, 3.0) - 0.25 * 0.0000000393527806 * pow(T_ref, 4.0);
                GfTH2 = Gf298H2 / T_ref - (-GammaH2 * (1.0 / T - 1.0 / T_ref) + aH2 * log(T / T_ref) + 0.5 * bH2 * (T - T_ref) + 0.1666 * cH2 * (pow(T, 2.0) - pow(T_ref, 2.0)) + 0.0833 * dH2 * (pow(T, 3.0) - pow(T_ref, 3.0)));

                aH2O = 33.5359370366856100;
                bH2O = -0.0050387621184615;
                cH2O = 0.0000201344912058;
                dH2O = -0.0000000075171861;
                Hf298H2O = -241814;                                                         /*Defining Ideal Gas Enthalpy of H2O at 298K*/
                Gf298H2O = -228590;                                                            /*Defining Ideal Gas Gibbs Free Energy of H2O at 298K*/
                GammaH2O = Hf298H2O - aH2O * T_ref - 0.5 * bH2O * pow(T_ref, 2.0) - 0.333 * cH2O * pow(T_ref, 3.0) - 0.25 * dH2O * pow(T_ref, 4.0);
                GfTH2O = Gf298H2O / T_ref - (-GammaH2O * (1.0 / T - 1.0 / T_ref) + aH2O * log(T / T_ref) + 0.5 * bH2O * (T - T_ref) + 0.1666 * cH2O * (pow(T, 2.0) - pow(T_ref, 2.0)) + 0.08333 * dH2O * (pow(T, 3.0) - pow(T_ref, 3.0)));

                aCH4 = -0.764102327;
                bCH4 = 0.094014712;
                cCH4 = -0.0000530366688;
                dCH4 = 0.0000000116898997;
                Hf298CH4 = -74520;                                                          /*Defining Ideal Gas Enthalpy of CH4 at 298K*/
                Gf298CH4 = -50490;                                                          /*Defining Ideal Gas Gibbs Free Energy of CH4 at 298K*/
                GammaCH4 = Hf298CH4 - aCH4 * 1.353922078413607 * T_ref - 0.5 * bCH4 * 1.353922078413607 * pow(T_ref, 2.0) - 0.333 * 1.353922078413607 * cCH4 * pow(T_ref, 3.0) - 0.25 * dCH4 * 1.353922078413607 * pow(T_ref, 4.0);
                GfTCH4 = Gf298CH4 / T_ref - (-GammaCH4 * (1.0 / T - 1.0 / T_ref) + aCH4 * 1.353922078413607 * log(T / T_ref) + 0.5 * bCH4 * 1.353922078413607 * (T - T_ref) + 0.1666 * cCH4 * 1.353922078413607 * (pow(T, 2.0) - pow(T_ref, 2.0)) + 0.08333 * dCH4 * 1.353922078413607 * (pow(T, 3.0) - pow(T_ref, 3.0)));

                aC2H6 = -0.764102327;
                bC2H6 = 0.094014712;
                cC2H6 = -0.0000530366688;
                dC2H6 = 0.0000000116898997;
               
                for (i = 0; i < 2; i++)
                {
                    Hf298P[i] = (-20.7068211128167 * (i + 4.05905037969529)) * 1000;
                    Gf298P[i] = (8.22925109627391 * (i - 4.0080592634736)) * 1000;
                }

                GfTP[0] = GfTCH4;

                GammaP[0] = Hf298P[0] - aC2H6 * (2.353922078413607) * T_ref - 0.5 * bC2H6 * (2.353922078413607) * pow(T_ref, 2.0) - 0.3333 * cC2H6 * (2.353922078413607) * pow(T_ref, 3.0) - 0.25 * dC2H6 * (2.353922078413607) * pow(T_ref, 4.0);
                GfTP[1] = Gf298P[0] / T_ref - (-GammaP[0] * (1.0 / T - 1.0 / T_ref) + aC2H6 * (2.353922078413607) * log(T / T_ref) + 0.5 * bC2H6 * (2.353922078413607) * (T - T_ref) + 0.1666 * cC2H6 * (2.353922078413607) * (pow(T, 2.0) - pow(T_ref, 2.0)) + 0.0833 * dC2H6 * (2.353922078413607) * (pow(T, 3.0) - pow(T_ref, 3.0)));

                GammaP[1] = Hf298P[1] - aC2H6 * (3.353922078413607) * T_ref - 0.5 * bC2H6 * (3.353922078413607) * pow(T_ref, 2.0) - 0.3333 * cC2H6 * (3.353922078413607) * pow(T_ref, 3.0) - 0.25 * dC2H6 * (3.353922078413607) * pow(T_ref, 4.0);
                GfTP[2] = Gf298P[1] / T_ref - (-GammaP[1] * (1.0 / T - 1.0 / T_ref) + aC2H6 * (3.353922078413607) * log(T / T_ref) + 0.5 * bC2H6 * (3.353922078413607) * (T - T_ref) + 0.1666 * cC2H6 * (3.353922078413607) * (pow(T, 2.0) - pow(T_ref, 2.0)) + 0.0833 * dC2H6 * (3.353922078413607) * (pow(T, 3.0) - pow(T_ref, 3.0)));

                GfTP[0] = GfTCH4;

                GfTCH2 = GfTP[2] - GfTP[1];                               /*propane-ethane*/                       /* Defining Gibbs free Energy of CH2*/
                GfTrxn1 = GfTH2O + GfTCH2 - GfTCO - 2.0 * GfTH2;                              /* Defining Gibbs free Energy of Reaction 1*/
                K_equilibriumCH2 = pow(2.71828, (-GfTrxn1 / R));                                       /*Defining Kequillibrium for Reaction 1 */

                GfTrxn2 = GfTCH4 - GfTCH2 - GfTH2;                                             /* GfTrxn2=GfTrxnCH4*/
                K_equilibriumCH4 = pow(2.71828, (-GfTrxn2 / R));

                for (i = 0; i < 7; i++)
                {
                    K_equilibriumP[i] = 1.0;
                }
                GfTP[0] = GfTCH4;
                for (i = 0; i < 24; i++)
                {
                    activity[i] = P* y[i] * phi[i] / P_ref;  
                    C_UDMI(c, t, i+79) = activity[i];

                }
              /*  printf("the activity[1]: %f", activity[i]);*/
                kf[1] = 0.0000025;                                                                          /*Defining forward rate constant for reaction 1*/
                kf[2] = 0.000225;                                                                           /*Defining forward rate constant for reaction 2*/
                kf[3] = 0.0019379;                                                                          /*Defining forward rate constant for reaction 3*/
                kf[4] = 0.2170533;
                kf[5] = 0.0192921;
                ACH2 = 2 * kf[5];
                BCH2 = (kf[1] / K_equilibriumCH2) * activity[3] + kf[2] * activity[0] + kf[3] * activity[4] + kf[4] * (activity[5] + activity[6] + activity[7] + activity[8] + activity[9] + activity[10] + activity[11]+ activity[12] + activity[13] + activity[14] + activity[15] + activity[16] + activity[16] + activity[17] + activity[18] + activity[19] + activity[20] + activity[21] + activity[22] + activity[23]);
                C = -(kf[1] * activity[1] * activity[0] * activity[0] + (kf[2] / K_equilibriumCH4) * activity[4] + (kf[3] / K_equilibriumP[0]) * activity[5] + (kf[4] / K_equilibriumP[1]) *( activity[6] +  activity[7] + activity[8] +  activity[9] + activity[10] +  activity[11] + activity[12] + activity[13] + activity[14] + activity[15] + activity[16] + activity[16] + activity[17] + activity[18] + activity[19] + activity[20] + activity[21] + activity[22] + activity[23]));
                activityCH2 = (-BCH2 + sqrt(BCH2 * BCH2 - 4 * ACH2 * C)) / (2 * ACH2);
                C_UDMI(c, t, 78) = activityCH2;

                density_catalyst = 700.0;
                porosity_catalyst = 0.44;
                 react1 = activity[0] * activity[0] * activity[1];
                prod1 = activity[3] * activityCH2;
                react2 = activity[0] * activityCH2;
                prod2 = activity[4];
                react3 = activity[4] * activityCH2;
                prod3 = activity[5];
                for (i = 4; i < 22; i++)
                {
                    react[i] = activity[i + 1] * activityCH2;
                    prod[i] = activity[i + 2];
                }

                rate1 = (kf[1] * (react1 - prod1 / K_equilibriumCH2)) * density_catalyst * porosity_catalyst * 0.001*10000;
                C_UDMI(c, t,12 ) = rate1;
              // printf("The rate of monomer formation: %f",C_UDMI(c,t,12));

                 rate2 = (kf[2] * (react2 - prod2 / K_equilibriumCH4)) * density_catalyst * porosity_catalyst * 0.001 * 10000;
                 C_UDMI(c, t, 13) = rate2;

                 rate3 = (kf[3] * (react3 - prod3)) * density_catalyst * porosity_catalyst * 0.001 * 10000;
                 C_UDMI(c, t, 14) = rate3;

                for (i = 4; i < 22; i++)
                {
                    rate[i] = (kf[4] * (react[i] - prod[i])) * density_catalyst * porosity_catalyst * 0.001 * 10000;
                    C_UDMI(c, t, i+137) = rate[i];
                }

            }
             end_c_loop_all(c, t)
        }


        if (DOMAIN_ID(d) == 3)
            thread_loop_c(t, d)
        {
            begin_c_loop_all(c, t)
            {
                for (n = 0; n < 20; n++)
                {
                    MWP[n] = 2.01590803640286 + 14.026579298921 * (n + 1);
                }
                for (n = 0; n < 4; n++)
                {
                    MW[n] = MWsyngas[n];
                }
                for (n = 0; n < 20; n++)
                {
                    MW[n + 4] = MWP[n];
                }

                tot_gas_moles = (C_YI(c, t, 0) / MW[3]) / tot_gas_moles +  C_YI(c, t, 1) / MW[8] + C_YI(c, t, 2) / MW[9] + C_YI(c, t, 3) / MW[10] + C_YI(c, t, 4) / MW[11] + C_YI(c, t, 5) / MW[12] + C_YI(c, t, 6) / MW[13] + C_YI(c, t, 7) / MW[14] + C_YI(c, t, 8) / MW[15] + C_YI(c, t, 9) / MW[16] + C_YI(c, t, 10) / MW[17] + C_YI(c, t, 11) / MW[18] + C_YI(c, t, 12) / MW[19] + C_YI(c, t, 13) / MW[20] + C_YI(c, t, 14) / MW[21] + C_YI(c, t, 15) / MW[22] + C_YI(c, t, 16) / MW[23];
                //    printf("\nweight frac: C_YI H2,CO,CO2,H2O,C1-C8 = %f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f", C_YI(c, t, 0), C_YI(c, t, 1), C_YI(c, t, 2), C_YI(c, t, 3), C_YI(c, t, 4), C_YI(c, t, 5), C_YI(c, t, 6), C_YI(c, t, 7), C_YI(c, t, 8), C_YI(c, t, 9), C_YI(c, t, 10), C_YI(c, t, 11));

                x[0] = 0;                                                                     /*H2*/
                x[1] = 0;                                                                     /*CO*/
                x[2] = 0;                                                                      /*CO2*/
                x[3] = (C_YI(c, t, 0) / MW[3]) / tot_gas_moles;                                /*H2O*/
                x[4] = 0;                                                                    /*CH4*/
                x[5] = 0;                                                                      /*C2*/
                x[6] = 0;                                                                     /*C3*/
                x[7] = 0;                                                                      /*C4*/
                x[8] = (C_YI(c, t, 1) / MW[8]) / tot_gas_moles;                            /*C5*/
                x[9] = (C_YI(c, t, 2) / MW[9]) / tot_gas_moles;                            /*C6*/
                x[10] = (C_YI(c, t, 3) / MW[10]) / tot_gas_moles;                         /*C7*/
                x[11] = (C_YI(c, t, 4) / MW[11]) / tot_gas_moles;                          /*C8*/
                x[12] = (C_YI(c, t, 5) / MW[12]) / tot_gas_moles;                            /*H2*/
                x[13] = (C_YI(c, t, 6) / MW[13]) / tot_gas_moles;                            /*CO*/
                x[14] = (C_YI(c, t, 7) / MW[14]) / tot_gas_moles;                            /*CO2*/
                x[15] = (C_YI(c, t, 8) / MW[15]) / tot_gas_moles;                            /*H2O*/
                x[16] = (C_YI(c, t, 9) / MW[16]) / tot_gas_moles;                            /*CH4*/
                x[17] = (C_YI(c, t, 10) / MW[17]) / tot_gas_moles;                            /*C2*/
                x[18] = (C_YI(c, t, 11) / MW[18]) / tot_gas_moles;                            /*C3*/
                x[19] = (C_YI(c, t, 12) / MW[19]) / tot_gas_moles;                            /*C4*/
                x[20] = (C_YI(c, t, 13) / MW[20]) / tot_gas_moles;                            /*C5*/
                x[21] = (C_YI(c, t, 14) / MW[21]) / tot_gas_moles;                            /*C6*/
                x[22] = (C_YI(c, t, 15) / MW[22]) / tot_gas_moles;                         /*C7*/
                x[23] = (C_YI(c, t, 16) / MW[23]) / tot_gas_moles;
                //  printf("\ntot_moles & mole frac: yH2,CO,CO2,H2O,C1-C8 = %f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f", tot_gas_moles,y[0], y[1], y[2], y[3], y[4], y[5], y[6], y[7], y[8], y[9], y[10], y[11]);
               
                C_UDMI(c, t, 272) = C_VOF(c, t);
                for (n = 0; n < 24; n++)
                {
                    C_UDMI(c, t, 245 + n) = x[n];
                   bi[n] = C_UDMI(c, t, 160 + n);
                   ai[n] = C_UDMI(c, t, 185 + n);
                }
                asum = 0.0;
                bsum = 0.0;
                for (n = 0; n < 24; n++)
                {
                    for (p = 0; p < 24; p++)
                    {
                        m[n][p] = pow((ai[n] * ai[p]), 0.5);

                        asum = asum + x[n] * x[p] * m[n][p];
                    }
                    bsum = bsum + x[n] * bi[n];
                }
               
                a = asum; b = bsum;
                T = 493.15;
                P = 2000000;
                A = a * P / (R * R * T * T);
                B = b * P / (R * T);
                //  printf("\n the Pressure = %f and Temperature = %f", P, T);
                
                n_iterations = N_ITER;
                if (n_iterations < 800)
                  

                for (i = 0; i < 24; i++)
                {
                    long_term[i][0] = 0.0;
                }

                for (i = 0; i < 24; i++)
                {
                    for (j = 0; j < 24; j++)
                    {
                        long_term[i][j + 1] = long_term[i][j] + m[i][j] * x[j];
                    }
                }

                Zl = C_UDMI(c, t, 23);
                //  Zv = 0.97;

                for (i = 0; i < 24; i++)
                {
                    q[i] = (A / B) * (2 * long_term[i][12] / A - bi[i] / b);
                    lnphi[i] = (bi[i] / b) * (Zl - 1) - log(Zl - B) - q[i] * I;
                    phi[i] = exp(lnphi[i]);
                    //   C_UDMI(c, t, i + 51) = q[i];
                    C_UDMI(c, t, i + 220) = phi[i];
                }

            }
            end_c_loop_all(c, t)
        }

    
    }

}


DEFINE_INIT(MY_INIT, md)
{
    int pt;
    cell_t c;
    Thread* t;
    Domain* d;
    sub_domain_loop(d, md, pt)
    {
        if (DOMAIN_ID(d) == 2)
            thread_loop_c(t, d)
        {
            begin_c_loop_all(c, t)
            {
                n = 0; s = 0; p = 0; i = 0; j = 0;
                alphax = 0.1;  Mol_Weight = 1e-6; wH2 = -0.215993; wCO = 0.045; wCO2 = 0.223621; wH2O = 0.344; react3 = 1e-6; prod3 = 1e-6;
                a = 1e-6; b = 1e-6; P = 1e-6; T = 1e-6; A = 1e-6; B = 1e-6; V_gas = 1e-6; V_liq = 1e-6; Zv = 1e-6; Zl = 1e-6; rho_gas1 = 1e-6; rho_liquid = 1e-6; rhol = 1e-6; asum = 0.0; bsum = 0.0; change = 1e-6; dfdz = 1e-6; tot_gas_moles = 1e-6;
                P_ref = 1e-6; T_ref = 1e-6; aCO = 1e-6; bCO = 1e-6; cCO = 1e-6; dCO = 1e-6; Hf298CO = 1e-6; Gf298CO = 1e-6; GammaCO = 1e-6; GfTCO = 1e-6; aH2 = 1e-6; bH2 = 1e-6; cH2 = 1e-6; dH2 = 1e-6; Hf298H2 = 1e-6; Gf298H2 = 1e-6;  GammaH2 = 1e-6; GfTH2 = 1e-6; aH2O = 1e-6; bH2O = 1e-6; cH2O = 1e-6; dH2O = 1e-6; Hf298H2O = 1e-6; Gf298H2O = 1e-6; GammaH2O = 1e-6;
                GfTrxn1 = 1e-6; GfTrxn2 = 1e-6; rate1 = 1e-6; rate2 = 1e-6; GfTH2O = 1e-6; aCH4 = 1e-6; bCH4 = 1e-6; cCH4 = 1e-6; dCH4 = 1e-6; Hf298CH4 = 1e-6; Gf298CH4 = 1e-6; GammaCH4 = 1e-6; GfTCH4 = 1e-6; aC2H6 = 1e-6; bC2H6 = 1e-6; cC2H6 = 1e-6; dC2H6 = 1e-6; GfTCH2 = 1e-6;  rate3 = 1e-6;
                prod1 = 1e-6; react1 = 1e-6; prod2 = 1e-6; react2 = 1e-6; K_equilibriumCH2 = 1e-6; K_equilibriumCH4 = 1e-6; K_equilibriumC2H6 = 1e-6; activityCH2 = 1e-6; density_catalyst = 1e-6; porosity_catalyst = 1e-6; ACH2 = 1e-6; BCH2 = 1e-6; C = 1e-6;
                X = 0.0; omegaA = 1e-6; omegaB = 1e-6; I = 0.0; PcH2 = 1313000.0; PcCO = 3494000.0; PcCO2 = 7383000.0; PcH2O = 22064000.0; TcH2 = 33.19; TcCO = 132.85;  TcCO2 = 304.21; TcH2O = 647.14;

                MWsyngas[0] = 2.01588;
                MWsyngas[1] = 28.0101;
                MWsyngas[2] = 44.0095;
                MWsyngas[3] = 18.0153;

                for (i = 0; i < 13; i++)
                {
                    MW[i] = 0.0;
                }
               // C_UDMI(c, t, 56) = B;
                C_UDMI(c, t, 57) = 0.0;
                C_UDMI(c, t, 58) = 0.0;

                for (i = 0; i < 12; i++)
                {
                    w[i] = 0.0;
                    fugacitycoeff[i] = 0.0;
                    ac[i] = 0.0;
                    k[i] = 0.0;
                    r[i] = 0.0;
                    Tr[i] = 0.0;
                    Pr[i] = 0.0;
                    D[i] = 0.0;
                    y[i] = 0.0; activity[i] = 0.0; alpha0[i] = 0.0; alpha1[i] = 0.0; alphaT[i] = 0.0; ai[i] = 0.0; bi[i] = 0.0; q[i] = 0.0; lnphi[i] = 0.0; phi[i] = 0.0; Tc[i] = 0.0; Pc[i] = 0.0;
                }
                for (i = 0; i < 8; i++)
                {
                    wP[i] = 0.0;
                    rate[i] = 0.0;
                }
                for (i = 0; i < 7; i++)
                {
                    GfTrxnP[i] = 0.0;
                }
                for (i = 0; i < 8; i++)
                {
                    MWP[i] = 0.0;
                    PcP[i] = 0.0;
                    TcP[i] = 0.0;
                }
                for (i = 0; i < 6; i++)
                {
                    GammaP[i] = 0.0;
                    GftP[i] = 0.0;
                    Hft298P[i] = 0.0;
                    GftRxns[i] = 0.0;
                    K_equilibriums[i] = 0.0;
                    Hf298P[i] = 0.0;
                    Gf298P[i] = 0.0;
                    GfTP[i] = 0.0;
                    prod[i] = 0.0;
                    react[i] = 0.0;
                    kf[i] = 0.0;
                    K_equilibriumP[i] = 0.0;
                }

                for (i = 0; i < 200; i++)
                {
                    f[i] = 0.0;
                    z[i] = 0.0;
                }
                for (i = 0; i < 12; i++)
                {
                    for (j = 0; j < 13; j++)
                    {
                        long_term[i][j] = 0.0;
                    }
                }
                for (i = 0; i < 12; i++)
                {
                    for (j = 0; j < 12; j++)
                    {
                        m[i][j] = 0.0;
                        Q[i][j] = 0.0;
                    }
                }

            }
            end_c_loop_all(c, t)
        }
        if (DOMAIN_ID(d) == 5)
            thread_loop_c(t, d)
        {
            begin_c_loop_all(c, t)
            {
                C_UDMI(c, t, 26) = 0.00001;
            }
            end_c_loop_all(c, t)
        }
    }
}


