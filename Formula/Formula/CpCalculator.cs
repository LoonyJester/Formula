using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace Formula
{
    public class CpCalculator
    {

        public double CalculatedFCR(double inletRate, double exchangeVolume, double exchangeVolumeChange,
    double exchangeHct, double exchangeFcr, double exchangeEndHct, double exchangeTargetFcr,
    double exchangeRfHct, int flag)
        {

            double therapeuticFcr = exchangeTargetFcr / exchangeFcr;
            double m = 0;
            double p;

            if (Math.Abs(exchangeVolumeChange) > 0.0000001)
            {


                m = 1 - Math.Log(therapeuticFcr) / (Math.Log(1 + (exchangeVolumeChange / exchangeVolume)));
                p = Math.Pow(therapeuticFcr, m / (m - 1));
            }
            else
            {
                p = therapeuticFcr;
            }

            var f = exchangeRfHct / exchangeHct * ((p - 1) / (p - (exchangeEndHct / exchangeHct)));
            double qrf;

            if (Math.Abs(exchangeVolumeChange) > 0.0000001)
            {

                qrf = inletRate * (m / (f * (m - 1)));
            }
            else
            {
                qrf = inletRate / f;
            }

            double qp = (f - 1) * qrf;
            double exchangeTime;

            if (Math.Abs(exchangeVolumeChange) > 0.0000001)
            {
                double flowRateDifference = qrf + qp - inletRate;
                m = (qrf + qp) / flowRateDifference;
                exchangeTime = exchangeVolume / flowRateDifference *
                                ((1 / Math.Pow(therapeuticFcr, (1 / (m - 1)))) - 1);
            }
            else
            {
                exchangeTime = -exchangeVolume / inletRate * Math.Log(therapeuticFcr);
            }

            if (flag == 1)
                return qrf;
            if (flag == 2)
                return qp;
            if (flag == 3)
                return exchangeTime;
            if (flag == 4)
                return therapeuticFcr;
            if (flag == 5)
                return p;
            if (flag == 6)
                return f;

            return 0;

        }

        public double CalculatedFCRVR(double exchangeHct, double exchangeRfHct, double exchangeEndHct,
double exchangeVolumeChange, double exchangeVolume, double exchangeReplacementVolume, double inletRate,
int flag)
        {
            if (flag != 1) return 0;

            double a = exchangeVolumeChange;
            double b = exchangeVolume;
            double c = exchangeReplacementVolume;
            double d = exchangeRfHct;
            double e = exchangeEndHct;
            double g = exchangeHct;

            double q1;
            double q5;
            double t;
            double fcrCalc = 0;
            double dq;
            double v;
            double m;

            if (g == e)
            {


                if (a == 0)
                {

                    q1 = (g / d) * inletRate;
                    q5 = inletRate - q1;
                    t = c / q1;
                    fcrCalc = Math.Exp(-inletRate * t / b);
                }


                if (a != 0)
                {

                    q1 = (g / d) * (1 / (1 - ((g / d) * (a / c)))) * inletRate;
                    q5 = ((a / c) - 1) * q1 + inletRate;
                    t = c / q1;
                    dq = q1 + q5 - inletRate;
                    m = (q1 + q5) / dq;
                    v = b + dq * t;
                    fcrCalc = Math.Pow((b / v), (m - 1));
                }

            }

            if (g != e)
            {



                if (a == 0)
                {
                    exchangeVolumeChange = 1;
                    a = 1;
                }
                var x = b / (a + b);


                var F0 = d / e;

                double n = 0;
                double F = 0;
                do
                {
                    F = F0;
                    double gF = Math.Pow((1 / (1 + a / b)), F * c / a) - (d - F * e) / (d - F * g);
                    double dgF = ((c * Math.Log(x) * Math.Pow(d - F * g, 2) * Math.Pow(x, c * F / a)) + (a * d * (e - g))) /
                                 (a * Math.Pow((d - (F * g)), 2));
                    F0 = F - (gF / dgF);
                    n = n + 1;
                } while (!(Math.Abs(F0 - F) < 0.0001 || n > 20));



                F = F0;

                q1 = (1 / (F - (a / c))) * inletRate;
                q5 = (F - 1) * q1;
                t = c / q1;
                dq = q1 + q5 - inletRate;
                m = (q1 + q5) / dq;
                v = b + dq * t;
                fcrCalc = Math.Pow(b / v, m - 1);

            }
            return fcrCalc;
        }


        public double[] RBCXTT(double inletRate, double qrf, double qp, double currentRfHct, double modeTime,
double initialModeVolume, double initialModeHct, double initialModeFcr)
        {

            double flowRateDifference = qrf + qp - inletRate;
            double m = 0;
            if (Math.Abs(flowRateDifference) > 0.0000001)
            {
                m = (qrf + qp) / flowRateDifference;
            }

            double f = (qrf + qp) / qrf;

            double currentPatientVolume = initialModeVolume + (flowRateDifference * modeTime);
            double currentFcr;

            if (Math.Abs(flowRateDifference) > 0.0000001)
            {


                currentFcr = initialModeFcr * Math.Pow((initialModeVolume / currentPatientVolume), (m - 1));
            }
            else
            {
                currentFcr = initialModeFcr * Math.Exp(-inletRate * modeTime / initialModeVolume);
            }
            double currentPatientHct;

            if (Math.Abs(flowRateDifference) > 0.0000001)
            {


                currentPatientHct = (currentRfHct / f) - (((currentRfHct / f) - initialModeHct) *
                                                            Math.Pow((initialModeVolume / currentPatientVolume), m));
            }
            else
            {
                currentPatientHct = (currentRfHct / f) -
                                      (((currentRfHct / f) - initialModeHct) *
                                       Math.Exp(-inletRate * modeTime / initialModeVolume));
            }

            return new[] { currentPatientVolume, currentPatientHct, currentFcr };

        }

        #region Enums
        public enum ProcedureType
        {
            Depletion,
            Exchange,
            DepletionExchange
        }

        public enum GenderType
        {
            Mail,
            Female
        }

        public enum WeightType
        {
            lbs,
            kg
        }

        public enum HeightType
        {
            inc,
            sm
        }

        public enum BloodVolumeType
        {
            Nadler,
            Override
        }

        public enum ReplacementVolumeType
        {
            Calculate,
            Override
        }
        #endregion Enums


        public double Min(double a, double b)
        {
            double functionReturnValue = 0;
            if (a < b)
            {
                functionReturnValue = a;
            }
            else
            {
                functionReturnValue = b;
            }
            return functionReturnValue;
        }


        public double Max(double a, double b)
        {
            double functionReturnValue = 0;
            if (a > b)
            {
                functionReturnValue = a;
            }
            else
            {
                functionReturnValue = b;
            }
            return functionReturnValue;
        }

        public object EP(ProcedureType Procedure_Type, GenderType Gender, double Weight, WeightType Weight_Units, double Height, 
            HeightType Height_Units, double Total_Blood_Volume, BloodVolumeType Total_Blood_Volume_Override, double FCR, double Initial_Hct,
      double Target_End_Hct, double Avg_RF_Hct, double Replacement_Volume, ReplacementVolumeType Replacement_Volume_Override, 
      double Target_Depletion_Hct, double Blood_Warmer_Volume, double Flag = 0)
        {
            object functionReturnValue = null;

           string Predicted_Error = "None";

            double IP_WB_Volume = 0;
            double IP_Replaced_Volume = 0;
            double PD_Replaced_Volume = 0;
            double Total_Waste_Processed = 0;
            double RF_Volume = 0;
            double FCR_Factor = 1.085;
            double RBC_Citrate_Concentration = 21.4;
            double Percent_AC_in_RBC = 2;

            if (Height_Units == HeightType.inc)
            {
                Height = Height / 0.3937008 / 100;
            }
            else
            {
                Height = Height / 100;
            }

            if (Weight_Units == WeightType.lbs)
            {
                Weight = Weight / 2.2046223;
            }

            Initial_Hct = Initial_Hct / 100;
            Avg_RF_Hct = Avg_RF_Hct / 100;
           double Blood_Prime_Hct = 0;
            double Estimation_Hct = 0;// Estimation_Hct Estimation_Hct = Estimation_Hct / 100;
            Target_Depletion_Hct = Target_Depletion_Hct / 100;
            Target_End_Hct = Target_End_Hct / 100;
            FCR = FCR / 100;

            if (Total_Blood_Volume_Override == BloodVolumeType.Nadler)
            {
               // double Estimation_Hct;
                  double Low_TBV_Factor = 80;
                if (Weight < 25)
                {
                    Total_Blood_Volume = Weight * Low_TBV_Factor;
                }
                else
                {
                    if (Gender == GenderType.Female)
                    {
                        Total_Blood_Volume = (0.3561 * Math.Pow(Height, 3) + 0.03308 * Weight + 0.1833) * 1000;
                    }
                    else
                    {
                        Total_Blood_Volume = (0.3669 * Math.Pow(Height, 3) + 0.03219 * Weight + 0.6041) * 1000;
                    }
                }
            }

            double Fluid_Balance = 0;
            double Total_Reinfusion_Volume = 0;


            double IP_Return_Path_Volume = 60 - 15 + Blood_Warmer_Volume;
            double PD_Return_Path_Volume = 60 + Blood_Warmer_Volume;

            double Calculated_CIR_Limit = 1.25 * Weight;

            double Cit_Con_AC = 21.4;

            double Entered_Patient_Volume = Total_Blood_Volume;
            double Entered_Patient_Hct = Initial_Hct;
            double Entered_Avg_RF_Hct = Avg_RF_Hct;
            double Entered_Fluid_Balance = Fluid_Balance;
            double Entered_Max_WB_Rate = 50;
            double Entered_AC_Ratio = 12;
            double Entered_CIR = 1.25;
          string  Entered_AC_Type = "ACD";
        string    Entered_Blood_Prime = "None";
            double Entered_End_Hct = Target_End_Hct;
string            Entered_Divert_Prime = "Yes";
      string      Entered_Reinfusion = "No";
            double Entered_FCR = Math.Exp(Math.Log(FCR) / FCR_Factor) * (Entered_End_Hct / Entered_Patient_Hct) * ((Entered_Patient_Volume + Entered_Fluid_Balance - Total_Reinfusion_Volume) / Entered_Patient_Volume);

            double Entered_Depletion_Hct = 0;
            if (Procedure_Type == ProcedureType.Depletion)
            {
                Entered_Depletion_Hct = Entered_End_Hct;
            }
            else
            {
                Entered_Depletion_Hct = Target_Depletion_Hct;
            }
            double Entered_Replacement_Volume =0;
            if (Replacement_Volume_Override == ReplacementVolumeType.Override)
            {
                Entered_Replacement_Volume = Replacement_Volume;
            }

            Predicted_Error = "None";
            double Predicted_Replacement_Volume = 0;
            double Predicted_Other_RF_Volume = 0;

            double RF_AC_Fraction = Percent_AC_in_RBC / 100;
            double Cit_Con_Rep_Fluid = RF_AC_Fraction * RBC_Citrate_Concentration;
            double CIR_Limited_Qrf = Calculated_CIR_Limit / Cit_Con_Rep_Fluid;

            Cit_Con_AC = 21.4;

            double EqQb = Calculated_CIR_Limit * (Entered_AC_Ratio + 1) / Cit_Con_AC;

            string Divert_Prime = Entered_Divert_Prime;
            double Max_WB_Rate = Entered_Max_WB_Rate;
            double AC_Ratio = Entered_AC_Ratio;
            double Estimation_Volume = Entered_Patient_Volume;
            Estimation_Hct = Entered_Patient_Hct;
            double Estimation_FCR = 1;
            double RF_Used = 0;
            Avg_RF_Hct = Entered_Avg_RF_Hct;
            Target_End_Hct = Entered_End_Hct;
            double Target_FCR = Entered_FCR;
            double Target_Volume = Entered_Patient_Volume + Entered_Fluid_Balance - Total_Reinfusion_Volume;

            if (Target_Volume < 0.9 * Entered_Patient_Volume)
            {
                Predicted_Error = "Increase Target Fluid Blance";
            }

            Target_Depletion_Hct = Entered_Depletion_Hct;
            double Target_Replacement_Volume = Entered_Replacement_Volume;
            double CIR_Limit = Calculated_CIR_Limit;

            if (Procedure_Type == ProcedureType.Exchange)
            {
                double Qrf_IP = 0;
                double Qwb_IP = Min(CIR_Limited_Qrf * (AC_Ratio + 1) / AC_Ratio, Max_WB_Rate);
                Qwb_IP = Min(Qwb_IP, 50);
                double Qac_IP = Qwb_IP / (AC_Ratio + 1);
                if (Qac_IP < 1.3)
                {
                    Predicted_Error = "AC Ratio Unachievable";
                }
                if (Estimation_Hct > Avg_RF_Hct)
                {
                    Qrf_IP = Qwb_IP - Qac_IP;
                }
                else
                {
                    Qrf_IP = (Estimation_Hct / Avg_RF_Hct) * (Qwb_IP - Qac_IP);
                }
                double Qp_IP = Qwb_IP - Qac_IP - Qrf_IP;
                double Remaining_IP_WB_Volume = 105 - IP_WB_Volume;
                double Est_Mode_Time = Remaining_IP_WB_Volume / Qwb_IP;
                double Total_Return_Volume = (Qrf_IP + Qp_IP) * Est_Mode_Time;
                double Mode_RF_Volume = Qrf_IP * Est_Mode_Time;
                double Mode_Other_RF_Volume = 0;
                if (IP_Replaced_Volume < IP_Return_Path_Volume)
                {
                    IP_Return_Path_Volume = IP_Return_Path_Volume - IP_Replaced_Volume;
                    if (IP_Return_Path_Volume > Total_Return_Volume)
                    {
                        IP_Return_Path_Volume = Total_Return_Volume;
                    }
                    double Partial_Time = IP_Return_Path_Volume / (Qwb_IP - Qac_IP);

                    double[] Estimation = RBCXTT(Qwb_IP - Qac_IP, Qrf_IP, Qp_IP, 0, Partial_Time, Estimation_Volume, Estimation_Hct, Estimation_FCR);

                    Estimation_Volume = Estimation[0];
                    Estimation_Hct = Estimation[1];
                    Estimation_FCR = Estimation[2];
                }
                else
                {
                    IP_Return_Path_Volume = 0;
                }

                // TODO: it line can be wrong 
                // oroginal kine is If Return_Path_Volume < Total_Return_Volume Then
                // but there is not any 
                if (IP_Return_Path_Volume < Total_Return_Volume)
                {
                    double Partial_Time = (Total_Return_Volume - IP_Return_Path_Volume) / (Qwb_IP - Qac_IP);

                    double[] Estimation = RBCXTT(Qwb_IP - Qac_IP, Qrf_IP, Qp_IP, Avg_RF_Hct, Partial_Time, Estimation_Volume, Estimation_Hct, Estimation_FCR);

                    Estimation_Volume = Estimation[0];
                    Estimation_Hct = Estimation[1];
                    Estimation_FCR = Estimation[2];
                }

                Predicted_Replacement_Volume = Mode_RF_Volume;
                Predicted_Other_RF_Volume = Mode_Other_RF_Volume;

                double Qwb_Ex = Max_WB_Rate;
                double Qac_Ex = Qwb_Ex / (AC_Ratio + 1);
                if (Qac_Ex < 1.3)
                {
                    Predicted_Error = "AC Ratio Unachievable";
                }
                if (Replacement_Volume_Override == ReplacementVolumeType.Override)
                {
                    if (Target_Replacement_Volume - RF_Used - Predicted_Replacement_Volume <= 0)
                    {
                        Predicted_Error = "Increase Replacement Volume";
                    }
                    var Calculated_FCR = CalculatedFCRVR(Estimation_Hct, Avg_RF_Hct, Target_End_Hct, Target_Volume - Estimation_Volume, Estimation_Volume, Target_Replacement_Volume - RF_Used - Predicted_Replacement_Volume, (Qwb_Ex - Qac_Ex), 1);
                    Target_FCR = Calculated_FCR * Estimation_FCR;
                    if (Target_FCR <= 0)
                    {
                        Predicted_Error = "Decrease Replacement Volume";
                    }
                    if (Target_FCR >= Estimation_FCR)
                    {
                        Predicted_Error = "Decrease Target End Hct";
                    }

                }
                else
                {

                    if (Target_FCR <= 0)
                    {
                        Predicted_Error = "Increase Target FCR";
                    }
                    if (Target_FCR >= Estimation_FCR)
                    {
                        Predicted_Error = "Decrease Target FCR";
                    }
                }

                if ((Target_FCR / Estimation_FCR) >= (0.99 + ((Target_Volume - Estimation_Volume) / Estimation_Volume)))
                {
                    Predicted_Error = "Increase Target Fluid Balance";
                }

                double Qrf_Ex = CalculatedFCR((Qwb_Ex - Qac_Ex), Estimation_Volume, (Target_Volume - Estimation_Volume), Estimation_Hct, Estimation_FCR, Target_End_Hct, Target_FCR, Avg_RF_Hct, 1);
                Qp_Ex = CalculatedFCR((Qwb_Ex - Qac_Ex), Estimation_Volume, (Target_Volume - Estimation_Volume), Estimation_Hct, Estimation_FCR, Target_End_Hct, Target_FCR, Avg_RF_Hct, 2);
                Est_Mode_Time = CalculatedFCR((Qwb_Ex - Qac_Ex), Estimation_Volume, (Target_Volume - Estimation_Volume), Estimation_Hct, Estimation_FCR, Target_End_Hct, Target_FCR, Avg_RF_Hct, 3);
                Max_Hct = max(Estimation_Hct, Target_End_Hct);
                Plasma_AC_Fraction = 1 / (((1 - Max_Hct) * AC_Ratio) + 1);
                Estimated_CIR = (Qrf_Ex * Cit_Con_Rep_Fluid) + (Qp_Ex * Cit_Con_AC * Plasma_AC_Fraction);

                if (Estimated_CIR > CIR_Limit)
                {
                    CIR_Adjustment = CIR_Limit / Estimated_CIR;
                    Qwb_Ex = Qwb_Ex * CIR_Adjustment;
                    Qac_Ex = Qac_Ex * CIR_Adjustment;
                    if (Qac_Ex < 1.3)
                    {
                        Predicted_Error = "AC Ratio Unachievable";
                    }
                    Qrf_Ex = CalculatedFCR((Qwb_Ex - Qac_Ex), Estimation_Volume, (Target_Volume - Estimation_Volume), Estimation_Hct, Estimation_FCR, Target_End_Hct, Target_FCR, Avg_RF_Hct, 1);
                    Qp_Ex = CalculatedFCR((Qwb_Ex - Qac_Ex), Estimation_Volume, (Target_Volume - Estimation_Volume), Estimation_Hct, Estimation_FCR, Target_End_Hct, Target_FCR, Avg_RF_Hct, 2);
                    Est_Mode_Time = CalculatedFCR((Qwb_Ex - Qac_Ex), Estimation_Volume, (Target_Volume - Estimation_Volume), Estimation_Hct, Estimation_FCR, Target_End_Hct, Target_FCR, Avg_RF_Hct, 3);
                }

                prbc_Hct = 0.8;
                AC_Dilution = AC_Ratio / (AC_Ratio + 1);
                Max_Hct = max(Estimation_Hct, Target_End_Hct);
                Hct_ACWB = AC_Dilution * Max_Hct;
                Qp_max = Qwb_Ex * (1 - (Hct_ACWB / prbc_Hct));

                if (Qp_Ex > Qp_max)
                {
                    Predicted_Error = "Unable to meet Hct";
                }

                if (Qp_Ex < 1.3)
                {
                    if (Qp_Ex < 0 | Qp_Ex > 1E-07)
                    {
                        Predicted_Error = "Decrease Target End Hct";
                    }
                    else
                    {
                        Qp_Ex = 0;
                    }
                }

                if (Qrf_Ex > 150)
                {
                    Predicted_Error = "Decrease Target Fluid Balance";
                }

                if (Qrf_Ex < 1.3)
                {
                    Predicted_Error = "Decrease Target FCR";
                }

                Mode_RF_Volume = Qrf_Ex * Est_Mode_Time;
                Mode_Other_RF_Volume = 0;

                if (Estimation_Hct > Target_End_Hct)
                {
                    Plasma_AC_Fraction_Start = 1 / (((1 - Estimation_Hct) * AC_Ratio) + 1);
                    Plasma_AC_Fraction_End = 1 / (((1 - Target_End_Hct) * AC_Ratio) + 1);
                    Plasma_AC_Fraction = (Plasma_AC_Fraction_Start + Plasma_AC_Fraction_End) / 2;
                }
                else
                {
                    Plasma_AC_Fraction = 1 / (((1 - Target_End_Hct) * AC_Ratio) + 1);
                }

                Mode_Time = Est_Mode_Time;

                Estimation = RBCXTT(Qwb_Ex - Qac_Ex, Qrf_Ex, Qp_Ex, Avg_RF_Hct, Est_Mode_Time, Estimation_Volume, Estimation_Hct, Estimation_FCR);

                Current_Patient_Volume = Estimation(0);
                Current_Patient_Hct = Estimation(1);
                Current_FCR = Estimation(2);

                Estimation_Volume = Current_Patient_Volume;
                Estimation_Hct = Current_Patient_Hct;
                Estimation_FCR = Current_FCR;

                Predicted_Replacement_Volume = Predicted_Replacement_Volume + Mode_RF_Volume;
                Predicted_Other_RF_Volume = Predicted_Other_RF_Volume + Mode_Other_RF_Volume;


            }
            else
            {
                Qwb_IP = min(Max_WB_Rate, 50);

                Qac_IP = Qwb_IP / (AC_Ratio + 1);
                if (Qac_IP < 1.3)
                {
                    Predicted_Error = "AC Ratio Unachievable";
                }
                Qrf_IP_Est = Qwb_IP - Qac_IP;

                Remaining_IP_WB_Volume = 105 - IP_WB_Volume;
                Est_Mode_Time = Remaining_IP_WB_Volume / Qwb_IP;

                Mode_RF_Volume = 0;
                Mode_Other_RF_Volume = (Qwb_IP - Qac_IP) * Est_Mode_Time;

                Estimation = RBCXTT(Qwb_IP - Qac_IP, Qwb_IP - Qac_IP, 0, 0, Est_Mode_Time, Estimation_Volume, Estimation_Hct, Estimation_FCR);

                Current_Patient_Volume = Estimation(0);
                Current_Patient_Hct = Estimation(1);
                Current_FCR = Estimation(2);

                Estimation_Volume = Current_Patient_Volume;
                Estimation_Hct = Current_Patient_Hct;
                Estimation_FCR = Current_FCR;

                Predicted_Replacement_Volume = Mode_RF_Volume;
                Predicted_Other_RF_Volume = Mode_Other_RF_Volume;

                Start_Depletion_Hct = Estimation_Hct;
                Volume_to_Remove = -Estimation_Volume * Log(Target_Depletion_Hct / Estimation_Hct);

                Qwb_Dep = min(Max_WB_Rate, EqQb);

                Qac_Dep = Qwb_Dep / (AC_Ratio + 1);

                if (Qac_Dep < 1.3)
                {
                    Predicted_Error = "AC Ratio Unachievable";
                }

                prbc_Hct = 0.8;
                AC_Dilution = AC_Ratio / (AC_Ratio + 1);
                Avg_Depletion_Hct = (Estimation_Hct + Target_Depletion_Hct) / 2;
                Avg_Hct_ACWB = AC_Dilution * Avg_Depletion_Hct;
                Qp_Ideal_Est = Qwb_Dep * (1 - (Avg_Hct_ACWB / prbc_Hct));

                Qrf_Dep_Est = Qwb_Dep - Qac_Dep - Qp_Ideal_Est;

                if (Procedure_Type == "Depletion/Exchange")
                {
                    Volume_to_Remove = Volume_to_Remove - PD_Return_Path_Volume;
                }

                if (Volume_to_Remove < 0)
                {
                    Predicted_Error = "Depletion Target Unachievable";
                }

                Est_Mode_Time = Volume_to_Remove / (Qwb_Dep - Qac_Dep);
                Est_Depletion_Time = Est_Mode_Time;

                Mode_RF_Volume = 0;
                Mode_Other_RF_Volume = Qrf_Dep_Est * Est_Mode_Time;

                Estimation = RBCXTT(Qwb_Dep - Qac_Dep, Qwb_Dep - Qac_Dep, 0, 0, Est_Mode_Time, Estimation_Volume, Estimation_Hct, Estimation_FCR);

                Current_Patient_Volume = Estimation(0);
                Current_Patient_Hct = Estimation(1);
                Current_FCR = Estimation(2);

                Estimation_Volume = Current_Patient_Volume;
                Estimation_Hct = Current_Patient_Hct;
                Estimation_FCR = Current_FCR;

                Predicted_Replacement_Volume = Predicted_Replacement_Volume + Mode_RF_Volume;
                Predicted_Other_RF_Volume = Predicted_Other_RF_Volume + Mode_Other_RF_Volume;

                Plasma_AC_Fraction_Start = 1 / (((1 - Start_Depletion_Hct) * AC_Ratio) + 1);
                Plasma_AC_Fraction_End = 1 / (((1 - Estimation_Hct) * AC_Ratio) + 1);
                Plasma_AC_Fraction = (Plasma_AC_Fraction_Start + Plasma_AC_Fraction_End) / 2;


                if (Procedure_Type == "Depletion/Exchange")
                {
                    Qwb_PD = Max_WB_Rate;

                    Qac_PD = Qwb_PD / (AC_Ratio + 1);
                    if (Qac_PD < 1.3)
                    {
                        Predicted_Error = "AC Ratio Unachievable";
                    }

                    if (Replacement_Volume_Override == "Override")
                    {
                        if (Target_Replacement_Volume - RF_Used - Predicted_Replacement_Volume <= 0)
                        {
                            Predicted_Error = "Increase_Replacement_Volume";
                        }
                        Calculated_FCR = CalculatedFCRVR(Estimation_Hct, Avg_RF_Hct, Target_End_Hct, Target_Volume - Estimation_Volume, Estimation_Volume, Target_Replacement_Volume - RF_Used - Predicted_Replacement_Volume, Qwb_PD - Qac_PD, 1);
                        PD_FCR = Calculated_FCR * Estimation_FCR;

                        if (PD_FCR <= 0)
                        {
                            Predicted_Error = "Decrease Replacement Volume";
                        }
                        if (PD_FCR >= Estimation_FCR)
                        {
                            Predicted_Error = "Decrease Target End Hct";
                        }
                    }
                    else
                    {
                        PD_FCR = Target_FCR;

                        if (PD_FCR <= 0)
                        {
                            Predicted_Error = "Increase Target FCR";
                        }
                        if (PD_FCR >= Estimation_FCR)
                        {
                            Predicted_Error = "Decrease Target FCR";
                        }
                    }

                    if ((PD_FCR / Estimation_FCR) >= (0.99 + ((Target_Volume - Estimation_Volume) / Estimation_Volume)))
                    {
                        Predicted_Error = "Increase_Target_Fluid_Balance";
                    }

                    Qrf_PD = CalculatedFCR(Qwb_PD - Qac_PD, Estimation_Volume, Target_Volume - Estimation_Volume, Estimation_Hct, Estimation_FCR, Target_End_Hct, PD_FCR, Avg_RF_Hct, 1);
                    Qp_PD = CalculatedFCR(Qwb_PD - Qac_PD, Estimation_Volume, Target_Volume - Estimation_Volume, Estimation_Hct, Estimation_FCR, Target_End_Hct, PD_FCR, Avg_RF_Hct, 2);

                    Max_Hct = max(Estimation_Hct, Target_End_Hct);
                    Plasma_AC_Fraction = 1 / (((1 - Max_Hct) * AC_Ratio) + 1);
                    Estimated_CIR = (Qrf_PD * Cit_Con_Rep_Fluid) + (Qp_PD * Cit_Con_AC * Plasma_AC_Fraction);
                    if (Estimated_CIR > CIR_Limit)
                    {
                        CIR_Adjustment = CIR_Limit / Estimated_CIR;
                        Qwb_PD = Qwb_PD * CIR_Adjustment;
                        Qac_PD = Qac_PD * CIR_Adjustment;
                        if (Qac_PD < 1.3)
                        {
                            Predicted_Error = "AC Ratio Unachievable";
                        }
                        Qrf_PD = CalculatedFCR(Qwb_PD - Qac_PD, Estimation_Volume, Target_Volume - Estimation_Volume, Estimation_Hct, Estimation_FCR, Target_End_Hct, PD_FCR, Avg_RF_Hct, 1);
                        Qp_PD = CalculatedFCR(Qwb_PD - Qac_PD, Estimation_Volume, Target_Volume - Estimation_Volume, Estimation_Hct, Estimation_FCR, Target_End_Hct, PD_FCR, Avg_RF_Hct, 2);
                    }

                    prbc_Hct = 0.8;
                    AC_Dilution = AC_Ratio / (AC_Ratio + 1);
                    Max_Hct = max(Estimation_Hct, Target_End_Hct);
                    Hct_ACWB = AC_Dilution * Max_Hct;
                    Qp_max = Qwb_PD * (1 - (Hct_ACWB / prbc_Hct));
                    if (Qp_PD > Qp_max)
                    {
                        Predicted_Error = "Increase Target End Hct";
                    }


                    if (Qp_PD < 1.3)
                    {
                        if (Qp_PD < 0 | Qp_PD > 1E-07)
                        {
                            Predicted_Error = "Decrease Target End Hct";
                        }
                        else
                        {
                            Qp_PD = 0;
                        }
                    }

                    if (Qrf_PD > 150)
                    {
                        Predicted_Error = "Decrease Target Fluid Balance";
                    }

                    if (Qrf_PD < 1.3)
                    {
                        Predicted_Error = "Decrease Target FCR";
                    }

                    PD_Volume_Remaining = PD_Return_Path_Volume - PD_Replaced_Volume;
                    Est_Mode_Time = PD_Volume_Remaining / (Qrf_PD + Qp_PD);
                    Mode_RF_Volume = Qrf_PD * Est_Mode_Time;
                    Mode_Other_RF_Volume = 0;

                    if (Estimation_Hct > Target_Depletion_Hct)
                    {
                        Plasma_AC_Fraction_Start = 1 / (((1 - Estimation_Hct) * AC_Ratio) + 1);
                        Plasma_AC_Fraction_End = 1 / (((1 - Target_Depletion_Hct) * AC_Ratio) + 1);
                        Plasma_AC_Fraction = (Plasma_AC_Fraction_Start + Plasma_AC_Fraction_End) / 2;
                    }
                    else
                    {
                        Plasma_AC_Fraction = 1 / (((1 - Target_End_Hct) * AC_Ratio) + 1);
                    }

                    Estimation = RBCXTT(Qwb_PD - Qac_PD, Qrf_PD, Qp_PD, 0, Est_Mode_Time, Estimation_Volume, Estimation_Hct, Estimation_FCR);

                    Current_Patient_Volume = Estimation(0);
                    Current_Patient_Hct = Estimation(1);
                    Current_FCR = Estimation(2);

                    Estimation_Volume = Current_Patient_Volume;
                    Estimation_Hct = Current_Patient_Hct;
                    Estimation_FCR = Current_FCR;

                    Predicted_Replacement_Volume = Predicted_Replacement_Volume + Mode_RF_Volume;
                    Predicted_Other_RF_Volume = Predicted_Other_RF_Volume + Mode_Other_RF_Volume;

                    Qwb_Ex = Max_WB_Rate;
                    Qac_Ex = Qwb_Ex / (AC_Ratio + 1);
                    if (Qac_Ex < 1.3)
                    {
                        Predicted_Error = "AC Ratio Unachievable";
                    }
                    if (Replacement_Volume_Override == "Override")
                    {
                        if (Target_Replacement_Volume - RF_Used - Predicted_Replacement_Volume <= 0)
                        {
                            Prdicted_Error = "Increase Replacement Volume";
                        }
                        Calculated_FCR = CalculatedFCRVR(Estimation_Hct, Avg_RF_Hct, Target_End_Hct, Target_Volume - Estimation_Volume, Estimation_Volume, Target_Replacement_Volume - RF_Used - Predicted_Replacement_Volume, (Qwb_Ex - Qac_Ex), 1);
                        Target_FCR = Calculated_FCR * Estimation_FCR;

                        if (Target_FCR < 0)
                        {
                            Predicted_Error = "Decrease Replacement Volume";
                        }
                        if (Target_FCR >= Estimation_FCR)
                        {
                            Predicted_Error = "Decrease Target End Hct";
                        }
                    }
                    else
                    {

                        if (Target_FCR <= 0)
                        {
                            Predicted_Error = "Increase_Target_FCR";
                        }
                        if (Target_FCR >= Estimation_FCR)
                        {
                            Predicted_Error = "Decrease Target FCR";
                        }

                    }

                    if ((Target_FCR / Estimation_FCR) >= (0.99 + ((Target_Volume - Estimation_Volume) / Estimation_Volume)))
                    {
                        Predicted_Error = "Increase Target Fluid Balance";
                    }

                    Qrf_Ex = CalculatedFCR((Qwb_Ex - Qac_Ex), Estimation_Volume, (Target_Volume - Estimation_Volume), Estimation_Hct, Estimation_FCR, Target_End_Hct, Target_FCR, Avg_RF_Hct, 1);
                    Qp_Ex = CalculatedFCR((Qwb_Ex - Qac_Ex), Estimation_Volume, (Target_Volume - Estimation_Volume), Estimation_Hct, Estimation_FCR, Target_End_Hct, Target_FCR, Avg_RF_Hct, 2);
                    Est_Mode_Time = CalculatedFCR((Qwb_Ex - Qac_Ex), Estimation_Volume, (Target_Volume - Estimation_Volume), Estimation_Hct, Estimation_FCR, Target_End_Hct, Target_FCR, Avg_RF_Hct, 3);
                    Max_Hct = max(Estimation_Hct, Target_End_Hct);
                    Plasma_AC_Fraction = 1 / (((1 - Max_Hct) * AC_Ratio) + 1);
                    Estimated_CIR = (Qrf_Ex * Cit_Con_Rep_Fluid) + (Qp_Ex * Cit_Con_AC * Plasma_AC_Fraction);

                    if (Estimated_CIR > CIR_Limit)
                    {
                        CIR_Adjustment = CIR_Limit / Estimated_CIR;
                        Qwb_Ex = Qwb_Ex * CIR_Adjustment;
                        Qac_Ex = Qac_Ex * CIR_Adjustment;
                        if (Qac_Ex < 1.3)
                        {
                            Predicted_Error = "AC Ratio Unachievable";
                        }
                        Qrf_Ex = CalculatedFCR((Qwb_Ex - Qac_Ex), Estimation_Volume, (Target_Volume - Estimation_Volume), Estimation_Hct, Estimation_FCR, Target_End_Hct, Target_FCR, Avg_RF_Hct, 1);
                        Qp_Ex = CalculatedFCR((Qwb_Ex - Qac_Ex), Estimation_Volume, (Target_Volume - Estimation_Volume), Estimation_Hct, Estimation_FCR, Target_End_Hct, Target_FCR, Avg_RF_Hct, 2);
                        Est_Mode_Time = CalculatedFCR((Qwb_Ex - Qac_Ex), Estimation_Volume, (Target_Volume - Estimation_Volume), Estimation_Hct, Estimation_FCR, Target_End_Hct, Target_FCR, Avg_RF_Hct, 3);
                    }

                    prbc_Hct = 0.8;
                    AC_Dilution = AC_Ratio / (AC_Ratio + 1);
                    Max_Hct = max(Estimation_Hct, Target_End_Hct);
                    Hct_ACWB = AC_Dilution * Max_Hct;
                    Qp_max = Qwb_Ex * (1 - (Hct_ACWB / prbc_Hct));

                    if (Qp_Ex > Qp_max)
                    {
                        Predicted_Error = "Unable to meet Hct";
                    }

                    if (Qp_Ex < 1.3)
                    {
                        if (Qp_Ex < 0 | Qp_Ex > 1E-07)
                        {
                            Predicted_Error = "Decrease Target End Hct";
                        }
                        else
                        {
                            Qp_Ex = 0;
                        }
                    }

                    if (Qrf_Ex > 150)
                    {
                        Predicted_Error = "Decrease Target Fluid Balance";
                    }

                    if (Qrf_Ex < 1.3)
                    {
                        Predicted_Error = "Decrease Target FCR";
                    }

                    Mode_RF_Volume = Qrf_Ex * Est_Mode_Time;
                    Mode_Other_RF_Volume = 0;

                    if (Estimation_Hct > Target_End_Hct)
                    {
                        Plasma_AC_Fraction_Start = 1 / (((1 - Estimation_Hct) * AC_Ratio) + 1);
                        Plasma_AC_Fraction_End = 1 / (((1 - Target_End_Hct) * AC_Ratio) + 1);
                        Plasma_AC_Fraction = (Plasma_AC_Fraction_Start + Plasma_AC_Fraction_End) / 2;
                    }
                    else
                    {
                        Plasma_AC_Fraction = 1 / (((1 - Target_End_Hct) * AC_Ratio) + 1);
                    }

                    Mode_Time = Est_Mode_Time;

                    Estimation = RBCXTT(Qwb_Ex - Qac_Ex, Qrf_Ex, Qp_Ex, Avg_RF_Hct, Est_Mode_Time, Estimation_Volume, Estimation_Hct, Estimation_FCR);

                    Current_Patient_Volume = Estimation(0);
                    Current_Patient_Hct = Estimation(1);
                    Current_FCR = Estimation(2);

                    Estimation_Volume = Current_Patient_Volume;
                    Estimation_Hct = Current_Patient_Hct;
                    Estimation_FCR = Current_FCR;

                    Predicted_Replacement_Volume = Predicted_Replacement_Volume + Mode_RF_Volume;
                    Predicted_Other_RF_Volume = Predicted_Other_RF_Volume + Mode_Other_RF_Volume;

                }
            }

            if (Predicted_Error == "None")
            {
                if (Procedure_Type == "Depletion")
                {
                    WB_Flow_Rate = Qwb_Dep;
                }
                else
                {
                    WB_Flow_Rate = Qwb_Ex;
                }
                if (Replacement_Volume_Override == "Override")
                {
                    FCR = Estimation_FCR;
                }
                else
                {
                    Replacement_Volume = Predicted_Replacement_Volume;
                }
                Est_Other_RF = Predicted_Other_RF_Volume;
                Predicted_Error = Predicted_Error;
            }
            else
            {
                WB_Flow_Rate = 0;

                if (Replacement_Volume_Override == "Override")
                {
                    FCR = 0;
                }
                else
                {
                    Replacement_Volume = 0;
                }

                Est_Other_RF = 0;
                Predicted_Error = Predicted_Error;
            }

            if (Flag == 0)
                functionReturnValue = (Math.Pow((Estimation_FCR * (Entered_Patient_Volume / Estimation_Volume) * (Entered_Patient_Hct / Estimation_Hct)), (FCR_Factor))) * 100;
            if (Flag == 1)
                functionReturnValue = Replacement_Volume;
            if (Flag == 2)
                functionReturnValue = Est_Other_RF;
            if (Flag == 3)
                functionReturnValue = Predicted_Error;
            return functionReturnValue;

        }



    }
}
