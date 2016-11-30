using System;

namespace Formula
{
    public class Temp
    {

        public double CalculatedFCR(double inletRate, double exchangeVolume, double exchangeVolumeChange,
            double exchangeHct, double exchangeFcr, double exchangeEndHct, double exchangeTargetFcr,
            double exchangeRfHct, int flag)
        {

            double therapeuticFcr = exchangeTargetFcr/exchangeFcr;
            double m = 0;
            double p;

            if (Math.Abs(exchangeVolumeChange) > 0.0000001)
            {


                m = 1 - Math.Log(therapeuticFcr)/(Math.Log(1 + (exchangeVolumeChange/exchangeVolume)));
                p = Math.Pow(therapeuticFcr, m/(m - 1));
            }
            else
            {
                p = therapeuticFcr;
            }

            var f = exchangeRfHct/exchangeHct*((p - 1)/(p - (exchangeEndHct/exchangeHct)));
            double qrf;

            if (Math.Abs(exchangeVolumeChange) > 0.0000001)
            {

                qrf = inletRate*(m/(f*(m - 1)));
            }
            else
            {
                qrf = inletRate/f;
            }

            double qp = (f - 1)*qrf;
            double exchangeTime;

            if (Math.Abs(exchangeVolumeChange) > 0.0000001)
            {
                double flowRateDifference = qrf + qp - inletRate;
                m = (qrf + qp)/flowRateDifference;
                exchangeTime = exchangeVolume/flowRateDifference*
                                ((1/Math.Pow(therapeuticFcr, (1/(m - 1)))) - 1);
            }
            else
            {
                exchangeTime = -exchangeVolume/inletRate*Math.Log(therapeuticFcr);
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
    }


}
