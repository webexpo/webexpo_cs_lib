
namespace Zygotine.Statistics.Distribution
{
    using System;

    public static partial class GammaDistribution
    {
        /* R file: dgamma.c
        */

        /* Pdf 
         */
        public static double DGammaFromScale(double x, double shape, double scale, bool giveLog)
        {
            double pr;

            if (double.IsNaN(x) || double.IsNaN(shape) || double.IsNaN(scale))
                return x + shape + scale;
            if (shape < 0 || scale <= 0)
                return double.NaN; // ML_ERR_return_NAN;
            if (x < 0)
                return (giveLog ? double.NegativeInfinity : 0.0);
            if (shape == 0) /* point mass at 0 */
                return (x == 0) ? double.NegativeInfinity : (giveLog ? double.NegativeInfinity : 0.0);
            if (x == 0)
            {
                if (shape < 1) return double.NegativeInfinity;
                if (shape > 1) return (giveLog ? double.NegativeInfinity : 0.0);
                /* else */
                return giveLog ? -Math.Log(scale) : 1 / scale;
            }

            if (shape < 1)
            {
                pr = PoissonDistribution.dpois_raw(shape, x / scale, giveLog);
                return giveLog ? pr + Math.Log(shape / x) : pr * shape / x;
            }
            /* else  shape >= 1 */
            pr = PoissonDistribution.dpois_raw(shape - 1, x / scale, giveLog);
            return giveLog ? pr - Math.Log(scale) : pr / scale;
        }


        public static double DGammaFromRate(double x, double shape, double rate, bool giveLog)
        {
            if (rate == 0)
            {
                return double.NaN;
            }

            double scale = 1 / rate;
            return DGammaFromScale(x, shape, scale, giveLog);
        }
    }

}
