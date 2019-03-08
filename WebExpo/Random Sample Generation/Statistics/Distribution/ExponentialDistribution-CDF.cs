namespace Zygotine.Statistics.Distribution
{
    using System;
    using Zygotine.Statistics.Util;

    public static partial class ExponentialDistribution
    {
        /*
         * R file
         * pexp.c
         * 
         */

        public static double PExp(double x, double scale = 1.0, bool lowerTail = true, bool logP = false)
        {
            return PExp(new double[] { x }, scale, lowerTail, logP)[0];
        }

        public static double[] PExp(double[] xA, double scale = 1.0, bool lowerTail = true, bool logP = false)
        {
            double[] rep = new double[xA.Length];

            for (int i = 0; i < xA.Length; i++)
            {
                double x = xA[i];
                if (double.IsNaN(x) || double.IsNaN(scale))
                {
                    rep[i] = x + scale;
                }

                if (scale < 0)
                {
                    rep[i] = double.NaN;
                }

                if (x <= 0.0)
                    rep[i] = (lowerTail ? (logP ? double.NegativeInfinity : 0.0) : (logP ? 0.0 : 1.0));
                /* same as weibull( shape = 1): */
                x = -(x / scale);
                if (lowerTail)
                    rep[i] = (logP
                    ? (x > -RVaria.M_LN2 ? Math.Log(-RVaria.expm1(x)) : RVaria.log1p(-Math.Exp(x)))
                    : -RVaria.expm1(x));
                /* else:  !lower_tail */
                else
                    rep[i] = (logP ? (x) : Math.Exp(x));
            }

            return rep;
        }
    }
}

