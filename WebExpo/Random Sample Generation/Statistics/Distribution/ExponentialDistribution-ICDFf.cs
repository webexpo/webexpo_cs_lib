
namespace Zygotine.Statistics.Distribution
{
    using System;
    using System.Linq;
    using Zygotine.Statistics.Util;

    public static partial class ExponentialDistribution
    {

        /* R file: qexp.c
         * 
         */

        /* Icdf 
         */

        public static double[] QExp(double[] p, double scale = 1.0, bool lowerTail = true, bool logP = false)
        {
            return p.Select(x => QExp(x, scale, lowerTail, logP)).ToArray();
        }

        public static double QExp(double p, double scale = 1.0, bool lowerTail = true, bool logP = false)
        {
            if (double.IsNaN(p) || double.IsNaN(scale))
                return p + scale;
            if (scale < 0) return double.NaN; //ML_ERR_return_NAN;

            if ((logP && p > 0) || (!logP && (p < 0 || p > 1)))
            //ML_ERR_return_NAN
            {
                return RVaria.R_NaN;
            }
            if (p == (lowerTail ? (logP ? double.NegativeInfinity : 0.0) : (logP ? 0.0 : 1.0)))
                return 0;

            return -scale * (lowerTail ? (logP ? ((p) > -RVaria.M_LN2 ? Math.Log(-RVaria.expm1(p)) : RVaria.log1p(-Math.Exp(p))) : RVaria.log1p(-p)) : (logP ? (p) : Math.Log(p)));
        }









    }
   

}

