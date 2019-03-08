namespace Zygotine.Statistics.Distribution
{
    using System;
    using Zygotine.Numerics;
    using Zygotine.Statistics.Util;
    using Zygotine.Util;

    public static partial class PoissonDistribution
    {

        /*
         * R rile
         * dpois.c
         */

        public static double DPoisson(double x, double lambda, bool giveLog)
        {
            if (double.IsNaN(x) || double.IsNaN(lambda))
                return x + lambda;

            if (lambda < 0)
                return double.NaN; // ML_ERR_return_NAN;

            if ((Math.Abs((x) - Math.Floor(x)) > 1e-7 * RVaria.fmax2(1.0, Math.Abs(x))))
            {
                //MATHLIB_WARNING("non-integer x = %f", x);	
                return (giveLog ? double.NegativeInfinity : 0.0);
            }

            if (x < 0 || !x.IsFinite())
                return (giveLog ? double.NegativeInfinity : 0.0);

            x = Math.Round(x);
            return (dpois_raw(x, lambda, giveLog));
        }

        internal static double dpois_raw(double x, double lambda, bool log_p)
        {

            /*       x >= 0 ; integer for dpois(), but not e.g. for pgamma()!
                lambda >= 0
            */
            if (lambda == 0) return ((x == 0) ? (log_p ? 0.0 : 1.0) : (log_p ? double.NegativeInfinity : 0.0));
            if (!lambda.IsFinite()) return (log_p ? double.NegativeInfinity : 0.0);
            if (x < 0) return ((log_p ? double.NegativeInfinity : 0.0));
            if (x <= lambda * RVaria.DBL_MIN) return ((log_p ? (-lambda) : Math.Exp(-lambda)));
            if (lambda < x * RVaria.DBL_MIN)
                return (log_p
                    ? (-lambda + x * Math.Log(lambda) - Functions.LogGamma(x + 1))
                    : Math.Exp(-lambda + x * Math.Log(lambda) - Functions.LogGamma(x + 1)));

            return (RVaria.R_D_fexp(RVaria.M_2PI * x, -Functions.StirlingError(x) - RVaria.bd0(x, lambda), log_p));
        }


    }
}
