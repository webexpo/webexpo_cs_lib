namespace Zygotine.Statistics.Distribution
{
    using System;
    using Zygotine.Statistics.Util;
    using Zygotine.Util;

    public static partial class NormalDistribution
    {
        /*
         * 
         * R file
         * dnorm.c
         * 
         */

        /*
         * Pdf 
         */
        public static double DNorm(double x, double mu = 0, double sigma = 1, bool give_log = false)
        {
            // In dnorm.c, the function name is dnorm4!
            if (double.IsNaN(x) || double.IsNaN(mu) || double.IsNaN(sigma))
                return RVaria.R_NaN; // return x + mu + sigma;
            if (!sigma.IsFinite()) return (give_log ? double.NegativeInfinity : 0.0);
            if (!x.IsFinite() && mu == x) return RVaria.R_NaN;/* x-mu is NaN */
            if (sigma <= 0)
            {
                if (sigma < 0)
                    // ML_ERR_return_NAN;
                    return RVaria.R_NaN;
                /* sigma == 0 */
                return (x == mu) ? double.PositiveInfinity : (give_log ? double.NegativeInfinity : 0.0);
            }
            x = (x - mu) / sigma;

            if (!x.IsFinite()) return (give_log ? double.NegativeInfinity : 0.0);

            x = Math.Abs(x);
            if (x >= 2 * Math.Sqrt(RVaria.DBL_MAX)) return (give_log ? double.NegativeInfinity : 0.0);
            if (give_log)
                return -(RVaria.M_LN_SQRT_2PI + 0.5 * x * x + Math.Log(sigma));
            //  M_1_SQRT_2PI = 1 / sqrt(2 * pi)

            //#ifdef MATHLIB_FAST_dnorm
            //    // and for R <= 3.0.x and R-devel upto 2014-01-01:
            //    return M_1_SQRT_2PI * exp(-0.5 * x * x) / sigma;
            //#else

            // more accurate, less fast :
            if (x < 5) return RVaria.M_1_SQRT_2PI * Math.Exp(-0.5 * x * x) / sigma;

            /* ELSE:

             * x*x  may lose upto about two digits accuracy for "large" x
             * Morten Welinder's proposal for PR#15620
             * https://bugs.r-project.org/bugzilla/show_bug.cgi?id=15620

             * -- 1 --  No hoop jumping when we underflow to zero anyway:

             *  -x^2/2 <         log(2)*.Machine$double.min.exp  <==>
             *     x   > sqrt(-2*log(2)*.Machine$double.min.exp) =IEEE= 37.64031
             * but "thanks" to denormalized numbers, underflow happens a bit later,
             *  effective.D.MIN.EXP <- with(.Machine, double.min.exp + double.ulp.digits)
             * for IEEE, DBL_MIN_EXP is -1022 but "effective" is -1074
             * ==> boundary = sqrt(-2*log(2)*(.Machine$double.min.exp + .Machine$double.ulp.digits))
             *              =IEEE=  38.58601
             * [on one x86_64 platform, effective boundary a bit lower: 38.56804]
             */
            if (x > Math.Sqrt(-2 * RVaria.M_LN2 * (RVaria.DBL_MIN_EXP + 1 - RVaria.DBL_MANT_DIG))) return 0.0;

            /* Now, to get full accurary, split x into two parts,
             *  x = x1+x2, such that |x2| <= 2^-16.
             * Assuming that we are using IEEE doubles, that means that
             * x1*x1 is error free for x<1024 (but we have x < 38.6 anyway).

             * If we do not have IEEE this is still an improvement over the naive formula.
             */
            double x1 = //  R_forceint(x * 65536) / 65536 =
            RVaria.ldexp(Math.Round(RVaria.ldexp(x, 16)), -16);
            double x2 = x - x1;
            return RVaria.M_1_SQRT_2PI / sigma *
            (Math.Exp(-0.5 * x1 * x1) * Math.Exp((-0.5 * x2 - x1) * x2));

        }

        public static double[] DNorm(double[] x, double mu = 0, double sigma = 1, bool give_log = false)
        {
            double[] rep = new double[x.Length];
            for(int i= 0; i<x.Length; i++)
            {
                rep[i] = DNorm(x[i], mu, sigma, give_log);
            }
            return rep;
        }

        private const double SIXTEN = 16; /* Cutoff allowing exact "*" and "/" */
    }
}
