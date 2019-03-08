namespace Zygotine.Numerics
{
    using System;
    using Zygotine.Statistics.Util;

    public static class Functions
    {

        /* R file: gamma.c
         * name used in R: gammafn
         */
        public static double Gamma(double x)
        {


            int i, n;
            double y;
            double sinpiy, value;


            const double
                xmin = -170.5674972726612,
                xmax = 171.61447887182298,
                xsml = 2.2474362225598545e-308,
                dxrel = 1.490116119384765696e-8;
            const int
                ngam = 22;


            if (double.IsNaN(x)) return x;

            /* If the argument is exactly zero or a negative integer
             * then return NaN. */
            if (x == 0 || (x < 0 && x == Math.Round(x)))
            {
                // ML_ERROR(ME_DOMAIN, "gammafn");
                return double.NaN;
            }

            y = Math.Abs(x);

            if (y <= 10)
            {

                /* Compute gamma(x) for -10 <= x <= 10
                 * Reduce the interval and find gamma(1 + y) for 0 <= y < 1
                 * first of all. */

                n = (int)x;
                if (x < 0) --n;
                y = x - n;/* n = floor(x)  ==>	y in [ 0, 1 ) */
                --n;
                value = RVaria.chebyshev_eval(y * 2 - 1, gamcs, ngam) + .9375;
                if (n == 0)
                    return value;/* x = 1.dddd = 1+y */

                if (n < 0)
                {
                    /* compute gamma(x) for -10 <= x < 1 */

                    /* exact 0 or "-n" checked already above */

                    /* The answer is less than half precision */
                    /* because x too near a negative integer. */
                    if ((x < -0.5) && (Math.Abs(x - (int)(x - 0.5) / x) < dxrel))
                    {
                        // ML_ERROR(ME_PRECISION, "gammafn");
                    }

                    /* The argument is so close to 0 that the result would overflow. */
                    if (y < xsml)
                    {
                        // ML_ERROR(ME_RANGE, "gammafn");
                        if (x > 0) return double.PositiveInfinity;
                        else return double.NegativeInfinity;
                    }

                    n = -n;

                    for (i = 0; i < n; i++)
                    {
                        value /= (x + i);
                    }
                    return value;
                }
                else
                {
                    /* gamma(x) for 2 <= x <= 10 */

                    for (i = 1; i <= n; i++)
                    {
                        value *= (y + i);
                    }
                    return value;
                }
            }
            else
            {
                /* gamma(x) for	 y = |x| > 10. */

                if (x > xmax)
                {			/* Overflow */
                    // ML_ERROR(ME_RANGE, "gammafn");
                    return double.PositiveInfinity;
                }

                if (x < xmin)
                {			/* Underflow */
                    // ML_ERROR(ME_UNDERFLOW, "gammafn");
                    return 0.0;
                }

                if (y <= 50 && y == (int)y)
                { /* compute (n - 1)! */
                    value = 1.0;
                    for (i = 2; i < y; i++) value *= i;
                }
                else
                { /* normal case */
                    value = Math.Exp((y - 0.5) * Math.Log(y) - y + RVaria.M_LN_SQRT_2PI +
                        ((2 * y == (int)2 * y) ? StirlingError(y) : lgammacor(y)));
                }
                if (x > 0)
                    return value;

                if (Math.Abs((x - (int)(x - 0.5)) / x) < dxrel)
                {

                    /* The answer is less than half precision because */
                    /* the argument is too near a negative integer. */

                    // ML_ERROR(ME_PRECISION, "gammafn");
                }

                sinpiy = RVaria.sinpi(y);
                if (sinpiy == 0)
                {		/* Negative integer arg - overflow */
                    // ML_ERROR(ME_RANGE, "gammafn");
                    return double.PositiveInfinity;
                }

                return -RVaria.M_PI / (y * sinpiy * value);
            }
        }

        /* R file: lgamma.c
         * name used in R: lgammafn
         */
        public static double LogGamma(double x)
        {
            int? sign = null;
            return lgammafn_sign(x, ref sign);
        }

        /* R file: stirlerr.c
         * name used in R: stirlerr
         */
        public static double StirlingError(double n)
        {

            /*
              error for 0, 0.5, 1.0, 1.5, ..., 14.5, 15.0.
            */

            double nn;

            if (n <= 15.0)
            {
                nn = n + n;
                if (nn == (int)nn) return (sferr_halves[(int)nn]);
                return (LogGamma(n + 1.0) - (n + 0.5) * Math.Log(n) + n - RVaria.M_LN_SQRT_2PI);
            }

            nn = n * n;
            if (n > 500) return ((S0 - S1 / nn) / n);
            if (n > 80) return ((S0 - (S1 - S2 / nn) / nn) / n);
            if (n > 35) return ((S0 - (S1 - (S2 - S3 / nn) / nn) / nn) / n);
            /* 15 < n <= 35 : */
            return ((S0 - (S1 - (S2 - (S3 - S4 / nn) / nn) / nn) / nn) / n);
        }

        // R file: gamma.c
        #region double[] gamcs
        private static double[] gamcs = new double[42] {
            +.8571195590989331421920062399942e-2,
            +.4415381324841006757191315771652e-2,
            +.5685043681599363378632664588789e-1,
            -.4219835396418560501012500186624e-2,
            +.1326808181212460220584006796352e-2,
            -.1893024529798880432523947023886e-3,
            +.3606925327441245256578082217225e-4,
            -.6056761904460864218485548290365e-5,
            +.1055829546302283344731823509093e-5,
            -.1811967365542384048291855891166e-6,
            +.3117724964715322277790254593169e-7,
            -.5354219639019687140874081024347e-8,
            +.9193275519859588946887786825940e-9,
            -.1577941280288339761767423273953e-9,
            +.2707980622934954543266540433089e-10,
            -.4646818653825730144081661058933e-11,
            +.7973350192007419656460767175359e-12,
            -.1368078209830916025799499172309e-12,
            +.2347319486563800657233471771688e-13,
            -.4027432614949066932766570534699e-14,
            +.6910051747372100912138336975257e-15,
            -.1185584500221992907052387126192e-15,
            +.2034148542496373955201026051932e-16,
            -.3490054341717405849274012949108e-17,
            +.5987993856485305567135051066026e-18,
            -.1027378057872228074490069778431e-18,
            +.1762702816060529824942759660748e-19,
            -.3024320653735306260958772112042e-20,
            +.5188914660218397839717833550506e-21,
            -.8902770842456576692449251601066e-22,
            +.1527474068493342602274596891306e-22,
            -.2620731256187362900257328332799e-23,
            +.4496464047830538670331046570666e-24,
            -.7714712731336877911703901525333e-25,
            +.1323635453126044036486572714666e-25,
            -.2270999412942928816702313813333e-26,
            +.3896418998003991449320816639999e-27,
            -.6685198115125953327792127999999e-28,
            +.1146998663140024384347613866666e-28,
            -.1967938586345134677295103999999e-29,
            +.3376448816585338090334890666666e-30,
            -.5793070335782135784625493333333e-31
        };
        #endregion

        #region R file: stirlerr.c

        private static readonly double[] sferr_halves = new double[31]
        {
                0.0, /* n=0 - wrong, place holder only */
	            0.1534264097200273452913848,  /* 0.5 */
	            0.0810614667953272582196702,  /* 1.0 */
	            0.0548141210519176538961390,  /* 1.5 */
	            0.0413406959554092940938221,  /* 2.0 */
	            0.03316287351993628748511048, /* 2.5 */
	            0.02767792568499833914878929, /* 3.0 */
	            0.02374616365629749597132920, /* 3.5 */
	            0.02079067210376509311152277, /* 4.0 */
	            0.01848845053267318523077934, /* 4.5 */
	            0.01664469118982119216319487, /* 5.0 */
	            0.01513497322191737887351255, /* 5.5 */
	            0.01387612882307074799874573, /* 6.0 */
	            0.01281046524292022692424986, /* 6.5 */
	            0.01189670994589177009505572, /* 7.0 */
	            0.01110455975820691732662991, /* 7.5 */
	            0.010411265261972096497478567, /* 8.0 */
	            0.009799416126158803298389475, /* 8.5 */
	            0.009255462182712732917728637, /* 9.0 */
	            0.008768700134139385462952823, /* 9.5 */
	            0.008330563433362871256469318, /* 10.0 */
	            0.007934114564314020547248100, /* 10.5 */
	            0.007573675487951840794972024, /* 11.0 */
	            0.007244554301320383179543912, /* 11.5 */
	            0.006942840107209529865664152, /* 12.0 */
	            0.006665247032707682442354394, /* 12.5 */
	            0.006408994188004207068439631, /* 13.0 */
	            0.006171712263039457647532867, /* 13.5 */
	            0.005951370112758847735624416, /* 14.0 */
	            0.005746216513010115682023589, /* 14.5 */
	            0.005554733551962801371038690  /* 15.0 */
        };


        private const double
            S0 = 0.083333333333333333333,       /* 1/12 */
            S1 = 0.00277777777777777777778,     /* 1/360 */
            S2 = 0.00079365079365079365079365,  /* 1/1260 */
            S3 = 0.000595238095238095238095238, /* 1/1680 */
            S4 = 0.0008417508417508417508417508;/* 1/1188 */

        #endregion 

        /* R file: lgamma.c
         */
        private static double lgammafn_sign(double x, ref int? sgn)
        {
            double ans, y, sinpiy;

            const double
                xmax = 2.5327372760800758e+305,
                dxrel = 1.490116119384765625e-8;

            if (sgn != null) sgn = 1;

            if (double.IsNaN(x))
                return x;

            if ((sgn != null) && (x < 0) && ((Math.Floor(-x) % 2.0) == 0)) // fmod était utilisé
                sgn = -1;

            if (x <= 0 && x == Math.Truncate(x))
            { /* Negative integer argument */
                // ML_ERROR(ME_RANGE, "lgamma");
                return double.PositiveInfinity;/* +Inf, since lgamma(x) = log|gamma(x)| */
            }

            y = Math.Abs(x);

            if (y < 1e-306) return -Math.Log(y); // denormalized range, R change
            if (y <= 10) return Math.Log(Math.Abs(Gamma(x)));
            /*
              ELSE  y = |x| > 10 ---------------------- */

            if (y > xmax)
            {
                // ML_ERROR(ME_RANGE, "lgamma");
                return double.PositiveInfinity;
            }

            if (x > 0)
            { /* i.e. y = x > 10 */
                if (x > 1e17)
                    return (x * (Math.Log(x) - 1.0));
                else if (x > 4934720.0)
                    return (RVaria.M_LN_SQRT_2PI + (x - 0.5) * Math.Log(x) - x);
                else
                    return RVaria.M_LN_SQRT_2PI + (x - 0.5) * Math.Log(x) - x + lgammacor(x);
            }
            /* else: x < -10; y = -x */
            sinpiy = Math.Abs(RVaria.sinpi(y));

            if (sinpiy == 0)
            { /* Negative integer argument ===
			  Now UNNECESSARY: caught above */
                // MATHLIB_WARNING(" ** should NEVER happen! *** [lgamma.c: Neg.int, y=%g]\n",y);
                return RVaria.R_NaN;
                // ML_ERR_return_NAN;
            }

            ans = RVaria.M_LN_SQRT_PId2 + (x - 0.5) * Math.Log(y) - x - Math.Log(sinpiy) - lgammacor(y);

            if (Math.Abs((x - Math.Truncate(x - 0.5)) * ans / x) < dxrel)
            {

                /* The answer is less than half precision because
                 * the argument is too near a negative integer. */

                // ML_ERROR(ME_PRECISION, "lgamma");
            }

            return ans;
        }

        // R file: lgammacor.c
        #region private static readonly double[] algmcs
        private static readonly double[] algmcs = new double[15]{
            +.1666389480451863247205729650822e+0,
            -.1384948176067563840732986059135e-4,
            +.9810825646924729426157171547487e-8,
            -.1809129475572494194263306266719e-10,
            +.6221098041892605227126015543416e-13,
            -.3399615005417721944303330599666e-15,
            +.2683181998482698748957538846666e-17,
            -.2868042435334643284144622399999e-19,
            +.3962837061046434803679306666666e-21,
            -.6831888753985766870111999999999e-23,
            +.1429227355942498147573333333333e-24,
            -.3547598158101070547199999999999e-26,
            +.1025680058010470912000000000000e-27,
            -.3401102254316748799999999999999e-29,
            +.1276642195630062933333333333333e-30
    };
        #endregion

        // R file: lgammacor.c
        private static double lgammacor(double x)
        {

            double tmp;

            /* For IEEE double precision DBL_EPSILON = 2^-52 = 2.220446049250313e-16 :
             *   xbig = 2 ^ 26.5
             *   xmax = DBL_MAX / 48 =  2^1020 / 3 
             */

            int nalgm = 5;
            double
                xbig = 94906265.62425156,
                xmax = 3.745194030963158e306;

            if (x < 10)
                // ML_ERR_return_NAN
                return RVaria.R_NaN;
            else if (x >= xmax)
            {
                // ML_ERROR(ME_UNDERFLOW, "lgammacor");
                /* allow to underflow below */
            }
            else if (x < xbig)
            {
                tmp = 10 / x;
                return RVaria.chebyshev_eval(tmp * tmp * 2 - 1, algmcs, nalgm) / x;
            }
            return 1 / (x * 12);
        }
    }


}
