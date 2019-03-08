using System;
//using Extreme.FloatingPoint;
using Zygotine.Util;

namespace Zygotine.Statistics.Util
{
    public static class RVaria
    {

        //internal static 
        internal const double M_LN2 = 0.693147180559945309417; // math.h
        internal const double DBL_MAX_EXP = 1024; // float.h ieee-754
        internal const double DBL_MAX = 1.7976931348623158e+308; // C float.h ieee-754
        internal const double DBL_MIN = 2.2250738585072014E-308; // C float.h  Smallest positive normalized double- voir FPExtreme
        internal const int DBL_MIN_EXP = -1021; // C float.h
        internal const int DBL_MANT_DIG = 53; // C float.h
        internal const double DBL_EPSILON = 2.220446049250313e-16; // float.h.  DBL_EPSILON est le plus petit double positif tel que 1.0 + x > 1.0.
        internal const double R_AccuracyInfo_eps = DBL_EPSILON;
        internal const double M_2PI = 6.283185307179586476925286766559; // Rmath.h
        internal const double M_1_SQRT_2PI = 0.398942280401432677939946059934; //Rmath.h ...	/* 1/sqrt(2pi) */
        internal const double M_LN_SQRT_2PI = 0.918938533204672741780329736406;	// log(sqrt(2*pi)) == log(2*pi)/2  // Rmath.h
        internal const double M_LN_SQRT_PId2 = 0.225791352644727432363097614947;	// log(sqrt(pi/2)) // RMath.h
        internal const double M_SQRT2 = 1.41421356237309504880;
        internal const double M_SQRT1_2 = 0.707106781186547524401; // 1/sqrt(2)  
        internal const double M_SQRT_32 = 5.656854249492380195206754896838; // Rmath.h
        internal const double M_PI = 3.141592653589793238462643383279502884197169399375; // Constants.h
        internal const double R_NaN = double.NaN;
        //ajouts pour digamma mai 2017
        internal const int CHAR_BIT = 8;
        internal const int INT_MAX = 2147483647;
        internal const int FLT_RADIX = 2;
        internal const double M_LOG10_2 = 0.301029995663981195213738894724;



        /*
         * R file: bd0.c
         * 
         *  AUTHOR
         *	Catherine Loader, catherine@research.bell-labs.com.
         *	October 23, 2000.
         *
         *  Merge in to R:
         *	Copyright (C) 2000, The R Core Team
         *
         *  This program is free software; you can redistribute it and/or modify
         *  it under the terms of the GNU General Public License as published by
         *  the Free Software Foundation; either version 2 of the License, or
         *  (at your option) any later version.
         *
         *  This program is distributed in the hope that it will be useful,
         *  but WITHOUT ANY WARRANTY; without even the implied warranty of
         *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
         *  GNU General Public License for more details.
         *
         *  You should have received a copy of the GNU General Public License
         *  along with this program; if not, a copy is available at
         *  http://www.r-project.org/Licenses/
         *
         *
         *  DESCRIPTION
         *	Evaluates the "deviance part"
         *	bd0(x,M) :=  M * D0(x/M) = M*[ x/M * log(x/M) + 1 - (x/M) ] =
         *		  =  x * log(x/M) + M - x
         *	where M = E[X] = n*p (or = lambda), for	  x, M > 0
         *
         *	in a manner that should be stable (with small relative error)
         *	for all x and M=np. In particular for x/np close to 1, direct
         *	evaluation fails, and evaluation is based on the Taylor series
         *	of log((1+v)/(1-v)) with v = (x-M)/(x+M) = (x-np)/(x+np).
         */


        static internal double bd0(double x, double np)
        {
            double ej, s, s1, v;
            int j;

            if (!x.IsFinite() || !np.IsFinite() || np == 0.0)
                //ML_ERR_return_NAN;
                return R_NaN;

            if (Math.Abs(x - np) < 0.1 * (x + np))
            {
                v = (x - np) / (x + np);  // might underflow to 0
                s = (x - np) * v;/* s using v -- change by MM */
                if (Math.Abs(s) < DBL_MIN) return s;
                ej = 2 * x * v;
                v = v * v;
                for (j = 1; j < 1000; j++)
                { /* Taylor series; 1000: no infinite loop
					as |v| < .1,  v^2000 is "zero" */
                    ej *= v;// = v^(2j+1)
                    s1 = s + ej / ((j << 1) + 1);
                    if (s1 == s) /* last term was effectively 0 */
                        return s1;
                    s = s1;
                }
            }
            /* else:  | x - np |  is not too small */
            return (x * Math.Log(x / np) + np - x);
        }

        // C language
        // returns x * Pow(2, n)
        static internal double ldexp(double x, int n)
        {
            return Tools.scalbln(x, n);
        }

        // R file: cospi.c
        internal static double sinpi(double x)
        {
            if (!x.IsFinite())
                return RVaria.R_NaN;

            // x = Math.mod(x, 2.0); // sin(pi(x + 2k)) == sin(pi x)  for all integer k
            // en C#
            x = x % 2.0;
            // map (-2,2) --> (-1,1] :
            if (x <= -1) x += 2.0; else if (x > 1.0) x -= 2.0;
            if (x == 0.0 || x == 1.0) return 0.0;
            if (x == 0.5) return 1.0;
            if (x == -0.5) return -1.0;
            // otherwise
            return Math.Sin(M_PI * x);
        }

        /* 
         * R file: arithmetic.c
         *
         * Méthode inutilisée - on conserve à des fins de documentation
         * Prendre note que le NA que nous utilisons diffère de celui de R.
         * Prendre égalament note que la traduction faite du code C ne correspond pas exactement à l'original.
         * 
         * int R_IsNA(double x) 
         * {
         *      if (isnan(x)) {
         *          ieee_double y;
         *          y.value = x;
         *          return (y.word[lw] == 1954);
         *      }
         *      return 0;
         * 
         * 
         */
        internal static bool R_IsNA(double x)
        {
            return (double.IsNaN(x) && (Tools.Mantissa(x) & 1954) == 1954);
        }

        // R file: arithmetic.c
        //Méthode inutilisée - on conserve à des fins de documentation 
        internal static bool R_IsNaN(double x)
        {
            return double.IsNaN(x);
        }

        // R file: arith.h
        // Méthode inutilisée - on conserve à des fins de documentation 
        internal static bool R_FINITE(double x)
        {
            /*
             * ISNA(x)        True for R's NA only
             * ISNAN(x)       True for R's NA and IEEE NaN
             * R_FINITE(x)    False for Inf, -Inf, NA, NaN
             */
            return x.IsFinite();
        }

        // R file: nmath.h
        // #define ISNAN(x) (isnan(x)!=0)
        // Méthode inutilisée - on conserve à des fins de documentation 
        internal static bool ISNAN(double x)
        {
            return double.IsNaN(x);
        }

        // R file: fsign.c
        internal static double fsign(double x, double y)
        {
            if (double.IsNaN(x) || double.IsNaN(y))
                return x + y;
            return ((y >= 0) ? Math.Abs(x) : -Math.Abs(x));
        }

        // R file: fmin2.c
        internal static double fmin2(double x, double y)
        {
            if (double.IsNaN(x) || double.IsNaN(y))
                return x + y;
            return (x < y) ? x : y;
        }

        // R file: fmax2.c
        internal static double fmax2(double x, double y)
        {
            if (double.IsNaN(x) || double.IsNaN(y))
                return x + y;
            return (x < y) ? y : x;
        }

        // R file: expm1.c
        internal static double expm1(double x)
        {
            double y, a = Math.Abs(x);

            if (a < DBL_EPSILON) return x;
            if (a > 0.697) return Math.Exp(x) - 1;  /* negligible cancellation */

            if (a > 1e-8)
                y = Math.Exp(x) - 1;
            else /* Taylor expansion, more accurate in this range */
                y = (x / 2 + 1) * x;

            /* Newton step for solving   log(1 + y) = x   for y : */
            /* WARNING: does not work for y ~ -1: bug in 1.5.0 */
            y -= (1 + y) * (log1p(y) - x);
            return y;
        }

        // R file: log1p.c
        internal static double log1p(double x)
        {
            /* series for log1p on the interval -.375 to .375
             *				     with weighted error   6.35e-32
             *				      log weighted error  31.20
             *			    significant figures required  30.93
             *				 decimal places required  32.01
             */


            int nlnrel = 22;
            double xmin = -0.999999985;

            //if (x == 0.0) return 0.0;/* speed */
            //if (x == -1) return(ML_NEGINF);
            //if (x  < -1) ML_ERR_return_NAN;

            if (Math.Abs(x) <= .375)
            {
                /* Improve on speed (only);
               again give result accurate to IEEE double precision: */
                if (Math.Abs(x) < .5 * DBL_EPSILON)
                    return x;

                if ((0 < x && x < 1e-8) || (-1e-9 < x && x < 0))
                    return x * (1 - .5 * x);
                /* else */
                return x * (1 - x * chebyshev_eval(x / .375, alnrcs, nlnrel));
            }
            /* else */
            if (x < xmin)
            {
                /* answer less than half precision because x too near -1 */
                throw new Exception("log1p: ME_PRECISION");
            }
            return Math.Log(1 + x);
        }


        // R file: log1p.c
        internal static double hypot(double a, double b)
        {
            double p, r, s, t, tmp, u;

            if (ISNAN(a) || ISNAN(b)) /* propagate Na(N)s: */
            {
                return a + b;
            }

            if (!R_FINITE(a) || !R_FINITE(b))
            {
                return double.PositiveInfinity;
            }
            p = fmax2(Math.Abs(a), Math.Abs(b));
            if (p != 0.0)
            {
                /* r = (min(|a|,|b|) / p) ^2 */
                tmp = fmin2(Math.Abs(a), Math.Abs(b)) / p;
                r = tmp * tmp;
                for (;;)
                {
                    t = 4.0 + r;
                    /* This was a test of 4.0 + r == 4.0, but optimizing
                    compilers nowadays infinite loop on that. */
                    if (Math.Abs(r) < 2 * DBL_EPSILON) break;
                    s = r / t;
                    u = 1.0 + 2.0 * s;
                    p *= u;

                    /* r = (s / u)^2 * r */
                    tmp = s / u;
                    r *= tmp * tmp;
                }
            }
            return p;
        }

        static readonly double[] alnrcs = new double[43]
        {
            +0.10378693562743769800686267719098e+1,
            -0.13364301504908918098766041553133e+0,
            +0.19408249135520563357926199374750e-1,
            -0.30107551127535777690376537776592e-2,
            +0.48694614797154850090456366509137e-3,
            -0.81054881893175356066809943008622e-4,
            +0.13778847799559524782938251496059e-4,
            -0.23802210894358970251369992914935e-5,
            +0.41640416213865183476391859901989e-6,
            -0.73595828378075994984266837031998e-7,
            +0.13117611876241674949152294345011e-7,
            -0.23546709317742425136696092330175e-8,
            +0.42522773276034997775638052962567e-9,
            -0.77190894134840796826108107493300e-10,
            +0.14075746481359069909215356472191e-10,
            -0.25769072058024680627537078627584e-11,
            +0.47342406666294421849154395005938e-12,
            -0.87249012674742641745301263292675e-13,
            +0.16124614902740551465739833119115e-13,
            -0.29875652015665773006710792416815e-14,
            +0.55480701209082887983041321697279e-15,
            -0.10324619158271569595141333961932e-15,
            +0.19250239203049851177878503244868e-16,
            -0.35955073465265150011189707844266e-17,
            +0.67264542537876857892194574226773e-18,
            -0.12602624168735219252082425637546e-18,
            +0.23644884408606210044916158955519e-19,
            -0.44419377050807936898878389179733e-20,
            +0.83546594464034259016241293994666e-21,
            -0.15731559416479562574899253521066e-21,
            +0.29653128740247422686154369706666e-22,
            -0.55949583481815947292156013226666e-23,
            +0.10566354268835681048187284138666e-23,
            -0.19972483680670204548314999466666e-24,
            +0.37782977818839361421049855999999e-25,
            -0.71531586889081740345038165333333e-26,
            +0.13552488463674213646502024533333e-26,
            -0.25694673048487567430079829333333e-27,
            +0.48747756066216949076459519999999e-28,
            -0.92542112530849715321132373333333e-29,
            +0.17578597841760239233269760000000e-29,
            -0.33410026677731010351377066666666e-30,
            +0.63533936180236187354180266666666e-31
        };

        // R file: chebyshev.c
        internal static int chebyshev_init(double[] dos, int nos, double eta)
        {
            int i, ii;
            double err;

            if (nos < 1)
                return 0;

            err = 0.0;
            i = 0;			/* just to avoid compiler warnings */
            for (ii = 1; ii <= nos; ii++)
            {
                i = nos - ii;
                err += Math.Abs(dos[i]);
                if (err > eta)
                {
                    return i;
                }
            }
            return i;
        }

        // R file: chebyshev.c
        internal static double chebyshev_eval(double x, double[] a, int n)
        {
            double b0, b1, b2, twox;
            int i;

            if (n < 1 || n > 1000) throw new Exception("chebyshev_eval: ML_ERR_return_NAN");

            if (x < -1.1 || x > 1.1) throw new Exception("chebyshev_eval: ML_ERR_return_NAN");

            twox = x * 2;
            b2 = b1 = 0;
            b0 = 0;
            for (i = 1; i <= n; i++)
            {
                b2 = b1;
                b1 = b0;
                b0 = twox * b1 - b2 + a[n - i];
            }
            return (b0 - b2) * 0.5;
        }

        // R file: dpq.h
        //#define R_D_Lval(p)	(lower_tail ? (p) : (0.5 - (p) + 0.5))
        internal static double R_D_Lval(double p, bool lower_tail)
        {
            return (lower_tail ? (p) : (0.5 - (p) + 0.5));
        }

        // R file: dpq.h
        //#define R_DT_Clog(p)	(lower_tail? R_D_LExp(p): R_D_log(p))
        internal static double R_DT_Clog(double p, bool lower_tail, bool log_p)
        {
            return lower_tail ?
                (log_p ?
                (
                    p > -M_LN2 ?
                    Math.Log(-expm1(p)) : log1p(-Math.Exp(p))) :
                    log1p(-p)) :
                (
                    log_p ? p :
                    Math.Log(p));
        }

        // R file: dpq.h
        //#define R_Q_P01_boundaries(p, _LEFT_, _RIGHT_)
        internal static double R_Q_P01_boundaries(double p, double left, double right)
        {

            bool log_p = true, lower_tail = false;
            if (log_p)
            {
                if (p > 0)
                    return double.NaN; // ML_ERR_return_NAN;				
                if (p == 0) /* upper bound*/
                    return lower_tail ? right : left;
                if (p == double.NegativeInfinity)
                    return lower_tail ? left : right;
            }
            else
            { /* !log_p */
                if (p < 0 || p > 1)
                    return double.NaN; // ML_ERR_return_NAN;				
                if (p == 0)
                    return lower_tail ? left : right;
                if (p == 1)
                    return lower_tail ? right : left;
            }
            return double.NaN;
        }

        // R file dpq.h 
        // #define R_D_fexp(f,x)     (give_log ? -0.5*log(f)+(x) : exp(x)/sqrt(f))
        internal static double R_D_fexp(double f, double x, bool give_log)
        {
            return (give_log ? -0.5 * Math.Log(f) + (x) : Math.Exp(x) / Math.Sqrt(f));
        }

        // R file dpq.h  
        // log(1 - exp(x))  in more stable form than log1p(- R_D_qIv(x)) :
        // #define R_Log1_Exp(x)   ((x) > -M_LN2 ? log(-expm1(x)) : log1p(-exp(x)))
        internal static double R_Log1_Exp(double x)
        {
            return ((x) > -M_LN2 ? Math.Log(-expm1(x)) : log1p(-Math.Exp(x)));
        }

        // R file dpq.h  
        // #define R_D_Cval(p)	(lower_tail ? (0.5 - (p) + 0.5) : (p))	/*  1 - p */
        internal static double R_D_Cval(double p, bool lower_tail)
        {
            return (lower_tail ? (0.5 - (p) + 0.5) : (p));	/*  1 - p */
        }

        // R file: dpq.h
        // #define R_DT_CIv(p)
        internal static double R_DT_CIv(double p, bool lower_tail, bool log_p)
        {
            return (log_p ? (lower_tail ? -expm1(p) : Math.Exp(p))
                   : R_D_Cval(p, lower_tail));
        }

        // R file: dpq.h
        // #define R_DT_qIv(p)
        internal static double R_DT_qIv(double p, bool lower_tail, bool log_p)
        {
            return (log_p ? (lower_tail ? Math.Exp(p) : -expm1(p))
               : R_D_Lval(p, lower_tail));
        }

    }
}
