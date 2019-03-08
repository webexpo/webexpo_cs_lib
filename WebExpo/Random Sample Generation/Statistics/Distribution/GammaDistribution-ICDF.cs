using System;
using Zygotine.Numerics;
using static Zygotine.Statistics.Util.RVaria;
using Zygotine.Util;

namespace Zygotine.Statistics.Distribution
{
    public static partial class GammaDistribution
    {
        /* 
         * R file: qgamma.c
         */

        /* 
         * Icdf 
         */
        public static double QGammaFromScale(double p, double shape, double scale, bool lowerTail = true, bool logP = false)
        {
            double
                i420 = 1.0 / 420.0,
                i2520 = 1.0 / 2520.0,
                i5040 = 1.0 / 5040.0;

            double p_, a, b, c, g, ch, ch0, p1;
            double p2, q, s1, s2, s3, s4, s5, s6, t, x;
            int i, max_it_Newton = 1;

            /* test arguments and initialise */

            if (double.IsNaN(p) || double.IsNaN(shape) || double.IsNaN(scale))
                return p + shape + scale;
            //R_Q_P01_boundaries(p, 0.0, ML_POSINF);
            if (logP)
            {
                if (p > 0)
                    return double.NaN; // ML_ERR_return_NAN;				
                if (p == 0) /* upper bound*/
                    return lowerTail ? double.NegativeInfinity : 0.0;
                if (p == double.NegativeInfinity)
                    return lowerTail ? 0.0 : double.NegativeInfinity;
            }
            else
            { /* !log_p */
                if (p < 0 || p > 1)
                    return double.NaN; // ML_ERR_return_NAN;				
                if (p == 0)
                    return lowerTail ? 0.0 : double.NegativeInfinity;
                if (p == 1)
                    return lowerTail ? double.NegativeInfinity : 0.0;
            }

            if (shape < 0 || scale <= 0)
                return double.NaN; //ML_ERR_return_NAN;

            if (shape == 0) /* all mass at 0 : */
                return 0.0;

            if (shape < 1e-10)
            {
                /* Warning seems unnecessary now: */
                //MATHLIB_WARNING("value of shape (%g) is extremely small: results may be unreliable", alpha);
                max_it_Newton = 7;/* may still be increased below */
            }

            p_ = R_DT_qIv(p, lowerTail, logP);/* lower_tail prob (in any case) */

            g = Functions.LogGamma(shape);/* log Gamma(v/2) */

            /*----- Phase I : Starting Approximation */
            ch = qchisq_appr(p, /* nu= 'df' =  */ 2 * shape, /* lgamma(nu/2)= */ g,
                     lowerTail, logP, /* tol= */ EPS1);
            if (!ch.IsFinite())
            {
                /* forget about all iterations! */
                max_it_Newton = 0; goto END;
            }
            if (ch < EPS2)
            {/* Corrected according to AS 91; MM, May 25, 1999 */
                max_it_Newton = 20;
                goto END;/* and do Newton steps */
            }

            /* FIXME: This (cutoff to {0, +Inf}) is far from optimal
             * -----  when log_p or !lower_tail, but NOT doing it can be even worse */
            if (p_ > pMAX || p_ < pMIN)
            {
                /* did return ML_POSINF or 0.;	much better: */
                max_it_Newton = 20;
                goto END;/* and do Newton steps */
            }


            /*----- Phase II: Iteration
             *	Call pgamma() [AS 239]	and calculate seven term taylor series
             */
            c = shape - 1;
            s6 = (120 + c * (346 + 127 * c)) * i5040; /* used below, is "const" */

            ch0 = ch;/* save initial approx. */
            for (i = 1; i <= MAXIT; i++)
            {
                q = ch;
                p1 = 0.5 * ch;
                p2 = p_ - pgamma_raw(p1, shape, /*lower_tail*/true, /*log_p*/false);
                if (!p2.IsFinite() || ch <= 0)
                { ch = ch0; max_it_Newton = 27; goto END; }/*was  return ML_NAN;*/

                t = p2 * Math.Exp(shape * M_LN2 + g + p1 - c * Math.Log(ch));
                b = t / ch;
                a = 0.5 * t - b * c;
                s1 = (210 + a * (140 + a * (105 + a * (84 + a * (70 + 60 * a))))) * i420;
                s2 = (420 + a * (735 + a * (966 + a * (1141 + 1278 * a)))) * i2520;
                s3 = (210 + a * (462 + a * (707 + 932 * a))) * i2520;
                s4 = (252 + a * (672 + 1182 * a) + c * (294 + a * (889 + 1740 * a))) * i5040;
                s5 = (84 + 2264 * a + c * (1175 + 606 * a)) * i2520;

                ch += t * (1 + 0.5 * t * s1 - b * c * (s1 - b * (s2 - b * (s3 - b * (s4 - b * (s5 - b * s6))))));
                if (Math.Abs(q - ch) < EPS2 * ch)
                    goto END;
                if (Math.Abs(q - ch) > 0.1 * ch)
                {/* diverging? -- also forces ch > 0 */
                    if (ch < q) ch = 0.9 * q; else ch = 1.1 * q;
                }
            }
            /* no convergence in MAXIT iterations -- but we add Newton now... */
            END:
            /* PR# 2214 :	 From: Morten Welinder <terra@diku.dk>, Fri, 25 Oct 2002 16:50
               --------	 To: R-bugs@biostat.ku.dk     Subject: qgamma precision

               * With a final Newton step, double accuracy, e.g. for (p= 7e-4; nu= 0.9)
               *
               * Improved (MM): - only if rel.Err > EPS_N (= 1e-15);
               *		    - also for lower_tail = FALSE	 or log_p = TRUE
               * 		    - optionally *iterate* Newton
               */
            x = 0.5 * scale * ch;
            if (max_it_Newton > 0)
            {
                /* always use log scale */
                if (!logP)
                {
                    p = Math.Log(p);
                    logP = true;
                }
                if (x == 0)
                {
                    const double _1_p = 1.0 + 1e-7;
                    const double _1_m = 1.0 - 1e-7;
                    x = DBL_MIN;
                    p_ = PGammaFromScaleParameter(x, shape, scale, lowerTail, logP);
                    if ((lowerTail && p_ > p * _1_p) ||
                       (!lowerTail && p_ < p * _1_m))
                        return (0.0);
                    /* else:  continue, using x = DBL_MIN instead of  0  */
                }
                else
                    p_ = PGammaFromScaleParameter(x, shape, scale, lowerTail, logP);
                if (p_ == double.NegativeInfinity) return 0; /* PR#14710 */
                for (i = 1; i <= max_it_Newton; i++)
                {
                    p1 = p_ - p;
                    if (Math.Abs(p1) < Math.Abs(EPS_N * p))
                        break;
                    /* else */
                    if ((g = DGammaFromScale(x, shape, scale, logP)) == (logP ? double.NegativeInfinity : 0.0))
                    {
                        break;
                    }
                    /* else :
                     * delta x = f(x)/f'(x);
                     * if(log_p) f(x) := log P(x) - p; f'(x) = d/dx log P(x) = P' / P
                     * ==> f(x)/f'(x) = f*P / P' = f*exp(p_) / P' (since p_ = log P(x))
                     */
                    t = logP ? p1 * Math.Exp(p_ - g) : p1 / g;/* = "delta x" */
                    t = lowerTail ? x - t : x + t;
                    p_ = PGammaFromScaleParameter(t, shape, scale, lowerTail, logP);
                    if (Math.Abs(p_ - p) > Math.Abs(p1) ||
                    (i > 1 && Math.Abs(p_ - p) == Math.Abs(p1)) /* <- against flip-flop */)
                    {
                        /* no improvement */
                        break;
                    } /* else : */
                      //# ifdef Harmful_notably_if_max_it_Newton_is_1
                      /* control step length: this could have started at
                         the initial approximation */
                    if (t > 1.1 * x) t = 1.1 * x;
                    else if (t < 0.9 * x) t = 0.9 * x;
                    //#endif
                    x = t;
                }
            }

            return x;
        }

        public static double[] QGammaFromScale(double[] p, double shape, double scale, bool lowerTail = true, bool logP = false)
        {
            Func<double, double> CloseQGamma = delegate (double p1) { return QGammaFromScale(p1, shape, scale, lowerTail: lowerTail, logP: logP); };

            double[] rep = new double[p.Length];
            for (int i = 0; i < rep.Length; i++)
            {
                rep[i] = CloseQGamma(p[i]);
            }
            return rep;
        }

        public static double QGammaFromRate(double p, double shape, double rate, bool lowerTail = true, bool logP = false)
        {
            if (rate == 0)
            {
                return double.NaN;
            }

            double scale = 1 / rate;
            return QGammaFromScale(p, shape, scale, lowerTail, logP);
        }

        public static double[] QGammaFromRate(double[] p, double shape, double rate, bool lowerTail = true, bool logP = false)
        {
            if (rate == 0)
            {
                return Tools.Rep(double.NaN, p.Length);
            }

            double scale = 1 / rate;
            return QGammaFromScale(p, shape, scale, lowerTail, logP);
        }

        private const double
           C7 = 4.67,
           C8 = 6.66,
           C9 = 6.73,
           C10 = 13.32;

        private static double qchisq_appr(
                double p,
                double nu,
                double g /* = log Gamma(nu/2) */,
                bool lower_tail,
                bool log_p,
                double tol /* EPS1 */)
        {

            double alpha, a, c, ch, p1;
            double p2, q, t, x;

            /* test arguments and initialise */

            if (double.IsNaN(p) || double.IsNaN(nu))
                return p + nu;

            if ((log_p && (p > 0)) ||
                (!log_p && (p < 0 || p > 1)))
                return double.NaN; //ML_ERR_return_NAN

            //R_Q_P01_check(p);

            if (nu <= 0)
                return double.NaN; //ML_ERR_return_NAN;

            alpha = 0.5 * nu;/* = [pq]gamma() shape */
            c = alpha - 1;
            //p1 = (lower_tail ? (log_p ? (p) : Math.Log(p)) : ((p) > -M_LN2 ? Math.Log(-expm1(p)) : log1p(-Math.Exp(p))));
            p1 = (lower_tail ? (log_p ? (p) : Math.Log(p)) : (log_p ? ((p) > -M_LN2 ? Math.Log(-expm1(p)) : log1p(-Math.Exp(p))) : log1p(-p)));
            if (nu < (-1.24 * p1))
            {   /* for small chi-squared */
                /* Math.Log(alpha) + g = Math.Log(alpha) + Math.Log(gamma(alpha)) =
                 *        = Math.Log(alpha*gamma(alpha)) = lgamma(alpha+1) suffers from
                 *  catastrophic cancellation when alpha << 1
                 */
                double lgam1pa = (alpha < 0.5) ? GammaDistribution.lgamma1p(alpha) : (Math.Log(alpha) + g);
                ch = Math.Exp((lgam1pa + p1) / alpha + M_LN2);

            }
            else if (nu > 0.32)
            {   /*  using Wilson and Hilferty estimate */

                x = NormalDistribution.QNorm(p, 0, 1, lower_tail, log_p);
                p1 = 2.0 / (9 * nu);
                ch = nu * Math.Pow(x * Math.Sqrt(p1) + 1 - p1, 3);

                /* approximation for p tending to 1: */
                if (ch > 2.2 * nu + 6)
                    ch = -2 * (R_DT_Clog(p, lower_tail, log_p) - c * Math.Log(0.5 * ch) + g);
            }
            else
            { /* "small nu" : 1.24*(-Math.Math.Log(p)) <= nu <= 0.32 */
                ch = 0.4;
                a = R_DT_Clog(p, lower_tail, log_p) + g + c * M_LN2;

                do
                {
                    q = ch;
                    p1 = 1.0 / (1 + ch * (C7 + ch));
                    p2 = ch * (C9 + ch * (C8 + ch));
                    t = -0.5 + (C7 + 2 * ch) * p1 - (C9 + ch * (C10 + 3 * ch)) / p2;
                    ch -= (1 - Math.Exp(a + 0.5 * ch) * p2 * p1) / t;
                } while (Math.Abs(q - ch) > tol * Math.Abs(ch));
            }

            return ch;
        }

        private static double
            EPS1 = 1e-2,
            EPS2 = 5e-7,/* final precision of AS 91 */
            EPS_N = 1e-15,/* precision of Newton step / iterations */
            MAXIT = 1000,/* was 20 */
            pMIN = 1e-100,  /* was 0.000002 = 2e-6 */
            pMAX = (1 - 1e-14);/* was (1-1e-12) and 0.999998 = 1 - 2e-6 */


    }

}
