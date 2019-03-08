namespace Zygotine.Statistics.Distribution
{
    using System;
    using Zygotine.Statistics.Util;
    using Zygotine.Util;

    public static partial class NormalDistribution
    {
        /* R file: pnorm.c
         */

        /* Cdf 
         */

        public static double PNorm(double x, double mu = 0.0, double sigma = 1.0, bool lower_tail = true, bool log_p = false)
        {
            // R function name is pnorm5
            double p, cp = 0.0;

            /* Note: The structure of these checks has been carefully thought through.
             * For example, if x == mu and sigma == 0, we get the correct answer 1.
             */
            if (double.IsNaN(x) || double.IsNaN(mu) || double.IsNaN(sigma))
                return x + mu + sigma;

            if (!x.IsFinite() && mu == x) return double.NaN;/* x-mu is NaN */
            if (sigma <= 0)
            {
                if (sigma < 0) return double.NaN;
                /* sigma = 0 : */
                return (x < mu) ? (lower_tail ? (log_p ? double.NegativeInfinity : 0.0) : (log_p ? 0.0 : 1.0)) : (lower_tail ? (log_p ? 0.0 : 1.0) : (log_p ? Double.NegativeInfinity : 0.0));
            }
            p = (x - mu) / sigma;
            if (!p.IsFinite())
                return (x < mu) ? (lower_tail ? (log_p ? double.NegativeInfinity : 0.0) : (log_p ? 0.0 : 1.0)) : (lower_tail ? (log_p ? 0.0 : 1.0) : (log_p ? double.NegativeInfinity : 0.0));
            x = p;

            pnorm_both(x, ref p, ref cp, !lower_tail, log_p);

            return (lower_tail ? p : cp);
        }

        public static double[] PNorm(double[] x, double mu = 0.0, double sigma = 1.0, bool lower_tail = true, bool log_p = false)
        {
            double[] rep = new double[x.Length];
            double p, cp = 0.0;

            for (int i1 = 0; i1 < x.Length; i1++)
            {
                rep[i1] = double.NaN;
            }

            if (double.IsNaN(mu) || double.IsNaN(sigma))
            {
                return rep;
            }

            for (int i2 = 0; i2 < x.Length; i2++)
            {
                if (double.IsNaN(x[i2]) || (!x[i2].IsFinite() && mu == x[i2]))
                {
                    rep[i2] = double.NaN;
                    continue;
                }

                if (sigma <= 0)
                {
                    if (sigma < 0)
                    {
                        rep[i2] = double.NaN;
                        break;
                    }
                    /* sigma = 0 : */
                    rep[i2] = (x[i2] < mu) ? (lower_tail ? (log_p ? double.NegativeInfinity : 0.0) : (log_p ? 0.0 : 1.0)) : (lower_tail ? (log_p ? 0.0 : 1.0) : (log_p ? Double.NegativeInfinity : 0.0));
                    continue;
                }

                p = (x[i2] - mu) / sigma;
                if (!p.IsFinite())
                {
                    rep[i2] = (x[i2] < mu) ? (lower_tail ? (log_p ? double.NegativeInfinity : 0.0) : (log_p ? 0.0 : 1.0)) : (lower_tail ? (log_p ? 0.0 : 1.0) : (log_p ? double.NegativeInfinity : 0.0));
                    break;
                }
                double x2 = p;
                pnorm_both(x2, ref p, ref cp, !lower_tail, log_p);
                rep[i2] = (lower_tail ? p : cp);

            }

            return rep;
        }

        public static double[] PNorm(double[] x, double[] mu, double[] sigma, bool lower_tail = true, bool log_p = false)
        {
            int n = x.Length > mu.Length ? x.Length : mu.Length;
            if (n < sigma.Length)
                n = sigma.Length;

            double[] rep = new double[n];

            int iX, iMu, iSigma;
            for (int i = 0; i < n; i++)
            {
                iX = i % x.Length;
                iMu = i % mu.Length;
                iSigma = i % sigma.Length;
                rep[i] = PNorm(x[iX], mu[iMu], sigma[iSigma], lower_tail, log_p);
            }
            return rep;
        }

        private static void pnorm_both(double x, ref double cum, ref double ccum, bool i_tail, bool log_p)
        {
            /* i_tail in {0,1,2} means: "lower", "upper", or "both" :
               if(lower) return  *cum := P[X <= x]
               if(upper) return *ccum := P[X >  x] = 1 - P[X <= x]
            */

            double xden, xnum, temp, del, eps, xsq, y;
            double min = RVaria.DBL_MIN;
            
            int i;
            bool lower, upper;

            if (double.IsNaN(x)) { cum = ccum = x; return; }

            /* Consider changing these : */
            eps = RVaria.DBL_EPSILON * 0.5;

            /* i_tail in {0,1,2} =^= {lower, upper, both} */
            lower = !i_tail; //  = i_tail != 1;
            upper = i_tail; // = i_tail != 0

            y = Math.Abs(x);
            if (y <= 0.67448975)
            { /* qnorm(3/4) = .6744.... -- earlier had 0.66291 */
                if (y > eps)
                {
                    xsq = x * x;
                    xnum = a[4] * xsq;
                    xden = xsq;
                    for (i = 0; i < 3; ++i)
                    {
                        xnum = (xnum + a[i]) * xsq;
                        xden = (xden + b[i]) * xsq;
                    }
                }
                else xnum = xden = 0.0;

                temp = x * (xnum + a[3]) / (xden + b[3]);
                if (lower) cum = 0.5 + temp;
                if (upper) ccum = 0.5 - temp;
                if (log_p)
                {
                    if (lower) cum = Math.Log(cum);
                    if (upper) ccum = Math.Log(ccum);
                }
            }
            else if (y <= RVaria.M_SQRT_32)
            {

                /* Evaluate pnorm for 0.674.. = qnorm(3/4) < |x| <= sqrt(32) ~= 5.657 */

                xnum = c[8] * y;
                xden = y;
                for (i = 0; i < 7; ++i)
                {
                    xnum = (xnum + c[i]) * y;
                    xden = (xden + d[i]) * y;
                }
                temp = (xnum + c[7]) / (xden + d[7]);

                xsq = Math.Truncate(y * SIXTEN) / SIXTEN;
                del = (y - xsq) * (y + xsq);
                if (log_p)
                {
                    cum = (-xsq * xsq * 0.5) + (-del * 0.5) + Math.Log(temp);
                    if ((lower && x > 0.0) || (upper && x <= 0.0))
                        ccum = RVaria.log1p(-Math.Exp(-xsq * xsq * 0.5) *
                              Math.Exp(-del * 0.5) * temp);
                }
                else
                {
                    cum = Math.Exp(-xsq * xsq * 0.5) * Math.Exp(-del * 0.5) * temp;
                    ccum = 1.0 - cum;
                }
                if (x > 0.0)
                {/* swap  ccum <--> cum */
                    temp = cum;
                    if (lower) cum = ccum; ccum = temp;
                }
            }

            else if ((log_p && y < 1e170) /* avoid underflow below */
                                          /*  ^^^^^ MM FIXME: can speedup for log_p and much larger |x| !
                                           * Then, make use of  Abramowitz & Stegun, 26.2.13, something like

                                           xsq = x*x;

                                           if(xsq * DBL_EPSILON < 1.)
                                              del = (1. - (1. - 5./(xsq+6.)) / (xsq+4.)) / (xsq+2.);
                                           else
                                              del = 0.;
                                           *cum = -.5*xsq - M_LN_SQRT_2PI - log(x) + log1p(-del);
                                           *ccum = log1p(-exp(*cum)); /.* ~ log(1) = 0 *./

                                           swap_tail;

                                           [Yes, but xsq might be infinite.]

                                          */
                || (lower && -37.5193 < x && x < 8.2924)
                || (upper && -8.2924 < x && x < 37.5193)
            )
            {

                /* Evaluate pnorm for x in (-37.5, -5.657) union (5.657, 37.5) */
                xsq = 1.0 / (x * x); /* (1./x)*(1./x) might be better */
                xnum = p[5] * xsq;
                xden = xsq;
                for (i = 0; i < 4; ++i)
                {
                    xnum = (xnum + p[i]) * xsq;
                    xden = (xden + q[i]) * xsq;
                }
                temp = xsq * (xnum + p[4]) / (xden + q[4]);
                temp = (RVaria.M_1_SQRT_2PI - temp) / y;

                // do_del(x);

                xsq = Math.Truncate(x * SIXTEN) / SIXTEN;
                del = (x - xsq) * (x + xsq);
                if (log_p)
                {
                    cum = (-xsq * xsq * 0.5) + (-del * 0.5) + Math.Log(temp);
                    if ((lower && x > 0.0) || (upper && x <= 0.0))
                        ccum = RVaria.log1p(-Math.Exp(-xsq * xsq * 0.5) *
                              Math.Exp(-del * 0.5) * temp);
                }
                else
                {
                    cum = Math.Exp(-xsq * xsq * 0.5) * Math.Exp(-del * 0.5) * temp;
                    ccum = 1.0 - cum;
                }
                if (x > 0.0)
                {/* swap  ccum <--> cum */
                    temp = cum;
                    if (lower) cum = ccum; ccum = temp;
                }
            }
            else
            { /* large x such that probs are 0 or 1 */
                if (x > 0) { cum = (log_p ? 0.0 : 1.0); ccum = (log_p ? double.NegativeInfinity : 0.0); }
                else { cum = (log_p ? double.NegativeInfinity : 0.0); ccum = (log_p ? 0.0 : 1.0); }
            }


            //#ifdef NO_DENORMS
            /* do not return "denormalized" -- we do in R */
            if (log_p)
            {
                if (cum > -min) cum = -0.0;
                if (ccum > -min) ccum = -0.0;
            }
            else
            {
                if (cum < min)
                    cum = 0.0;
                if (ccum < min)
                    ccum = 0.0;
            }
            //#endif
            return;
        }
    }
}
