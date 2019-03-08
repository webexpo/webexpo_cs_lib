namespace Zygotine.Statistics.Distribution
{
    using System;
    using Zygotine.Statistics.Util;
    using Zygotine.Util;

    public static partial class PoissonDistribution
    {

        public static double RPois(double mu)
        {
            {
                double retVal;
                if (!TestParameters(mu, out retVal))
                    return retVal;
            }
            return poisson_rand(mu);
        }

        public static double[] RPois_N(int n, double mu)
        {
            {
                double retVal;
                if (!TestParameters(mu, out retVal))
                    return Tools.GetConstantArray(n, retVal);
            }
            double[] rep = new double[n];
            for (int i = 0; i < n; i++)
                rep[i] = poisson_rand(mu);
            return rep;
        }

        private static void F(double mu)
        {
        /* Step_F: /* 'subroutine' F : calculation of px,py,fx,fy. */

            if (pois < 10)
            { /* use factorials from table fact[] */
                px = -mu;
                py = Math.Pow(mu, pois) / fact[(int)pois];
            }
            else
            {
                /* Case pois >= 10 uses polynomial approximation
                   a0-a7 for accuracy when advisable */
                del = one_12 / fk;
                del = del * (1.0 - 4.8 * del * del);
                v = difmuk / fk;
                if (Math.Abs(v) <= 0.25)
                    px = fk * v * v * (((((((a7 * v + a6) * v + a5) * v + a4) *
                              v + a3) * v + a2) * v + a1) * v + a0)
                    - del;
                else /* |v| > 1/4 */
                    px = fk * Math.Log(1.0 + v) - difmuk - del;
                py = RVaria.M_1_SQRT_2PI / Math.Sqrt(fk);
            }
            x = (0.5 - difmuk) / s;
            x *= x;/* x^2 */
            fx = -0.5 * x;
            fy = omega * (((c3 * x + c2) * x + c1) * x + c0);
        }

        public static int countDirectToStepF = 0;
        public static int glblCount = 0;
        /* These are static --- persistent between calls for same mu : */
        private static int l, m;
        private static double b1, b2, c, c0, c1, c2, c3;
        private static double[] pp = new double[36];
        private static double p0, p, q, s, d, omega;
        private static double big_l;/* integer "w/o overflow" */
        private static double muprev = 0.0, muprev2 = 0.0;/*, muold	 = 0.*/

        /* Factorial Table (0:9)! */
        private static readonly double[] fact = new double[10]
    {
    1.0, 1.0, 2.0, 6.0, 24.0, 120.0, 720.0, 5040.0, 40320.0, 362880.0
    };
        private const double a0 = -0.5;
        private static readonly double a1 = 0.3333333;
        private static readonly double a2 = -0.2500068;
        private static readonly double a3 = 0.2000118;
        private static readonly double a4 = -0.1661269;
        private static readonly double a5 = 0.1421878;
        private static readonly double a6 = -0.1384794;
        private static readonly double a7 = 0.1250060;
        private static readonly double one_7 = 0.1428571428571428571;
        private static readonly double one_12 = 0.0833333333333333333;
        private static readonly double one_24 = 0.0416666666666666667;

        /* Local Vars  [initialize some for -Wall]: */
        private static double del, difmuk = 0.0, E = 0.0, fk = 0.0, fx, fy, g, px, py, t, u = 0.0, v, x;
        private static double pois = -1.0;
        private static int k, kflag;
        private static bool big_mu, new_big_mu = false;

        internal static bool TestParameters(double mu, out double retValue)
        {
            retValue = double.NaN; //ML_ERR_return_NAN;
            if (!mu.IsFinite() || mu <= 0.0)
            {
                if (mu == 0.0)
                    retValue = 0.0;
                return false;
            }
            return true;
        }

        internal static double poisson_rand(double mu)
        {
            glblCount++;
            /* Local Vars  [initialize some for -Wall]: */
            //double del;
            difmuk = 0.0; E = 0.0; fk = 0.0;
            //double fx, fy, g, px, py, t;
            u = 0.0;
            //double v, x;
            pois = -1.0;
            //int k, kflag;
            //bool big_mu;
            bool directToStepF = false;
            new_big_mu = false;


            big_mu = mu >= 10.0;
            if (big_mu)
                new_big_mu = false;

            if (!(big_mu && mu == muprev))
            {/* maybe compute new persistent par.s */

                if (big_mu)
                {
                    new_big_mu = true;
                    /* Case A. (recalculation of s,d,l	because mu has changed):
                     * The poisson probabilities pk exceed the discrete normal
                     * probabilities fk whenever k >= m(mu).
                     */
                    muprev = mu;
                    s = Math.Sqrt(mu);
                    d = 6.0 * mu * mu;
                    big_l = Math.Floor(mu - 1.1484);
                    /* = an upper bound to m(mu) for all mu >= 10.*/
                }
                else
                { /* Small mu ( < 10) -- not using normal approx. */

                    /* Case B. (start new table and calculate p0 if necessary) */

                    /*muprev = 0.;-* such that next time, mu != muprev ..*/
                    if (mu != muprev)
                    {
                        muprev = mu;
                        m = Math.Max(1, (int)mu);
                        l = 0; /* pp[] is already ok up to pp[l] */
                        q = p0 = p = Math.Exp(-mu);
                    }

                    for (;;)
                    {
                        /* Step U. uniform sample for inversion method */
                        u = UniformDistribution.RUnif();
                        if (u <= p0)
                            return 0.0;

                        /* Step T. table comparison until the end pp[l] of the
                           pp-table of cumulative poisson probabilities
                           (0.458 > ~= pp[9](= 0.45792971447) for mu=10 ) */
                        if (l != 0)
                        {
                            for (k = (u <= 0.458) ? 1 : Math.Min(l, m); k <= l; k++)
                                if (u <= pp[k])
                                    return (double)k;
                            if (l == 35) /* u > pp[35] */
                                continue;
                        }
                        /* Step C. creation of new poisson
                           probabilities p[l..] and their cumulatives q =: pp[k] */
                        l++;
                        for (k = l; k <= 35; k++)
                        {
                            p *= mu / k;
                            q += p;
                            pp[k] = q;
                            if (u <= q)
                            {
                                l = k;
                                return (double)k;
                            }
                        }
                        l = 35;
                    } /* end(repeat) */
                }/* mu < 10 */

            } /* end {initialize persistent vars} */

            /* Only if mu >= 10 : ----------------------- */

            /* Step N. normal sample */
            g = mu + s * NormalDistribution.RNorm();/* norm_rand() ~ N(0,1), standard normal */

            if (g >= 0.0)
            {
                pois = Math.Floor(g);
                /* Step I. immediate acceptance if pois is large enough */
                if (pois >= big_l)
                    return pois;
                /* Step S. squeeze acceptance */
                fk = pois;
                difmuk = mu - fk;
                u = UniformDistribution.RUnif(); /* ~ U(0,1) - sample */
                if (d * u >= difmuk * difmuk * difmuk)
                    return pois;
            }

            /* Step P. preparations for steps Q and H.
               (recalculations of parameters if necessary) */

            if (new_big_mu || mu != muprev2)
            {
                /* Careful! muprev2 is not always == muprev
               because one might have exited in step I or S
               */
                muprev2 = mu;
                omega = RVaria.M_1_SQRT_2PI / s;
                /* The quantities b1, b2, c3, c2, c1, c0 are for the Hermite
                 * approximations to the discrete normal probabilities fk. */

                b1 = one_24 / mu;
                b2 = 0.3 * b1 * b1;
                c3 = one_7 * b1 * b2;
                c2 = b2 - 15.0 * c3;
                c1 = b1 - 6.0 * b2 + 45.0 * c3;
                c0 = 1.0 - b1 + 3.0 * b2 - 15.0 * c3;
                c = 0.1069 / mu; /* guarantees majorization by the 'hat'-function. */
            }

            if (g >= 0.0)
            {
                /* 'Subroutine' F is called (kflag=0 for correct return) */
                kflag = 0;
                //goto Step_F;
                //F(mu);
                directToStepF = true;
                countDirectToStepF++;
            }


            for (;;)
            {
                if (!directToStepF)
                {
                    /* Step E. Exponential Sample */

                    E = ExponentialDistribution.RExp();	/* ~ Exp(1) (standard exponential) */

                    /*  sample t from the laplace 'hat'
                        (if t <= -0.6744 then pk < fk for all mu >= 10.) */
                    u = 2 * UniformDistribution.RUnif() - 1.0;
                    t = 1.8 + RVaria.fsign(E, u);
                }

                if ((t > -0.6744) || directToStepF)
                {
                    if (!directToStepF)
                    {
                        pois = Math.Floor(mu + s * t);
                        fk = pois;
                        difmuk = mu - fk;
                        kflag = 1;
                    }
                    /* 'subroutine' F is called (kflag=1 for correct return) */

                    F(mu);
                    directToStepF = false;
                    //Step_F: /* 'subroutine' F : calculation of px,py,fx,fy. */

                    //    if (pois < 10)
                    //    { /* use factorials from table fact[] */
                    //        px = -mu;
                    //        py = Math.Pow(mu, pois) / fact[(int)pois];
                    //    }
                    //    else
                    //    {
                    //        /* Case pois >= 10 uses polynomial approximation
                    //           a0-a7 for accuracy when advisable */
                    //        del = one_12 / fk;
                    //        del = del * (1.0 - 4.8 * del * del);
                    //        v = difmuk / fk;
                    //        if (Math.Abs(v) <= 0.25)
                    //            px = fk * v * v * (((((((a7 * v + a6) * v + a5) * v + a4) *
                    //                      v + a3) * v + a2) * v + a1) * v + a0)
                    //            - del;
                    //        else /* |v| > 1/4 */
                    //            px = fk * Math.Log(1.0 + v) - difmuk - del;
                    //        py = R.M_1_SQRT_2PI / Math.Sqrt(fk);
                    //    }
                    //    x = (0.5 - difmuk) / s;
                    //    x *= x;/* x^2 */
                    //    fx = -0.5 * x;
                    //    fy = omega * (((c3 * x + c2) * x + c1) * x + c0);
                    if (kflag > 0)
                    {
                        /* Step H. Hat acceptance (E is repeated on rejection) */
                        if (c * Math.Abs(u) <= py * Math.Exp(px + E) - fy * Math.Exp(fx + E))
                            break;
                    }
                    else
                        /* Step Q. Quotient acceptance (rare case) */
                        if (fy - u * fy <= py * Math.Exp(px - fx))
                        break;
                }/* t > -.67.. */
            }
            return pois;
        }
    }
}
