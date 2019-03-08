namespace Zygotine.Statistics.Distribution
{
    using System;
    using static Zygotine.Statistics.Util.RVaria;
    using Zygotine.Util;

    public static partial class GammaDistribution
    {
        /* R file: rgamma.c
        */


        /* random
         */
        public static double RGammaFromScale(double shape, double scale)
        {
            // shape parameter k (ici a) and a scale parameter θ (ici scale).
            // a est la moyenne lorsque scale = 1, sinon c'est a * scale.
            {
                double retVal;
                if (!TestParameters(shape, scale, out retVal))
                    return retVal;
            }
            return gamma_rand(shape, scale);
        }

        public static double[] RGammaFromScale(int n, double shape, double scale)
        {
            {
                double retVal;
                if (!TestParameters(shape, scale, out retVal))
                    return Zygotine.Util.Tools.GetConstantArray(n, retVal);
            }

            double[] rep = new double[n];
            for (int i = 0; i < n; i++)
                rep[i] = gamma_rand(shape, scale);
            return rep;
        }

        public static double RGammaFromRate(double shape, double rate)
        {
            return RGammaFromScale(shape, 1 / rate);
        }

        public static double[] RGammaFromRate(int n, double shape, double rate)
        {
            return RGammaFromScale(n, shape, 1 / rate);
        }

        private static bool TestParameters(double shape, double scale, out double retValue)
        {
            retValue = double.NaN;
            if (!shape.IsFinite() || !scale.IsFinite() || shape < 0.0 || scale <= 0.0)
            {
                if (scale == 0.0)
                {
                    retValue = 0.0;
                    return false;
                }
                retValue = double.NaN; //ML_ERR_return_NAN;
                return false;
            }
            return true;
        }

        private static double gamma_rand(double shape, double scale)
        {
            // On a et scale sont fini, a>=0 et scale>0
            //if (!R_FINITE(a) || !R_FINITE(scale) || a < 0.0 || scale <= 0.0)
            //{
            //    if (scale == 0.0) return 0.0;
            //    return double.NaN; //ML_ERR_return_NAN;
            //}
            /* State variables [FIXME for threading!] :*/

            double aa = 0.0;
            double aaa = 0.0;
            double s = 0, s2 = 0, d = 0;    /* no. 1 (step 1) */
            double q0 = 0, b = 0, si = 0, c = 0;/* no. 2 (step 4) */
            double e, p, q, r, t, u, v, w, x, ret_val;

            if (shape < 1.0)
            { /* GS algorithm for parameters a < 1 */
                if (shape == 0.0)
                    return 0.0;
                e = 1.0 + exp_m1 * shape;
                for (;;)
                {
                    p = e * UniformDistribution.RUnif();
                    if (p >= 1.0)
                    {
                        x = -Math.Log((e - p) / shape);
                        if (ExponentialDistribution.RExp() >= (1.0 - shape) * Math.Log(x))
                            break;
                    }
                    else
                    {
                        x = Math.Exp(Math.Log(p) / shape);
                        if (ExponentialDistribution.RExp() >= x)
                            break;
                    }
                }
                return scale * x;
            }

            /* --- a >= 1 : GD algorithm --- */

            /* Step 1: Recalculations of s2, s, d if a has changed */
            if (shape != aa)
            {
                aa = shape;
                s2 = shape - 0.5;
                s = Math.Sqrt(s2);
                d = sqrt32 - s * 12.0;
            }
            /* Step 2: t = standard normal deviate,
                       x = (s,1/2) -normal deviate. */

            /* immediate acceptance (i) */
            t = NormalDistribution.RNorm();
            x = s + 0.5 * t;
            ret_val = x * x;
            if (t >= 0.0)
                return scale * ret_val;

            /* Step 3: u = 0,1 - uniform sample. squeeze acceptance (s) */
            u = UniformDistribution.RUnif();
            if (d * u <= t * t * t)
                return scale * ret_val;

            /* Step 4: recalculations of q0, b, si, c if necessary */

            if (shape != aaa)
            {
                aaa = shape;
                r = 1.0 / shape;
                q0 = ((((((q7 * r + q6) * r + q5) * r + q4) * r + q3) * r
                       + q2) * r + q1) * r;

                /* Approximation depending on size of parameter a */
                /* The constants in the expressions for b, si and c */
                /* were established by numerical experiments */

                if (shape <= 3.686)
                {
                    b = 0.463 + s + 0.178 * s2;
                    si = 1.235;
                    c = 0.195 / s - 0.079 + 0.16 * s;
                }
                else if (shape <= 13.022)
                {
                    b = 1.654 + 0.0076 * s2;
                    si = 1.68 / s + 0.275;
                    c = 0.062 / s + 0.024;
                }
                else
                {
                    b = 1.77;
                    si = 0.75;
                    c = 0.1515 / s;
                }
            }
            /* Step 5: no quotient test if x not positive */

            if (x > 0.0)
            {
                /* Step 6: calculation of v and quotient q */
                v = t / (s + s);
                if (Math.Abs(v) <= 0.25)
                    q = q0 + 0.5 * t * t * ((((((a7 * v + a6) * v + a5) * v + a4) * v
                                  + a3) * v + a2) * v + a1) * v;
                else
                    q = q0 - s * t + 0.25 * t * t + (s2 + s2) * Math.Log(1.0 + v);


                /* Step 7: quotient acceptance (q) */
                if (Math.Log(1.0 - u) <= q)
                    return scale * ret_val;
            }
            while (true)
            {
                /* Step 8: e = standard exponential deviate
                 *	u =  0,1 -uniform deviate
                 *	t = (b,si)-double exponential (laplace) sample */
                e = ExponentialDistribution.RExp();
                u = UniformDistribution.RUnif();
                u = u + u - 1.0;
                if (u < 0.0)
                    t = b - si * e;
                else
                    t = b + si * e;
                /* Step	 9:  rejection if t < tau(1) = -0.71874483771719 */
                if (t >= -0.71874483771719)
                {
                    /* Step 10:	 calculation of v and quotient q */
                    v = t / (s + s);
                    if (Math.Abs(v) <= 0.25)
                        q = q0 + 0.5 * t * t *
                            ((((((a7 * v + a6) * v + a5) * v + a4) * v + a3) * v
                              + a2) * v + a1) * v;
                    else
                        q = q0 - s * t + 0.25 * t * t + (s2 + s2) * Math.Log(1.0 + v);
                    /* Step 11:	 hat acceptance (h) */
                    /* (if q not positive go to step 8) */
                    if (q > 0.0)
                    {
                        w = expm1(q);
                        /*  ^^^^^ original code had approximation with rel.err < 2e-7 */
                        /* if t is rejected sample again at step 8 */
                        if (c * Math.Abs(u) <= w * Math.Exp(e - 0.5 * t * t))
                            break;
                    }
                }
            } /* repeat .. until  `t' is accepted */
            x = s + 0.5 * t;
            return scale * x * x;
        }

        private const double
            sqrt32 = 5.656854,
            exp_m1 = 0.36787944117144232159,/* exp(-1) = 1/e */
            q1 = 0.04166669,
            q2 = 0.02083148,
            q3 = 0.00801191,
            q4 = 0.00144121,
            q5 = -7.388e-5,
            q6 = 2.4511e-4,
            q7 = 2.424e-4,

            a1 = 0.3333333,
            a2 = -0.250003,
            a3 = 0.2000062,
            a4 = -0.1662921,
            a5 = 0.1423657,
            a6 = -0.1367177,
            a7 = 0.1233795;


    }
}
