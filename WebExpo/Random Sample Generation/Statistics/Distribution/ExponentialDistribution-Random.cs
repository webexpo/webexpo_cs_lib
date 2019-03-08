namespace Zygotine.Statistics.Distribution
{
    using Zygotine.Util;

    /* R files: sexp.c, rexp.c
     */
    public static partial class ExponentialDistribution
    {
        public static double RExp()
        {
            return exp_rand();
        }

        public static double RExpFromScale(double scale)
        {
            {
                double retVal = double.NaN;
                if (!TestParameters(scale, out retVal))
                    return retVal;
            }

            return scale * exp_rand(); // --> in ./sexp.c
        }

        public static double[] RExpFromScale(int n, double scale)
        {

            {
                double retVal = double.NaN;
                if (!TestParameters(scale, out retVal))
                    return Tools.GetConstantArray(n, retVal);
            }


            double[] rep = new double[n];

            for (int i = 0; i < n; i++)
                rep[i] = scale * exp_rand(); // --> in ./sexp.c
            return rep;
        }

        public static double RExpFromRate(double rate)
        {
            return RExpFromScale(1 / rate);
        }

        public static double[] RExpFromRate(int n, double rate)
        {
            return RExpFromScale(n, 1 / rate);
        }


        private static bool TestParameters(double scale, out double retValue)
        {
            retValue = double.NaN;//ML_ERR_return_NAN;
            if (!scale.IsFinite() || scale <= 0.0)
            {
                if (scale == 0.0)
                    retValue = 0.0;
                return false;
            }
            return true;
        }

        private static double[] q =
        {
                0.6931471805599453,
                0.9333736875190459,
                0.9888777961838675,
                0.9984959252914960,
                0.9998292811061389,
                0.9999833164100727,
                0.9999985691438767,
                0.9999998906925558,
                0.9999999924734159,
                0.9999999995283275,
                0.9999999999728814,
                0.9999999999985598,
                0.9999999999999289,
                0.9999999999999968,
                0.9999999999999999,
                1.0000000000000000
                };

        private static double exp_rand()
        {
            /* q[k-1] = sum(log(2)^k / k!)  k=1,..,n, */
            /* The highest n (here 16) is determined by q[n-1] = 1.0 */
            /* within standard precision */


            double a = 0.0;
            double u = UniformDistribution.RUnif();    /* precaution if u = 0 is ever returned */
            // unif_rand based on MersenneTwister returns values in the open interval (0,1).
            // so following line is superflous
            // while (u <= 0.0 || u >= 1.0) u = unif_rand(); 
            for (;;)
            {
                u += u;
                if (u > 1.0)
                    break;
                a += q[0];
            }
            u -= 1.0;

            if (u <= q[0])
                return a + u;

            int i = 0;
            double ustar = UniformDistribution.RUnif(), umin = ustar;
            do
            {
                ustar = UniformDistribution.RUnif();
                if (umin > ustar)
                    umin = ustar;
                i++;
            } while (u > q[i]);
            return a + umin * q[0];
        }



    }

}

