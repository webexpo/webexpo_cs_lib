
namespace Zygotine.Statistics.Distribution
{
    using System;
    using Zygotine.Util;

    public static partial class NormalDistribution
    {
        /*
         * R files
         * 
         * snorm.c
         * rnorm
         */
        public static double RNorm()
        {
            return norm_rand();
        }

        public static double RNorm(double mu = 0.0, double sigma = 1.0)
        {
            {
                double retVal;
                if (!TestParameters(mu, sigma, out retVal))
                    return retVal;
            }

            return mu + sigma * norm_rand();
        }

        public static double[] RNorm(int n, double mu = 0.0, double sigma = 1.0)
        {
            {
                double retVal;
                if (!TestParameters(mu, sigma, out retVal))
                    return Tools.GetConstantArray(n, retVal);
            }

            double[] rep = new double[n];
            for (int i = 0; i < n; i++)
                rep[i] = mu + sigma * norm_rand();
            return rep;
        }

        public static double[] RNorm(int n, double[] mu, double[] sigma)
        {
            if (n <= 0)
                return new double[0];
            int lMu = mu.Length, iMu = 0;
            int lSigma = sigma.Length, iSigma = 0;
            double[] rep = new double[n];
            for (int i = 0; i < n; i++)
            {
                rep[i] = RNorm(mu[iMu], sigma[iSigma]);
                iMu = (iMu + 1) % lMu;
                iSigma = (iSigma + 1) % lSigma;
            }
            return rep;
        }

        private static bool TestParameters(double mu, double sigma, out double retValue)
        {
            retValue = double.NaN;

            if (double.IsNaN(mu) || !sigma.IsFinite() || sigma < 0.0)
            {
                retValue = double.NaN; // ML_ERR_return_NAN;
                return false;
            }
            if (sigma == 0.0 || !mu.IsFinite())
            {
                retValue = mu; /* includes mu = +/- Inf with finite sigma */
                return false;
            }
            return true;
        }

        private static double norm_rand()
        {
            // In C, case INVERSION: ...

            double u1;
            const double BIG = 134217728; /* 2^27 */
            /* unif_rand() alone is not of high enough precision */
            u1 = UniformDistribution.RUnif();
            u1 = (int)(BIG * u1) + UniformDistribution.RUnif();
            return QNorm(u1 / BIG); //, 0.0, 1.0, true, false);
        }

        private static double[] a = new double[] {
            2.2352520354606839287,
            161.02823106855587881,
            1067.6894854603709582,
            18154.981253343561249,
            0.065682337918207449113
            };

        private static double[] b = {
            47.20258190468824187,
            976.09855173777669322,
            10260.932208618978205,
            45507.789335026729956
            };
        private static double[] c = {
            0.39894151208813466764,
            8.8831497943883759412,
            93.506656132177855979,
            597.27027639480026226,
            2494.5375852903726711,
            6848.1904505362823326,
            11602.651437647350124,
            9842.7148383839780218,
            1.0765576773720192317e-8
            };
        private static double[] d = {
            22.266688044328115691,
            235.38790178262499861,
            1519.377599407554805,
            6485.558298266760755,
            18615.571640885098091,
            34900.952721145977266,
            38912.003286093271411,
            19685.429676859990727
            };
        private static double[] p = {
            0.21589853405795699,
            0.1274011611602473639,
            0.022235277870649807,
            0.001421619193227893466,
            2.9112874951168792e-5,
            0.02307344176494017303
            };
        private static double[] q = {
            1.28426009614491121,
            0.468238212480865118,
            0.0659881378689285515,
            0.00378239633202758244,
            7.29751555083966205e-5
            };

        internal static int PNorm(object p, bool log_p)
        {
            throw new NotImplementedException();
        }
    }
}
