namespace Zygotine.Statistics.Distribution
{
    using Zygotine.Util;

    public static partial class UniformDistribution
    {

        /*
         * 
         * We do not use R unif_rand
         * Instead we refer to MersenneTwister genrand_real3
         * 
         */

        //internal static RVaria.FN_Dbl_void RUnif = RNG.mtPRNG.genrand_real3;// .genrand_real3;

        //public static double RUnif()
        //{
        //    return unif_rand();
        //}

        public static double RUnif()
        {
            return RNG.MTPrng.genrand_real3();
        }

        public static double RUnif(double a, double b)
        {
            if (!a.IsFinite() || !b.IsFinite() || b < a)
                return double.NaN;

            if (a == b)
                return a;
            else
            {
                return a + (b - a) * RUnif();
            }
        }

        public static double[] RUnif(int n)
        {
            double[] rep = new double[n];
            for (int i = 0; i < n; i++)
                rep[i] = RUnif();
            return rep;
        }

        public static double[] RUnif53(int n)
        {
            double[] rep = new double[n];
            for (int i = 0; i < n; i++)
            {
                rep[i] = RNG.MTPrng.genrand_res53();
            }

            return rep;
        }

        public static double[] RUnif(int n, double[] a, double[] b)
        {
            double[] rep;
            rep = RUnif(n);
            int ia = 0, ib = 0, irep = 0;

            for (irep = 0; irep < n; irep++)
            {
                if (!a[ia].IsFinite() || !b[ib].IsFinite() || b[ib] < a[ia])
                    rep[irep] = double.NaN;
                else
                    rep[irep] = a[ia] + (b[ib] - a[ia]) * rep[irep];

                ia++;
                if (ia == a.Length) ia = 0;
                ib++;
                if (ib == b.Length) ib = 0;
            }
            return rep;
        }

        public static double[] RUnif(int n, double a, double b)
        {
            return RUnif(n, new double[] { a }, new double[] { b });
        }


    }
}
