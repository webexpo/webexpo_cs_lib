namespace Zygotine.WebExpo
{
    using System;
    using System.Linq;
    using System.Collections.Generic;
    using Zygotine.Numerics;
    using Zygotine.Statistics.Distribution;
    using Zygotine.Util;

    internal class Sigma_woLM_GenObject : GenObject
    {
        internal Sigma_woLM_GenObject(int N)
        {
            //de this.A  seront ignorés lNormMu et lNormSd
            this.A = new SGNFnAParam();
            this.A.N = N;
            this.A.M = 0;
        }

        internal override double F(double sigma, SGNFnAParam A)
        {
            return Math.Exp(-A.N * Math.Log(sigma) - A.B / Math.Pow(sigma, 2) - A.M);
        }

        internal override double LogF(double sigma, SGNFnAParam A)
        {
            return -A.N * Math.Log(sigma) - A.B / Math.Pow(sigma, 2);
        }

        internal override double LogFPrime(double sigma, SGNFnAParam A)
        {
            return -A.N / sigma + 2.0 * A.B / Math.Pow(sigma, 3);
        }

        internal override double LogFSecond(double sigma, SGNFnAParam A)
        {
            double s2 = Math.Pow(sigma, 2);
            double s4 = Math.Pow(s2, 2);
            return A.N / s2 - 6.0 * A.B / s4;
        }

        internal override List<double> LogFInvRemote(double target, SGNFnAParam A)
        {
            List<double> roots = new List<double>();
            double c3 = -A.N;
            double c2 = A.N - target;
            double c0 = -A.B;
            roots.AddRange(CubicEquation.GetRealRoots(Tools.Combine(c0, 0, c2, c3), l: 0));
            roots.Add(-target / A.N);
            return roots;
        } //# end of log.f.inv.remote

        internal override double[] Start(SGNFnAParam A)
        {
            return Tools.Combine(Math.Sqrt(2.0 * A.B / A.N));
        }
    } //Sigma_woLM_GenObject
}
