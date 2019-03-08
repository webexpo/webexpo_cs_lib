namespace Zygotine.WebExpo
{
    using System;
    using System.Collections.Generic;
    using Zygotine.Numerics;
    using Zygotine.Statistics.Distribution;
    using Zygotine.Util;

    internal class ME_CV_GenObject : GenObject
    {
        internal ME_CV_GenObject(int N, bool logNormalDistrn)
        {
            this.A = new SGNFnAParam();
            this.A.N = N;
            this.A.M = 0;
            this.LogNormalDistrn = logNormalDistrn;
            this.PotentiallyBimodal = false;
            this.MathLowerLimit = 0;
        }

        internal override double F(double v, SGNFnAParam A)
        {
            return Math.Exp(
                            -A.N * Math.Log(v)
                            - A.B / Math.Pow(v, 2)
                            - A.N * NormalDistribution.PNorm(1.0 / v, log_p: true)
                            - A.M);
        }

        internal override double LogF(double v, SGNFnAParam A)
        {
            return
                -A.N * Math.Log(v)
                - A.B / Math.Pow(v, 2)
                - A.N * NormalDistribution.PNorm(1.0 / v, log_p: true);
        }

        internal override double LogFPrime(double v, SGNFnAParam A)
        {
            LPhiObject l = WebExpoFunctions3.LPhi(1 / v);
            return
                -A.N / v
                + 2.0 * A.B / Math.Pow(v, 3)
                + A.N * l.R[0] / Math.Pow(v, 2);
        }

        internal override double LogFSecond(double v, SGNFnAParam A)
        {
            LPhiObject l = WebExpoFunctions3.LPhi(1 / v);
            double v2 = Math.Pow(v, 2);
            double v4 = Math.Pow(v2, 2);
            return
                A.N / v2
                - 6.0 * A.B / v4
                + A.N / v4 * (l.R[0] * (1 / v - 2.0 * v) + l.R2[0]);
        }

        internal override List<double> LogFInvRemote(double target, SGNFnAParam A)
        {
            List<double> roots = new List<double>();
            double D = -A.B;
            double B = -A.N * (Math.Log(A.Mode) - 1.0) - target;
            double qA = -A.N / A.Mode;

            roots.AddRange(CubicEquation.GetRealRoots(Util.Tools.Combine(D, 0, B, qA), l: 0, u: A.Mode, ROrder: true));// # modif_0.11
            roots.Add(Math.Exp(Math.Log(2.0) - target / A.N));
            return roots;
        }

        internal override double[] Start(SGNFnAParam A)
        {
            double s0 = Math.Sqrt(2.0 * A.B / A.N);
            if ((s0 < A.Range[0]) || (s0 > A.Range[1]))
            {
                s0 = A.Range.Mean();
            }
            return Tools.Combine(s0);
        }
    } //ME_CV_GenObject
}
