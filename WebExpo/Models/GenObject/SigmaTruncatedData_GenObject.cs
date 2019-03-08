namespace Zygotine.WebExpo
{
    using System;
    using System.Collections.Generic;
    using Zygotine.Statistics.Distribution;
    using Zygotine.Util;

    //sigma.truncatedData.gen.object
    internal class SigmaTruncatedData_GenObject : GenObject
    {
        internal SigmaTruncatedData_GenObject(int n)
        {
            this.A = new SGNFnAParam();
            this.A.N = n;
            this.A.M = 0;
        }

        internal override double F(double s, SGNFnAParam A)
        {
            return Math.Exp(-A.N * Math.Log(s) - A.B / Math.Pow(s, 2) - A.N * NormalDistribution.PNorm(A.Mu / s, log_p: true) - A.M);
        }

        internal override double LogF(double s, SGNFnAParam A)
        {
            return -A.N * Math.Log(s) - A.B / Math.Pow(s, 2) - A.N * NormalDistribution.PNorm(A.Mu / s, log_p: true);
        }

        internal override double LogFPrime(double s, SGNFnAParam A)
        {
            double z = A.Mu / s;
            LPhiObject l = WebExpoFunctions3.LPhi(z);
            return -A.N / s + 2.0 * A.B / Math.Pow(s, 3) + A.N * z * l.R[0] / s;
        }

        internal override double LogFSecond(double s, SGNFnAParam A)
        {
            double z = A.Mu / s;
            LPhiObject l = WebExpoFunctions3.LPhi(z);
            return A.N / Math.Pow(s, 2) - 6.0 * A.B / Math.Pow(s, 4) + A.N * z * ((Math.Pow(z, 2) - 2.0) * l.R[0] + z * l.R2[0]) / Math.Pow(s, 2);
        }

        internal override List<double> LogFInvRemote(double target, SGNFnAParam A)
        {
            double c3 = -A.N;
            double c2 = A.N - target;
            double c0 = -A.B;
            List<double> roots = new List<double>();
            double[] t = WebExpoFunctions3.LPQAC3;
            t = t.Multiply(Tools.Combine(1, A.Mu, Math.Pow(A.Mu, 2), Math.Pow(A.Mu, 3))).Multiply(-A.N);
            double[] theta = Tools.Combine(c0, 0, c2, c3).Add(t);
            roots.AddRange(Numerics.CubicEquation.GetRealRoots(theta, l: 0));
            roots.Add(2 * Math.Exp(-target / A.N));
            return roots;
        }


        internal override double[] Start(SGNFnAParam A)
        {
            double s0 = Math.Sqrt(2 * A.B / A.N);
            double g0 = NormalDistribution.PNorm(A.Mu / s0, log_p: true);
            double gp0 = -NormalDistribution.DNorm(A.Mu / s0) * A.Mu / Math.Pow(s0, 2) / NormalDistribution.PNorm(A.Mu / s0);
            double c3 = A.N * (gp0 - 1);
            double c2 = A.N * (1 - g0 + gp0 * s0);
            double c0 = -A.B;
            double[] theta = Tools.Combine(c0, 0, c2, c3);
            double[] roots = Numerics.CubicEquation.GetRealRoots(theta, l: 0).ToArray();
            if (roots.Length == 0)
            {
                roots = Tools.Combine(s0);
            }

            return roots;
        }
    }
}
