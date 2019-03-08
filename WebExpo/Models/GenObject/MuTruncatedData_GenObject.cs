namespace Zygotine.WebExpo
{
    using System;
    using System.Collections.Generic;
    using Zygotine.Numerics;
    using Zygotine.Statistics.Distribution;
    using Zygotine.Util;

    //mu.truncatedData.gen.object
    internal class MuTruncatedData_GenObject : GenObject
    {
        internal MuTruncatedData_GenObject(int n)
        {
            this.A = new SGNFnAParam();
            this.A.N = n;
            this.A.M = 0;
        }

        internal override double F(double mu, SGNFnAParam A)
        {
            return Math.Exp(-A.N / 2.0 * Math.Pow((mu - A.MuMean) / A.S, 2) - A.N * NormalDistribution.PNorm(mu / A.S, log_p: true) - A.M);
        }

        internal override double LogF(double mu, SGNFnAParam A)
        {
            return -A.N / 2.0 * Math.Pow((mu - A.MuMean) / A.S, 2) - A.N * NormalDistribution.PNorm(mu / A.S, log_p: true);
        }

        internal override double LogFPrime(double mu, SGNFnAParam A)
        {
            double z = mu / A.S;
            LPhiObject lphio = WebExpoFunctions3.LPhi(z);
            return -A.N * (mu - A.MuMean) / A.S2 - A.N * lphio.R[0] / A.S;
        }

        internal override double LogFSecond(double mu, SGNFnAParam A)
        {
            double z = mu / A.S;
            LPhiObject lphio = WebExpoFunctions3.LPhi(z);
            return A.N / A.S2 * (z * lphio.R[0] + lphio.R2[0] - 1);
        }

        internal override List<double> LogFInvRemote(double target, SGNFnAParam A)
        {
            List<double> roots = new List<double>();
            double k = -A.N / 2.0 / A.S2;
            //# right-side roots
            double B = -2.0 * A.MuMean;
            double C = Math.Pow(A.MuMean, 2);
            double[] coeff = Util.Tools.Combine(C, B, 1.0).Multiply(k);
            roots.AddRange(QuadraticEquation.GetRealRoots(coeff, target, ROrder: true));
            //# left-side roots
            double[] lpqac = WebExpoFunctions3.LPQAC2;
            coeff = coeff.Substract(Tools.Combine(A.N * lpqac[0], A.N * lpqac[1] / A.S, A.N * lpqac[2] / A.S2));
            roots.AddRange(QuadraticEquation.GetRealRoots(coeff, target, ROrder: true));
            return roots;
        } //# end of log.f.inv.remote

        internal override double[] Start(SGNFnAParam A)
        {
            // Le throw ne devrait pas se produire. Il faudrait revoir l'organisation des objects GenObject.  octobre 2018.
            throw new NotImplementedException();
        }
    }
}
