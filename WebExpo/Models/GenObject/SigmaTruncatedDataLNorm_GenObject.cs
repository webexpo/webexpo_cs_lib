namespace Zygotine.WebExpo
{
    using System;
    using System.Linq;
    using System.Collections.Generic;
    using Zygotine.Statistics.Distribution;
    using Zygotine.Util;

    //sigma.truncatedData.gen.lnorm.object de SEG.InformedVar ...
    internal class SigmaTruncatedDataLNorm_GenObject : GenObject
    {

        internal SigmaTruncatedDataLNorm_GenObject(int N, double lNormMu, double lNormSigma2, double m = 0)
        {
            this.A = new SGNFnAParam();
            this.A.N = N;
            this.A.LM = lNormMu;
            this.A.LS2 = lNormSigma2 * lNormSigma2;
            this.A.M = m;
            this.MathLowerLimit = 0;
        }





        internal override double F(double s, SGNFnAParam A)
        {
            return Math.Exp(
            -(A.N + 1) * Math.Log(s)
            - A.B / Math.Pow(s, 2)
            - A.N * NormalDistribution.PNorm(A.Mu / s, log_p: true)
            - (Math.Pow(Math.Log(s) - A.LM, 2)) / (2.0 * A.LS2)
            - A.M);
        }

        internal override double LogF(double s, SGNFnAParam A)
        {
            return
            -(A.N + 1) * Math.Log(s)
            - A.B / Math.Pow(s, 2)
            - A.N * NormalDistribution.PNorm(A.Mu / s, log_p: true)
            - (Math.Pow(Math.Log(s) - A.LM, 2)) / (2.0 * A.LS2);
        }

        internal override double LogFPrime(double s, SGNFnAParam A)
        {
            double z = A.Mu / s;
            LPhiObject l = WebExpoFunctions3.LPhi(z);
            return
                -(A.N + 1.0) / s
                + 2.0 * A.B / Math.Pow(s, 3)
                - (Math.Log(s) - A.LM) / A.LS2 / s
                + A.N * z * l.R[0] / s;
        }

        internal override double LogFSecond(double s, SGNFnAParam A)
        {
            double z = A.Mu / s;
            LPhiObject l = WebExpoFunctions3.LPhi(z);
            return
                (A.N + 1.0) / Math.Pow(s, 2)
                - 6.0 * A.B / Math.Pow(s, 4)
                + (Math.Log(s) - A.LM - 1) / A.LS2 / Math.Pow(s, 2)
                + A.N * z * ((Math.Pow(z, 2) - 2.0) * l.R[0] + z * l.R2[0]) / Math.Pow(s, 2);
        }

        //# new_0.11
        internal override List<double> LogFInvRemote(double target, SGNFnAParam A)
        {
            List<double> roots = new List<double>();
            double c3 = -(A.N + 1.0);
            double c2 = A.N + 1.0 - target;
            double c1 = 0;
            double c0 = -A.B;
            double[] t = WebExpoFunctions3.LPQAC3;
            t = Tools.Combine(-A.N * t[0], -A.N * t[1] * A.Mu, -A.N * t[2] * Math.Pow(A.Mu, 2), -A.N * t[3] * Math.Pow(A.Mu, 3));
            double[] theta = Tools.Combine(c0, c1, c2, c3).Add(t);
            roots.AddRange(Numerics.CubicEquation.GetRealRoots(theta, l: 0, ROrder: true));
            double C = -Math.Pow(A.LM, 2) / 2.0 / A.LS2 + A.N * Math.Log(2);
            double B = A.LM / A.LS2 - (A.N + 1.0);
            double A1 = -1.0 / 2.0 / A.LS2;
            roots.AddRange(Numerics.QuadraticEquation.GetRealRoots(Tools.Combine(C, B, A1), target, ROrder: true).Select(a => Math.Exp(a)));
            return roots;
        } //# end of Math.Log.f.inv.remote

        internal override double[] Start(SGNFnAParam A)
        {
            // Le throw ne devrait as se produire. Il faudrait revoir l'organisation des objects GenObject.  octobre 2018.
            throw new NotImplementedException();
        }

    } //# end of sigma.truncatedData.gen.lnorm.object

}
