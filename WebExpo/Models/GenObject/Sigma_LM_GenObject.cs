namespace Zygotine.WebExpo
{
    using System;
    using System.Collections.Generic;
    using Zygotine.Numerics;
    using Zygotine.Util;

    internal class Sigma_LM_GenObject : GenObject
    {
        internal Sigma_LM_GenObject(int N, double lNormMu, double lNormSd)
        {
            this.A = new SGNFnAParam();
            this.A.N = N;
            this.A.M = 0;
            this.A.LM = lNormMu;
            this.A.LS2 = lNormSd * lNormSd;
            this.PotentiallyBimodal = false;
            this.MathLowerLimit = 0.0;
        }

        internal override double F(double sigma, SGNFnAParam A)
        {
            //          exp(-(A$N+1)*log(sigma) - A$b/sigma^2 - (log(sigma)-A$lm)^2/(2*A$ls2) - A$M)}
            return Math.Exp(-(A.N + 1.0) * Math.Log(sigma) - A.B / Math.Pow(sigma, 2) - Math.Pow(Math.Log(sigma) - A.LM, 2) / (2.0 * A.LS2) - A.M);
        }

        internal override double LogF(double sigma, SGNFnAParam A)
        {
            return -(A.N + 1.0) * Math.Log(sigma) - A.B / Math.Pow(sigma, 2) - Math.Pow(Math.Log(sigma) - A.LM, 2) / (2.0 * A.LS2);
        } //

        internal override double LogFPrime(double sigma, SGNFnAParam A)
        {
            return (-(A.N + 1.0) / sigma) + (2.0 * A.B / Math.Pow(sigma, 3)) - ((Math.Log(sigma) - A.LM) / (sigma * A.LS2));
        }//LogFPrime

        internal override double LogFSecond(double sigma, SGNFnAParam A)
        {
            //(A$N + 1)/ sigma ^ 2 - 6 * A$b / sigma ^ 4 + (log(sigma) - A$lm - 1)/ (A$ls2* sigma^ 2)
            double s2 = Math.Pow(sigma, 2);
            double s4 = Math.Pow(s2, 2);
            double x = (A.N + 1.0) / s2 - (6.0 * A.B / s4) + ((Math.Log(sigma) - A.LM - 1.0) / (A.LS2 * s2));
            return x;
        } //LogFSecond

        internal override List<double> LogFInvRemote(double target, SGNFnAParam A)
        {
            List<double> roots = new List<double>();
            //# small values

            double c3 = -(A.N + 1.0) + A.LM / A.LS2;
            double c2 = A.N + 1 - Math.Pow(A.LM, 2) / 2.0 / A.LS2 - target;
            double c0 = -A.B;
            roots.AddRange(CubicEquation.GetRealRoots(Tools.Combine(c3, c2, 0, c0), l: 0));
            //# large values
            double C = -Math.Pow(A.LM, 2) / 2.0 / A.LS2 - target;
            double B = A.LM / A.LS2 - (A.N + 1);
            double vA = -1.0 / 2.0 / A.LS2;
            double[] tmp = QuadraticEquation.GetRealRoots(Tools.Combine(C, B, vA), l: 0, ROrder: true).ToArray().Exp();
            roots.AddRange(tmp);
            return roots;
        } //LogFInvRemote

        internal override double[] Start(SGNFnAParam A)
        {
            double s0 = Math.Exp(A.LM);
            double s1;
            double tmp = 2.0 * A.B / (A.N + 1.0 - A.LM / A.LS2);
            if (tmp > 0)
            {
                s1 = Math.Sqrt(tmp);
                return Tools.Combine(s0, s1);
            }
            else
            {
                return Tools.Combine(s0);
            }
        } //Start

        internal override string AAsString(SGNFnAParam A)
        {
            return string.Format("B={0}, LM={1}, LS2={2}, M={3}, N={4}, LowerLimit={5}, UpperLimit={6}", A.B, A.LM, A.LS2, A.M, A.N, A.LowerLimit, A.UpperLimit);
        }
    } //Sigma_LM_GenObject
}
