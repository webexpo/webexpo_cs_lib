namespace Zygotine.WebExpo
{
    using System;
    using System.Linq;
    using System.Collections.Generic;
    using Zygotine.Numerics;
    using Zygotine.Statistics.Distribution;
    using Zygotine.Util;

    internal class TrueValue_CV_Norm_GenObject : GenObject
    {
        public TrueValue_CV_Norm_GenObject()
        {
            this.A = new SGNFnAParam();
            this.A.M = 0;
            this.Range = Tools.Combine(0, double.PositiveInfinity);
            this.PotentiallyBimodal = true;
            this.MathLowerLimit = 0;
        }

        /*
    f <- function(t, A){exp(-log(t) - (A$y-t)^2/2/A$cv2/t^2 - (t-A$mu)^2/2/A$sigma2 - A$M)}
    log.f <- function(t, A){-log(t) - (A$y-t)^2/2/A$cv2/t^2 - (t-A$mu)^2/2/A$sigma2}       
    log.f.prime <- function(t, A){-1/t - (A$y/t^2-A$y2/t^3)/A$cv2 + (A$mu-t)/A$sigma2}   
    log.f.second <- function(t, A){1/t^2 + (2*A$y/t^3 - 3*A$y2/t^4)/A$cv2 -1/A$sigma2}
         */

        internal override double F(double t, SGNFnAParam A)
        {
            return Math.Exp(
                            -Math.Log(t)
                            - Math.Pow(A.Y - t, 2) / 2.0 / A.CV2 / Math.Pow(t, 2)
                            - Math.Pow(t - A.Mu, 2) / 2.0 / A.Sigma2 - A.M);
        }

        internal override double LogF(double t, SGNFnAParam A)
        {
            return
                    -Math.Log(t)
                    - Math.Pow(A.Y - t, 2) / 2.0 / A.CV2 / Math.Pow(t, 2)
                    - Math.Pow(t - A.Mu, 2) / 2.0 / A.Sigma2;
        }

        internal override double LogFPrime(double t, SGNFnAParam A)
        {
            return
                    -1.0 / t
                    - (A.Y / Math.Pow(t, 2) - A.Y2 / Math.Pow(t, 3)) / A.CV2
                    + (A.Mu - t) / A.Sigma2;
        }

        internal override double LogFSecond(double t, SGNFnAParam A)
        {
            return
                1.0 / Math.Pow(t, 2)
                + (2.0 * A.Y / Math.Pow(t, 3)
                - 3.0 * A.Y2 / Math.Pow(t, 4)) / A.CV2
                - 1.0 / A.Sigma2;
        }


        internal override List<double> LogFInvRemote(double target, SGNFnAParam A)
        {
            List<double> roots = new List<double>();
            roots.Add(A.Y / (1 + Math.Sqrt(-2.0 * A.CV2 * target)));

            //# We try a solution with a small t
            double D = -A.Y2 / 2.0 / A.CV2;
            double C = A.Y / A.CV2;
            /*
             * 3/2 - 1/2/A$cv2 - A$mu^2/2/A$sigma2 - target
            */
            double B =
                        3.0 / 2.0
                        - 1.0 / 2.0 / A.CV2
                        - Math.Pow(A.Mu, 2) / 2.0 / A.Sigma2
                        - target;
            double qA = -2.0 + A.Mu / A.Sigma2;
            double[] theta = Tools.Combine(D, C, B, qA);
            roots.AddRange(CubicEquation.GetRealRoots(theta, ROrder: true));
            // Look also for a solution with a large value for t
            double tmp = -2.0 * A.Sigma2 * (target + 1 / 2.0 / A.CV2);
            if (tmp > 0)
            {
                double racine = A.Mu + Math.Sqrt(tmp);
                if (racine > 0.0)
                {
                    roots.Add(racine);
                }
            }

            return roots;
        }

        internal override double[] Start(SGNFnAParam A)
        {
            List<double> roots = new List<double>();
            //# solutions in small values for t
            double C = A.Y2 / A.CV2;
            double B = -A.Y / A.CV2;
            roots.AddRange(QuadraticEquation.GetRealRoots(Tools.Combine(C, B, -1), l: 0, ROrder: true));
            //# solutions in large values for t
            double D = A.Y / A.CV2;
            B = -A.Mu / A.Sigma2;
            double qA = 1 / A.Sigma2;
            double[] theta = Tools.Combine(D, 1, B, qA);
            roots.AddRange(CubicEquation.GetRealRoots(theta, l: 0, ROrder: true));
            return roots.ToArray();
        } //# end of start
    }
}
