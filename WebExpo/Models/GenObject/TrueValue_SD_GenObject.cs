namespace Zygotine.WebExpo
{
    using System;
    using System.Linq;
    using System.Collections.Generic;
    using Zygotine.Numerics;
    using Zygotine.Statistics.Distribution;
    using Zygotine.Util;

    internal class TrueValue_SD_GenObject : GenObject
    {
        internal TrueValue_SD_GenObject() : base()
        {
            this.A = new SGNFnAParam();
            this.A.M = 0;
            this.Range = Tools.Combine(double.NegativeInfinity, double.PositiveInfinity);
        }
        /*
            f <- function(s, A){exp(- pnorm(exp(s)/A$ksi, log.p=T) - (A$y-exp(s))^2/2/A$ksi2 - (s-A$mu)^2/2/A$sigma2 - A$M)}
            log.f <- function(s, A){- pnorm(exp(s)/A$ksi, log.p=T) - (A$y-exp(s))^2/2/A$ksi2 - (s-A$mu)^2/2/A$sigma2} 
            log.f.prime <- function(s, A){z <- exp(s)/A$ksi; l <- lphi(z); -l$r*z + A$y*z/A$ksi - z^2 - (s-A$mu)/A$sigma2}
            log.f.second <- function(s, A){z <- exp(s)/A$ksi; l <- lphi(z); l$r*(z^3-z) + (l$r*z)^2 + (A$y*exp(s)-2*exp(2*s))/A$ksi2 - 1/A$sigma2}

        */
        internal override double F(double s, SGNFnAParam A)
        {
            return Math.Exp(
                -NormalDistribution.PNorm(Math.Exp(s) / A.Ksi, log_p: true)
                - Math.Pow((A.Y - Math.Exp(s)), 2) / 2.0 / A.Ksi2
                - Math.Pow((s - A.Mu), 2) / 2.0 / A.Sigma2
                - A.M);
        }

        internal override double LogF(double s, SGNFnAParam A)
        {
            return
                -NormalDistribution.PNorm(Math.Exp(s) / A.Ksi, log_p: true)
                - Math.Pow((A.Y - Math.Exp(s)), 2) / 2.0 / A.Ksi2
                - Math.Pow((s - A.Mu), 2) / 2.0 / A.Sigma2;
        }

        internal override double LogFPrime(double s, SGNFnAParam A)
        {
            //    log.f.prime <- function(s, A){z <- exp(s)/A$ksi; l <- lphi(z); -l$r*z + A$y*z/A$ksi - z^2 - (s-A$mu)/A$sigma2}
            double z = Math.Exp(s) / A.Ksi;
            LPhiObject lphio = WebExpoFunctions3.LPhi(z);
            return 
                    -lphio.R[0] * z 
                    + A.Y * z / A.Ksi 
                    - Math.Pow(z, 2) - (s - A.Mu) / A.Sigma2;
        }

        internal override double LogFSecond(double s, SGNFnAParam A)
        {
            //    log.f.prime <- function(s, A){z <- exp(s)/A$ksi; l <- lphi(z); -l$r*z + A$y*z/A$ksi - z^2 - (s-A$mu)/A$sigma2}
            double[] z = new double[] { Math.Exp(s) / A.Ksi };
            LPhiObject lphio = WebExpoFunctions3.LPhi(z);
            double r0 = lphio.R[0];
            double z0 = z[0];
            return r0 * (Math.Pow(z0, 3) - z[0]) + Math.Pow(r0 * z0, 2) + (A.Y * Math.Exp(s) - 2.0 * Math.Exp(2 * s)) / A.Ksi2 - 1.0 / A.Sigma2;
        }

        internal override List<double> LogFInvRemote(double target, SGNFnAParam A)
        {
            /*
                   # solutions for large s
                  tmp <- A$y + c(-1,1)*sqrt(2*abs(target)*A$ksi2) # vector of length 2
                  tmp <- tmp[tmp>0]
                  roots <- log(tmp)

             */
            // solutions for large s
            List<double> roots = new List<double>();
            double z = Math.Sqrt(2 * Math.Abs(target) * A.Ksi2);
            double pp = A.Y - z;
            double pg = A.Y + z;

            if (pp > 0.0)
            {
                roots.Add(Math.Log(pp));
                roots.Add(Math.Log(pg));
            }
            else if (pg > 0.0)
            {
                roots.Add(Math.Log(pg));
            }

            /*
                  C <- log(2) -A$y2/2/A$ksi2 - A$mu^2/2/A$sigma2
                  B <- A$mu/A$sigma2
                  qA <- -1/2/A$sigma2
                  tmp <- quadratic.solution(c(C, B, qA), target=target)
                  roots <- c(roots, tmp)

             */

            //# solutions for small s (large negative values)
            double C = Math.Log(2) - A.Y2 / 2.0 / A.Ksi2 - Math.Pow(A.Mu, 2) / 2.0 / A.Sigma2;
            double B = A.Mu / A.Sigma2;
            double qA = -1.0 / 2.0 / A.Sigma2;
            roots.AddRange(QuadraticEquation.GetRealRoots(Tools.Combine(C, B, qA), target, ROrder: true));

            /*
                  C <- log(2) -A$y2/2/A$ksi2 - A$mu^2/2/A$sigma2
                  B <- A$mu/A$sigma2
                  qA <- -1/2/A$sigma2
                  tmp <- quadratic.solution(c(C, B, qA), target=target)
                  roots <- c(roots, tmp)
             */

            //# solutions around 0
            C = Math.Log(1.0 / A.Ksi) - Math.Pow(A.Y - 1.0, 2) / 2.0 / A.Ksi2 - Math.Pow(A.Mu, 2) / 2.0 / A.Sigma2;
            B = A.Mu / A.Sigma2;
            qA = -1.0 / 2.0 / A.Sigma2;
            roots.AddRange(QuadraticEquation.GetRealRoots(Tools.Combine(C, B, qA), target, ROrder: true));
            return roots;
        } // # end of log.f.inv.remote

        internal override double[] Start(SGNFnAParam A)
        {
            return Tools.Combine(A.LogY, A.Mu);
        }
    }// TrueValue4SD
}
