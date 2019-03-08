namespace Zygotine.WebExpo
{
    using System;
    using System.Linq;
    using System.Collections.Generic;
    using Zygotine.Numerics;
    using Zygotine.Statistics.Distribution;
    using Zygotine.Util;

    internal class TrueValue_CV_LogN_GenObject : GenObject
    {
        // ME ~ CV && distrn ~ LogNormal
        internal TrueValue_CV_LogN_GenObject()
        {
            this.A = new SGNFnAParam();
            this.A.M = 0;
            this.Range = Tools.Combine(double.NegativeInfinity, double.PositiveInfinity);
            this.PotentiallyBimodal = true;
            this.MathLowerLimit = double.NegativeInfinity;
            this.LogNormalDistrn = true;
            this.ThroughCV = true;
            this.ThroughSD = false;
        }
 
        /*
            f <- function(s, A){exp(-s - (A$y*exp(-s)-1)^2/2/A$cv2 - (s-A$mu)^2/2/A$sigma2 - A$M)}
            log.f <- function(s, A){-s - (A$y*exp(-s)-1)^2/2/A$cv2 - (s-A$mu)^2/2/A$sigma2}
            log.f.prime <- function(s, A){-1 + (A$y2*exp(-2*s) - A$y*exp(-s))/A$cv2 - (s-A$mu)/A$sigma2}
            log.f.second <- function(s, A){-2*A$y2*exp(-2*s)/A$cv2 + A$y*exp(-s)/A$cv2 - 1/A$sigma2}
        */
        internal override double F(double s, SGNFnAParam A)
        {
            return Math.Exp(
                -s
                - Math.Pow(A.Y * Math.Exp(-s) - 1, 2) / 2.0 / A.CV2
                - Math.Pow(s - A.Mu, 2) / 2.0 / A.Sigma2
                - A.M);
        }

        internal override double LogF(double s, SGNFnAParam A)
        {
            return
                -s
                - Math.Pow(A.Y * Math.Exp(-s) - 1.0, 2) / 2.0 / A.CV2
                - Math.Pow(s - A.Mu, 2) / 2.0 / A.Sigma2;
        }

        internal override double LogFPrime(double s, SGNFnAParam A)
        {
            return -1.0
                    + (
                          A.Y2 * Math.Exp(-2 * s)
                        - A.Y * Math.Exp(-s)
                      ) / A.CV2
                   - (s - A.Mu) / A.Sigma2;
        }

        internal override double LogFSecond(double s, SGNFnAParam A)
        {
            return
                -2.0 * A.Y2 * Math.Exp(-2 * s) / A.CV2
                + A.Y * Math.Exp(-s) / A.CV2
                - 1 / A.Sigma2;
        }
        /*
        log.f.inv.remote <- function(target, A)
        {
          tmp <- 1 + c(-1,1)*sqrt(2*A$cv2*abs(target))
          tmp <- tmp[tmp>0]
          roots <- log(tmp/A$y)
    
          C <- -(A$y-1)^2/2/A$cv2 - A$mu^2/2/A$sigma2
          B <- -1 + (A$y-1)*A$y/A$cv2 + A$mu/A$sigma2
          qA <- -A$y2/2/A$cv2 - 1/2/A$sigma2
          tmp <- quadratic.solution(c(C, B, qA), target=target)
      
          roots <- c(roots, tmp)      
          roots
        } # end of log.f.inv.remote
      }        
        */
        internal override List<double> LogFInvRemote(double target, SGNFnAParam A)
        {
            List<double> roots = new List<double>();
            double z = Math.Sqrt(2.0 * A.CV2 * Math.Abs(target));
            double pPetit = 1.0 - z;
            double pGrand = 1.0 + z;
            if (pPetit > 0.0)
            {
                roots.Add(Math.Log(pPetit / A.Y));
                roots.Add(Math.Log(pGrand / A.Y));
            }
            else if (pGrand > 0.0)
            {
                roots.Add(Math.Log(pGrand / A.Y));
            }

            double C =
                        -Math.Pow(A.Y - 1.0, 2) / 2.0 / A.CV2
                        - Math.Pow(A.Mu, 2) / 2.0 / A.Sigma2;
            double B =
                        -1.0
                        + (A.Y - 1.0) * A.Y / A.CV2
                        + A.Mu / A.Sigma2;
            double qA =
                        -A.Y2 / 2.0 / A.CV2
                        - 1.0 / 2.0 / A.Sigma2;

            roots.AddRange(QuadraticEquation.GetRealRoots(Tools.Combine(C, B, qA), target, ROrder: true));
            return roots;
        }

        /*
            start <- function(A)
            {
              C <- -1 + A$y2/A$cv2 - A$y/A$cv2 + A$mu/A$sigma2
              B <- -2*A$y2/A$cv2 + A$y/A$cv2 - 1/A$sigma2
              qA <- 2*A$y2/A$cv2 - A$y/2/A$cv2
              tmp <- quadratic.solution(c(C, B, qA))
      
              start <- c(tmp, A$mu - A$sigma2)   
              start
            } # end of start
        */

        internal override double[] Start(SGNFnAParam A)
        {
            double C = -1.0 + A.Y2 / A.CV2 - A.Y / A.CV2 + A.Mu / A.Sigma2;
            double B = -2.0 * A.Y2 / A.CV2 + A.Y / A.CV2 - 1.0 / A.Sigma2;
            double qA = 2.0 * A.Y2 / A.CV2 - A.Y / 2.0 / A.CV2;
            List<double> roots = QuadraticEquation.GetRealRoots(Tools.Combine(C, B, qA), ROrder: true);
            roots.Add(A.Mu - A.Sigma2);
            return roots.ToArray();
        }
    } // 
}
