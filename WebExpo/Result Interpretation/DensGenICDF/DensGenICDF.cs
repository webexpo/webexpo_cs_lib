
namespace Zygotine.WebExpo
{
    using System;
    using System.Linq;
    using Zygotine.Util;

    internal class ICDFFunctions
    {
        internal static double PingpongTieBreaker(
            Func<double, SGNFnAParam, double> area,
            SGNFnAParam A,
            double start,
            double areaX,
            double target,
            double mathLowerLimit,
            double[] refRange,
            bool inestimableLowerLimit,
            double epsilon,
            bool distrnLeftSide = false,
            int maxCount = 100)
        {

            int hpMult = distrnLeftSide ? -1 : 1;
            bool cntn = true;
            int count = 0;
            double x = start;
            Visits visitedX = new Visits(maxCount);
            visitedX.Add(x);
            bool converged = false;
            bool caughtInLoop = false;
            while (cntn)
            {
                double hp = A.F(x, A) * hpMult;
                double change = (target - areaX) / hp;
                x = x + change;
                // new_0.11
                if (x < refRange[0])
                {
                    x = (refRange[0] == mathLowerLimit) && inestimableLowerLimit ? (x - change + mathLowerLimit) / 2.0 : refRange[0];
                }

                areaX = WebExpoFunctions3.SmoothedAreaEstimate(area, x, A, mathLowerLimit, refRange, inestimableLowerLimit, hpMult: hpMult);
                count++;
                converged = Math.Abs(areaX - target) < epsilon;
                visitedX.Add(x);
                caughtInLoop = visitedX.CaughtInLoop; // La propiété CaughtInLoop est mise à jour à chaque ajout d'un x.
                                                      // x doit être fini, et avoir déjà été visité.   
                cntn = !converged && (visitedX.Count <= maxCount) && !caughtInLoop;
            }

            if (!converged && caughtInLoop)
            {
                //  We have found the series of points that are repeatedly visited:
                // recalibrate and try Newton-Raphson again
                double[] tail = visitedX.GetFromTail(x);

                x = (distrnLeftSide) ? x = tail.Max() : tail.Min();
                areaX = WebExpoFunctions3.SmoothedAreaEstimate(area, x, A, mathLowerLimit, refRange, inestimableLowerLimit, hpMult: hpMult);

                if (distrnLeftSide)
                {
                    A.UpperLimit = x;
                }
                else
                {
                    A.LowerLimit = x;
                }

                target = target - areaX;
                areaX = 0.0;
                count = 0;
                cntn = true;
                while (cntn)
                {
                    double hp = A.F(x, A) * hpMult;
                    double change = (target - areaX) / hp;
                    x = x + change;
                    areaX = area(x, A);
                    count = count + 1;
                    converged = Math.Abs(areaX - target) < epsilon;
                    cntn = !converged && (count <= maxCount);
                }
            }

            if (!converged)
            {
                //TODO Exception
                throw new WEException("Newton-Raphson algorithm did not converge. Sorry.\n");
            }

            return x;
        } // end of ping.pong.tie.breaker
    }
}

