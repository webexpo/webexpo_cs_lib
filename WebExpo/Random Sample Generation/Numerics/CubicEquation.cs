namespace Zygotine.Numerics
{
    using System;
    using System.Collections.Generic;
    using System.Linq;
    using System.Numerics;
    using Zygotine.Util;

    public static class CubicEquation
    {
        public static List<double> GetRealRoots(double[] theta, double target = 0, double tolerance = 1E-6, double l = double.NegativeInfinity, double u = double.PositiveInfinity, bool ROrder = false, bool useRPoly = false)
        {
            // La fonction de R polyroot est une implémentation de CPOLY (https://en.wikipedia.org/wiki/Jenkins%E2%80%93Traub_algorithm)
            // Ici nous utilisons une traduction en C# de RPOLY.

            List<Complex> cRoots = GetRoots(theta: theta, target: target, ROrder: ROrder, useRPoly:useRPoly);
            IEnumerable<Complex> newRoots;
            if (tolerance > 0)
            {
                newRoots = cRoots.Where(c => Math.Abs(c.Imaginary) < tolerance);
            }
            else
            {
                newRoots = cRoots;
            }
            
            IEnumerable<double> realRoots = newRoots.Select(c => c.Real).Where(d => (d >= l) && (d <= u));
            return realRoots.ToList();
        } //# end of real.cubic.roots

        public static List<Complex> GetRoots(double[] theta, double target = 0, bool ROrder = false, bool useRPoly = true)
        {
            double a = ROrder ? theta[3] - target : theta[0] - target;
            if (a == 0)
            {
                if (ROrder)
                {
                    return QuadraticEquation.GetRoots(new double[] { theta[0], theta[1], theta[2] }, target, ROrder);
                }
                else
                {
                    return QuadraticEquation.GetRoots(new double[] { theta[1], theta[2], theta[3] }, target, ROrder);
                }
            }
            else
            {
                List<Complex> lc = new List<Complex>();
                double[] th = theta.Copy();
                Func<double[], List<Complex>> rootFinder;
                if( useRPoly)
                {
                    rootFinder = RealPolynomialRootFinder.FindRoots;
                }
                else
                {
                    rootFinder = CPoly.GetRoots;
                }

                if (ROrder)
                {
                    Array.Reverse(th);
                    th[3] = th[3] - target;
                    lc.AddRange(rootFinder(th));
                }
                else
                {
                    th[0] = th[0] - target;
                    lc.AddRange(rootFinder(th));
                }
                return lc;
            }
        }
    }
}