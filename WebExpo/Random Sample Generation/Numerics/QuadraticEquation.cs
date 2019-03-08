using System;
using System.Numerics;
using System.Collections.Generic;

namespace Zygotine.Numerics
{
    public static class QuadraticEquation
    {
        public static List<Complex> GetRoots(double[] th, double target = 0.0, bool ROrder = false)
        {
            double A, C;
            double B = th[1];
            if (ROrder)
            {
                A = th[2];
                C = th[0] - target;
            }
            else
            {
                A = th[0];
                C = th[2] - target;
            }

            List<Complex> soln = new List<Complex>();
            if (A == 0)
            {
                if (B == 0)
                {
                    return soln; // vide
                }
                soln.Add(-C / B);
                return soln;
            }

            double delta = (B * B) - 4 * A * C;
            double AA = 2 * A;
            if (delta < 0)
            {
                Complex sqrtDelta = Complex.Sqrt(delta);
                Complex racine = racine = (-B + sqrtDelta) / AA;
                soln.Add(racine);
                racine = (-B - sqrtDelta) / AA;
                soln.Add(racine);
            }
            else if (delta > 0)
            {
                double sqrtDelta = Math.Sqrt(delta);
                soln.Add((-B + sqrtDelta) / AA);
                soln.Add((-B - sqrtDelta) / AA);
            }
            else //delta == 0
            {
                soln.Add(-B);
            }

            return soln;
        }

        internal static List<double> GetRealRoots(double[] theta, double target = 0.0, double l = double.NegativeInfinity, double u = double.PositiveInfinity, bool ROrder = false)
        {
            List<double> soln = new List<double>();
            List<Complex> cplxSoln = GetRoots(theta, target, ROrder);
            
            foreach(Complex c in cplxSoln)
            {
                if((c.Imaginary == 0) && (c.Real > l) && (c.Real < u))
                {
                    soln.Add(c.Real);
                }
            }
            return soln;
        }
    }
}
