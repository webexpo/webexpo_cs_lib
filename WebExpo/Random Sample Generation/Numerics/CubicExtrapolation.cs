namespace Zygotine.Numerics
{
    using System.Collections.Generic;
    using System.Linq;
    using Zygotine.WebExpo;
    using Zygotine.Util;

    internal class CubicExtrapolation
    {

        private CubicExtrapolation(double target, Bounds bounds, bool monotonic, Range range)
        {
            //# target: scalar
            //# bounds: list

            HigherSide = bounds.HigherSide;
            Bounds = bounds;
            Monotonic = monotonic;
            Range = range;
        }

        internal static double Extrapolate(double target, Bounds bounds, bool monotonic, double[] range)
        {
            //# Approximate h by a cubic curve
            //# (with appropriate values and slopes at both ends)
            double[] tmp = new double[0];

            double[] bX = bounds.X;
            double[] bX2 = bX.Sqr();
            double[] bX3 = bX2.Multiply(bX);
            double[] X = new double[16]
            {
                1, 1, 0, 0,
                bX[0], bX[1], 1, 1,
                bX2[0], bX2[1], 2.0 * bX[0], 2.0 * bX[1],
                bX3[0], bX3[1], 3.0 * bX2[0], 3.0 * bX2[1]
            };

            double[] y = new double[4] { bounds.H[0], bounds.H[1], bounds.Hp[0], bounds.Hp[1] };
            double[] theta = XInverse4Y.Compute(X, y);

            bool useLinearExtrapolation = false;

            //TODO bien vérifier l'équivalence des conditions C# vs R
            //if (any(is.nan(theta)) | any(is.infinite(theta)))
            if (theta.Where(a => !a.IsFinite()).Count() > 0)
            {
                useLinearExtrapolation = true;
            }
            else
            {
                tmp = CubicEquation.GetRealRoots(theta, target, ROrder: true).ToArray();
                useLinearExtrapolation = tmp.Length == 0;// # modif_0.11
            }


            if (useLinearExtrapolation)
            {
                //# use a linear extrapolation instead since X is probably not invertible

                X = new double[4] { bX[0], bX[1], 1, 1 }; // for a 2 x 2 matrix
                theta = XInverse4Y.Compute(X, bounds.H);
                tmp = new double[] { (target - theta[1]) / theta[0] };
            }


            var inBoundSolns = tmp.Where(a => (a > bounds.X[0]) && (a <= bounds.X[1]));
            var inRangeSolns = tmp.Where(a => (a > range[0]) && (a <= range[1]));

            int m;
            double soln;
            if (bounds.IncludeSoln)
            {
                m = inBoundSolns.Count();
                if (m == 0)
                {
                    //# Use a linear extrapolation instead
                    X = new double[] { bX[0], bX[1], 1, 1 };  // 2 x 2 matrix
                    theta = XInverse4Y.Compute(X, bounds.H);
                    soln = (target - theta[1]) / theta[0];
                }
                else if (m == 1)
                {
                    soln = inBoundSolns.First();
                }
                else
                {
                    soln = inBoundSolns.ToArray().Mean();
                }
            }
            else
            {
                tmp = inRangeSolns.ToArray();
                m = tmp.Length;
                if (m == 0)
                {
                    soln = Tools.NA;
                }
                else if (m == 1)
                {
                    soln = tmp[0];
                }
                else if (monotonic)
                {

                    IEnumerable<double> tmpVar;
                    if (bounds.HigherSide == 1)
                    {
                        tmpVar = tmp.Where(a => a < bX[0]);

                        if (tmpVar.Count() > 0)
                        {
                            soln = tmpVar.Max();
                        }
                        else
                        {
                            soln = Tools.NA;
                        }
                    }
                    else
                    {
                        tmpVar = tmp.Where(a => a > bX[1]);
                        if (tmpVar.Count() > 0)
                        {
                            soln = tmp.Min();
                        }
                        else
                        {
                            soln = Tools.NA;
                        }
                    }
                }
                else
                {
                    //# keep the closest solution to any bound (or center of bounds, equivalently)
                    double[] d = tmp.Substract(bX.Mean()).Abs();
                    int? w = d.IndexOfMin();
                    //TODO verifier le code R, on assume qu'on a toujours une valeur.
                    soln = tmp[w.Value];
                }
            }

            return soln;
        } //# end of cubic.extrapolation

        internal int HigherSide { get; set; }
        internal bool Monotonic { get; set; }
        internal Bounds Bounds { get; set; }
        internal Range Range { get; set; }
    }
}
