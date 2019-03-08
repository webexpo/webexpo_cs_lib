
namespace Zygotine.Numerics
{
    using System;

    public static class Matrix2x2
    {
        private static bool Is2X2(double[][] m)
        {
            int nRow = m.Length;
            if ((m.Length == 2) && (m[0].Length == 2) && (m[1].Length == 2))
            {
                return true;
            }
            return false;
        }

        private static bool Is2X2(double[,] m)
        {
            int nRow = m.Length;
            if ((m.GetLength(0) == 2) && (m.GetLength(1) == 2))
            {
                return true;
            }
            return false;
        }

        internal static double[,] Solve2x2(double[,] m)
        {
            if (!Is2X2(m))
            {
                //TODO exception
                throw new Exception("Not a square matrix.");
            }

            double det = (m[0, 0] * m[1, 1]) - (m[0, 1] * m[1, 0]);
            if (det == 0)
            {
                return null;
            }

            double[,] mi = new double[2, 2];
            mi[0, 0] = m[1, 1] / det;
            mi[1, 1] = m[0, 0] / det;
            mi[1, 0] = -m[1, 0] / det;
            mi[0, 1] = -m[0, 1] / det;
            return mi;
        }

        internal static double[,] Solve2x2(double[] v, bool byCol = true)
        {
            double[,] vals = Create(v, byCol);
            return Solve2x2(vals);
        }

        internal static double[,] Create(double[] v, bool byCol = true)
        {
            double[,] rep = new double[2, 2];
            rep[0, 0] = v[0];
            rep[1, 1] = v[3];

            if (byCol)
            {
                rep[1, 0] = v[1];
                rep[0, 1] = v[2];
            }
            else
            {
                rep[0, 1] = v[1];
                rep[1, 0] = v[2];
            }

            return rep;
        }

        internal static double[] Product(double[,] m, double[] v)
        {
            return new double[2] { m[0, 0] * v[0] + m[0, 1] * v[1], m[1, 0] * v[0] + m[1, 1] * v[1] };
        }

    }
}

