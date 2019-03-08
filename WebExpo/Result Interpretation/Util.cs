using System;
using System.Collections.Generic;

namespace Zygotine.Util
{
    public static partial class Tools
    {


        public static double[] Combine(params double[] values)
        {
            return values;
        }

        public static double[] Combine(params double[][] values)
        {

            if (values == null)
            {
                return null;
            }

            if (values.Length == 0)
            {
                return new double[0];
            }

            int len = 0;
            foreach (double[] t in values)
            {
                len += t.Length;
            }

            double[] rep = new double[len];
            int p = 0;
            for (int i = 0; i < values.Length; i++)
            {
                if (values[i]?.Length != 0)
                {
                    Array.Copy(values[i], 0, rep, p, values[i].Length);
                    p += values[i].Length;
                }
            }

            return rep;
        }


        public static int[] Combine(params int[][] values)
        {

            if (values == null)
            {
                return null;
            }

            if (values.Length == 0)
            {
                return new int[0];
            }

            int len = 0;
            foreach (int[] t in values)
            {
                len += t.Length;
            }

            int[] rep = new int[len];
            int p = 0;
            for (int i = 0; i < values.Length; i++)
            {
                if (values[i]?.Length != 0)
                {
                    Array.Copy(values[i], 0, rep, p, values[i].Length);
                    p += values[i].Length;
                }
            }

            return rep;
        }



        public static string DebugShow(params object[] values)
        {

            string s = "";
            foreach (object value in values)
            {
                s += value.ToString() + " - ";
            }
            return s;
        }

        public static double[] GetConstantArray(int n, double value)
        {
            double[] rep = new double[n];
            if (value != 0.0)
                for (int i = 0; i < n; i++)
                    rep[i] = value;
            return rep;
        }

        public static int[] IntegralSequence(int from, int to, int by)
        {
            if (
                ((by > 0) && (to <= from)) ||
                ((by < 0) && (to >= from)))
                return new int[1] { from };

            int[] y = new int[((to - from) / by) + 1];
            for (int i = 0; i < y.Length; i++)
            {
                y[i] = from + (i * by);
            }
            return y;
        }

        public static double RoundToSignificantDigits(this double d, int digits)
        {
            if (d == 0)
            {
                return 0;
            }

            double scale = Math.Pow(10, Math.Floor(Math.Log10(Math.Abs(d))) + 1);
            return scale * Math.Round(d / scale, digits);
        }

        public static double[] Extract(this double[] src, int[] index)
        {
            if(src == null)
            {
                return null;
            }

            if(src.Length == 0)
            {
                return new double[0];
            }

            List<double> rep = new List<double>();

            foreach (int i in index)
            {
                if(i>src.Length)
                {
                    rep.Add(Tools.NA);
                }
                else
                {
                    rep.Add(src[i]);
                }
            }

            return rep.ToArray();
        }

        public static double[] Sequence(double from, double to, double by)
        {

            double del = to - from, dd;
            double n = del / by;

            int nn = (int)(n + 1e-10);
            double[] rep = new double[nn + 1];
            dd = Math.Abs(del) / Math.Max(Math.Abs(to), Math.Abs(from));
            for (int i = 0; i <= nn; i++)
                rep[i] = from + (double)i * by;

            if (nn > 0)
                if ((by > 0 && rep[nn] > to) || (by < 0 && rep[nn] < to))
                    rep[nn] = to;
            return rep;
        }

  
    }
}
