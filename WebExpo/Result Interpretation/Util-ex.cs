using System;
using System.Collections.Generic;
using System.Linq;
using System.Numerics;
using System.Text;

namespace Zygotine.Util
{

    public static partial class Tools
    {


        public static double[] Abs(this double[] x)
        {
            return Apply(x, Math.Abs);
        }

        public static double[] Add(this double[] x, double[] y)
        {
            return DoBinaryOperation(x, y, Addition);
        }

        public static double[] Add(this double[] x, double y)
        {
            return DoBinaryOperation(x, new double[1] { y }, Addition);
        }

        public static double[] Add(this double x, double[] y)
        {
            return DoBinaryOperation(new double[1] { x }, y, Addition);
        }

        //public static double[] AdditiveInverse(this double[] x)
        //{
        //    if (x == null)
        //    {
        //        return null;
        //    }
        //    double[] y = new double[x.Length];
        //    for (int i = 0; i < x.Length; i++)
        //        y[i] = -x[i];
        //    return y;
        //}

        //public static bool Any(this double[] x, Func<double, bool> test)
        //{
        //    if ((x == null) || (x.Length == 0))
        //        return false;
        //    for (int i = 0; i < x.Length; i++)
        //    {
        //        if (test(x[i]))
        //            return true;
        //    }
        //    return false;
        //}

        internal static double[] Apply(this double[] x, Func<double, double> f)
        {
            if (x == null)
                return null;
            if (x.Length == 0)
                return new double[0];

            double[] rep = new double[x.Length];
            for (int i = 0; i < x.Length; i++)
                rep[i] = f(x[i]);
            return rep;
        }

        public static double[] Combine(this double x, double y)
        {
            return new double[2] { x, y };
        }

        public static double[] Combine(this double[] x, double y)
        {
            if ((x == null) || (x.Length == 0))
                return new double[1] { y };

            double[] rep = new double[x.Length + 1];
            Array.Copy(x, rep, x.Length);
            rep[rep.Length - 1] = y;
            return rep;
        }

        public static double[] Combine(this double x, double[] y)
        {
            if ((y == null) || (y.Length == 0))
                return new double[1] { x };
            double[] rep = new double[y.Length + 1];
            Array.Copy(y, 0, rep, 1, y.Length);
            rep[0] = x;
            return rep;
        }

        public static double[] Combine(this double[] x, double[] y)
        {
            double[] rep;
            if ((x == null) || (x.Length == 0))
            {
                if ((y == null) || (y.Length == 0))
                {
                    return new double[0];
                }

                rep = new double[y.Length];
                Array.Copy(y, rep, y.Length);
                return rep;
            }
            else
            {
                if ((y == null) || (y.Length == 0))
                {
                    rep = new double[x.Length];
                    Array.Copy(x, rep, x.Length);
                    return rep;
                }
                rep = new double[x.Length + y.Length];
                Array.Copy(x, rep, x.Length);
                Array.Copy(y, 0, rep, x.Length, y.Length);
                return rep;
            }
        }

        public static string[] Copy(this string[] x)
        {
            if (x == null)
            {
                return null;
            }

            string[] rep = new string[x.Length];
            Array.Copy(x, 0, rep, 0, x.Length);
            return rep;
        }

        public static double[] Copy(this double[] x)
        {
            if (x == null)
            {
                return null;
            }
            double[] rep = new double[x.Length];
            Array.Copy(x, 0, rep, 0, x.Length);
            return rep;
        }

        public static double[] Divide(this double[] x, double[] y)
        {
            return DoNonCommutativeBinaryOperation(x, y, Division);
        }

        public static double[] Divide(this double[] x, int[] y)
        {
            if (y == null)
                return new double[0];
            double[] y1 = new double[y.Length];
            for (int i = 0; i < y.Length; i++)
                y1[i] = (double)y[i];
            return DoBinaryOperation(x, y1, Division);
        }


        public static double[] Divide(this double[] x, double divisor)
        {
            return DoNonCommutativeBinaryOperation(x, new double[1] { divisor }, Division);
        }

        public static double[] Divide(this double x, double[] y)
        {
            return DoNonCommutativeBinaryOperation(new double[] { x }, y, Division);
        }

        public static double[] Divide(this double x, int[] y)
        {
            double[] yy = y.Select(Y => (double)Y).ToArray();
            return DoNonCommutativeBinaryOperation(new double[] { x }, yy, Division);
        }

        public static double[] Exp(this double[] x)
        {
            return Apply(x, Math.Exp);
        }

        public static string FormatDouble(this double x)
        {
            return IsNA(x) ? "NA" : x.ToString();
        }

        public static double[] FromIndices(this double[] x, int[] indices)
        {
            if ((x == null) || (x.Length == 0))
            {
                return new double[0];
            }

            List<double> ld = new List<double>();
            foreach (int i in indices)
            {
                if (i >= 0 & i < x.Length)
                    ld.Add(x[i]);
            }
            return ld.ToArray();
        }

        internal static int? IndexOfMin(this double[] x)
        {
            switch (x.Length)
            {
                case 0:
                    return null;
                case 1:
                    return 0;
                default:
                    int pos = 0;
                    for (int i = 1; i < x.Length; i++)
                    {
                        if (x[i] < x[pos])
                        {
                            pos = i;
                        }
                    }
                    return pos;
            }
        }

        internal static int? IndexOfMax(this double[] x)
        {
            switch (x.Length)
            {
                case 0:
                    return null;
                case 1:
                    return 0;
                default:
                    int pos = 0;
                    for (int i = 1; i < x.Length; i++)
                    {
                        if (x[i] > x[pos])
                        {
                            pos = i;
                        }
                    }
                    return pos;
            }
        }

        /// <summary>
        /// Checks if the given <paramref name="x"/> has finite value.
        /// </summary>
        /// <param name="number">A floating-point number.</param>
        /// <returns><c>true</c> if <paramref name="x"/> is finite, <c>false</c> otherwise.</returns>
        /// <remarks>
        /// <para>
        /// A double is finite if it is neither double.NaN, double.NegativeInfinity nor double.PositiveInfinity
        /// </para>
        /// </remarks>
        internal static bool IsFinite(this double x)
        {
            // Si l'exposant est composé de 1, on a Double.NegativeInfinity || Double.PositiveInfinity || Double.NaN
            return ((BitConverter.DoubleToInt64Bits(x) & _exponentMask) != _exponentMask);
        }

        public static bool IsNA(this double d)
        {
            // la comparaison n'est pas parfaite.
            // Il faudrait exiger que les octets de 1 à 5 soient à zéro et le bit du signe à 1.
            // On ne teste donc pas l'égalité de d et de NA.
            return !double.IsNaN(d) ? false : BitConverter.GetBytes(d)[6] == 252;

        }

        public static bool IsNaN(this double d)
        {
            return double.IsNaN(d);
        }

        public static bool IsND(this double d)
        {
            // la comparaison n'est pas parfaite.
            // Il faudrait exiger que les octets de 1 à 5 soient à zéro et le bit du signe à 1.
            // On ne teste donc pas l'égalité de d et de NA.
            return !double.IsNaN(d) ? false : BitConverter.GetBytes(d)[6] == 254;
        }

        public static double[] Log(this double[] x)
        {
            return Apply(x, Math.Log);
        }

        internal static long Mantissa(this double x)
        {
            return BitConverter.DoubleToInt64Bits(x) & _mantissaMask;
        }

        public static double Max(this double[] x)
        {
            double maximum = double.NegativeInfinity;
            if (x == null)
            {
                return maximum;
            }

            for (int i = 0; i < x.Length; i++)
                if (x[i] > maximum)
                {
                    maximum = x[i];
                }

            return maximum;
        }


        public static double Mean(this IEnumerable<double> x, bool skipNaN = false)
        {
            //en R, mean(NULL) == NA
            if (x == null)
            {
                return Tools.NA;
            }

            double cumul = 0;
            int n = 0;
            foreach (double d in x)
            {
                if (double.IsNaN(d) && skipNaN)
                {
                    continue;
                }
                cumul += d;
                n++;
            }

            //en R, mean(numeric(0) == NaN)
            return n == 0 ? double.NaN : cumul / n;
        }

        public static double Min(this double[] x)
        {
            double minimum = double.PositiveInfinity;
            if (x == null)
            {
                return minimum;
            }

            for (int i = 0; i < x.Length; i++)
                if (x[i] < minimum)
                {
                    minimum = x[i];
                }
            return minimum;
        }

        public static double[] MinMax(this double[] x)
        {
            double[] y = new double[2] { double.PositiveInfinity, double.NegativeInfinity };
            for (int i = 0; i < x.Length; i++)
            {
                if (x[i] < y[0])
                {
                    y[0] = x[i];
                }
                if (x[i] > y[1])
                {
                    y[1] = x[i];
                }
            }
            return y;
        }

        public static double Prod(this double[] x)
        {
            if (x == null) 
            {
                return 1.0;
            }

            double rep = 1.0;
            foreach(double d in x)
            {
                rep *= d;
            }

            return rep;
        }

        public static double[] Reciprocal(this double[] x)
        {
            if (x == null)
            {
                return null;
            }
            double[] y = new double[x.Length];
            for (int i = 0; i < x.Length; i++)
                y[i] = 1 / x[i];
            return y;
        }

        public static double[] Multiply(this double[] x, double[] y)
        {
            return DoBinaryOperation(x, y, Multiplication);
        }

        public static double[] Multiply(this double[] x, int[] y)
        {
            if (y == null)
                return new double[0];
            double[] y1 = new double[y.Length];
            for (int i = 0; i < y.Length; i++)
            {
                y1[i] = (double)y[i];
            }
            return DoBinaryOperation(x, y1, Multiplication);

        }

        public static double[] Multiply(this double[] x, double y)
        {
            return DoBinaryOperation(x, new double[1] { y }, Multiplication);
        }

        public static double[] Multiply(this double[] x, int y)
        {
            return DoBinaryOperation(x, new double[1] { y }, Multiplication);
        }

        public static double[] Multiply(this double x, double[] y)
        {
            return DoBinaryOperation(new double[1] { x }, y, Multiplication);
        }

        public static double[] ParallelMax(this double[] x, double scalar)
        {
            if (x == null)
                return null;
            if (x.Length == 0)
                return new double[0];
            return x.Select(a => (a > scalar ? a : scalar)).ToArray();
        }

        public static double[] ParallelMax(this double[] x, double[] y)
        {
            double[] AShorter, ALonger;
            if (x.Length <= y.Length)
            {
                AShorter = x;
                ALonger = y;
            }
            else
            {
                AShorter = y;
                ALonger = x;
            }
            double[] rep = new double[ALonger.Length];
            for (int iL = 0, iS = 0; iL < ALonger.Length; iL++, iS = (iS + 1) % AShorter.Length)
            {
                double l = ALonger[iL];
                double s = AShorter[iS];
                if (double.IsNaN(l) || double.IsNaN(s))
                {
                    rep[iL] = l + s; // NA ou NaN
                }
                else
                {
                    rep[iL] = l > s ? l : s;
                }
            }
            return rep;
        }

        public static double[] ParallelMin(this double[] x, double scalar)
        {
            if (x == null)
                return null;
            if (x.Length == 0)
                return new double[0];
            return x.Select(a => (a < scalar ? a : scalar)).ToArray();
        }

        public static double[] ParallelMin(this double[] x, double[] y)
        {
            double[] AShorter, ALonger;
            if (x.Length <= y.Length)
            {
                AShorter = x;
                ALonger = y;
            }
            else
            {
                AShorter = y;
                ALonger = x;
            }

            double[] rep = new double[ALonger.Length];
            for (int iL = 0, iS = 0; iL < ALonger.Length; iL++, iS = (iS + 1) % AShorter.Length)
            {
                double l = ALonger[iL];
                double s = AShorter[iS];
                if (double.IsNaN(l) || double.IsNaN(s))
                {
                    rep[iL] = l + s; // NA ou NaN
                }
                else
                {
                    rep[iL] = l < s ? l : s;
                }

            }
            return rep;
        }

        public static double[] Remove(this double[] x, int position)
        {
            // Handles negative ends.

            if (position < 0)
            {
                position = x.Length + position;
            }
            if ((position < 0) || (position >= x.Length))
            {
                throw new ArgumentOutOfRangeException(nameof(Remove));
            }
            // Return new array.
            double[] res = new double[x.Length - 1];
            int j = 0;
            for (int i = 0; i < x.Length; i++)
            {
                if (i == position)
                    continue;

                res[j] = x[i];
                j++;
            }
            return res;
        }

        public static double[] Remove(this double[] x, int[] positions)
        {
            if (x == null)
                return null;
            if (x.Length == 0)
            {
                return new double[0];
            }

            if (positions.Any(i => i < 0))
            {
                throw new ArgumentOutOfRangeException(nameof(SubArray));
            }

            List<double> l = new List<double>();
            int previousPosition = -1;
            foreach (int p in positions.OrderByDescending(s => s))
            {
                if (p != previousPosition)
                {
                    l.Add(x[p]);
                    previousPosition = p;
                }
            }
            l.Reverse();
            return l.ToArray();
        }

        public static double[] Rep(this double value, int times)
        {
            if (times < 0)
            {
                return null;
            }

            if (times == 0)
            {
                return new double[0];
            }

            double[] tmp = new double[times];
            for (int i = 0; i < times; i++)
            {
                tmp[i] = value;
            }

            return tmp;
        }

        public static bool[] Rep(this bool value, int times)
        {
            if (times < 0)
            {
                return null;
            }

            if (times == 0)
            {
                return new bool[0];
            }

            bool[] tmp = new bool[times];
            for (int i = 0; i < times; i++)
            {
                tmp[i] = value;
            }

            return tmp;
        }

        public static double[] Sqr(this double[] x)
        {
            return Apply(x, Square);
        }

        public static double[] Sqrt(this double[] x)
        {
            return Apply(x, Math.Sqrt);
        }

        public static double StandardDev(this double[] x, bool skipNaN = false)
        {
            return x.MeanStandardDeviation(skipNaN).Item2;
        }


        public static double StandardDeviation(this IEnumerable<double> source)
        {
            return (double)Math.Sqrt((double)source.Variance());
        }


        public static Tuple<double, double> MeanStandardDeviation(this double[] x, bool skipNaN = false)
        {
            double mean = x.Mean(skipNaN);

            if (double.IsNaN(mean))
                return new Tuple<double, double>(Tools.NA, Tools.NA);

            double somme = 0;

            for (int i = 0; i < x.Length; i++)
            {
                double tmp = x[i] - mean;
                somme += (tmp * tmp);
            }

            return new Tuple<double, double>(mean, Math.Sqrt(somme / (x.Length - 1)));
        }

        public static string Show(this double[] x)
        {
            if (x == null)
            {
                return "null";
            }
            if (x.Length == 0)
            {
                return "()";
            }
            string[] zz = x.Select(z => z.ToString("R")).ToArray();
            return "(" + string.Join(",", zz) + ")";
        }

        public static string Show(this int[] x)
        {
            if (x == null)
            {
                return "null";
            }
            if (x.Length == 0)
            {
                return "()";
            }

            return "(" + string.Join(", ", x) + ")";
        }

        public static string ShowR(this double[] x)
        {
            if (x == null)
            {
                return "null";
            }

            if (x.Length == 0)
            {
                return "c()";
            }

            string[] zz = x.Select(z => z.ToString("R")).ToArray();
            return "c(" + string.Join(",", zz) + ")";
        }

        public static string ShowR(this double x)
        {
            return x.ToString("R");
        }

        public static double[] SubArray(this double[] x, int start, int length)
        {
            // Return new array.
            if (x == null)
                return null;
            if (x.Length == 0)
            {
                return new double[0];
            }

            if (start + length > x.Length)
            {
                length = x.Length - start;
            }
            double[] res = new double[length];

            for (int i = start, j = 0; j < length; i++, j++)
            {
                res[j] = x[i];
            }
            return res;
        }

        public static double[] Roots(this Complex[] roots, double tolerance = 1E-6, double lowerBound = double.NegativeInfinity, double upperBound = double.PositiveInfinity)
        {
            return roots.Select(c => c.Real).Where(r => (Math.Abs(r) < tolerance) && (r >= lowerBound) && (r <= upperBound)).ToArray();
        }

        public static double[] SubArray(this double[] x, int[] positions)
        {
            if (x == null)
                return null;
            if (x.Length == 0)
            {
                return new double[0];
            }

            if (positions.Any(i => i < 0))
            {
                throw new ArgumentOutOfRangeException(nameof(SubArray));
            }
            List<double> l = new List<double>();
            int previousPosition = -1;
            foreach (int p in positions.OrderBy(s => s))
            {
                if (p != previousPosition)
                {
                    l.Add(x[p]);
                    previousPosition = p;
                }
            }
            return l.ToArray();
        }

        public static double[] Substract(this double[] x, double[] y)
        {
            return DoNonCommutativeBinaryOperation(x, y, Substraction);
        }

        public static double[] Substract(this double[] x, double y)
        {
            return DoNonCommutativeBinaryOperation(x, new double[1] { y }, Substraction);
        }

        public static double[] Substract(this double x, double[] y)
        {
            return DoNonCommutativeBinaryOperation(new double[1] { x }, y, Substraction);
        }

        public static int[] ToInt32(this double[] x)
        {
            int[] rep = new int[x.Length];
            for (int i = 0; i < x.Length; i++)
                rep[i] = (Int32)x[i];
            return rep;
        }

        public static string ToR(this double[] x)
        {
            if (x == null)
            {
                return "NULL";
            }

            if (x.Length == 0)
            {
                return "numeric(0)";
            }

            StringBuilder sb = new StringBuilder();
            sb.Append("c(");
            sb.Append(x[0]);
            for (int i = 1; i < x.Length; i++)
            {
                sb.Append(", ").Append(x[i]);
            }

            sb.Append(")");
            return sb.ToString();
        }

        public static string ToString(double x)
        {
            if (IsFinite(x) || x == double.PositiveInfinity || x == double.NegativeInfinity)
            {
                return x.ToString();
            }
            if (Tools.IsNA(x))
            {
                return "NA";
            }
            return "?";
        }

        public static string ToString(double[] x)
        {
            const int len = 5;
            if (x == null)
            {
                return "null";
            }
            if (x.Length == 0)
            {
                return "[]";
            }

            int l = Math.Min(len, x.Length);
            string[] a = new string[l];
            for (int i = 0; i < l; i++)
            {
                a[i] = ToString(x[i]);
            }

            string rep = string.Join(", ", a);
            rep = "[" + rep + ((x.Length > l) ? " ...]" : "]");
            return rep;
        }

        public static double[] Unique(this double[] x)
        {
            // Probablement un peu couteux
            if (x == null || x.Length == 0)
                return new double[0];
            SortedSet<double> ss = new SortedSet<double>();
            List<double> ll = new List<double>();
            foreach (double x1 in x)
            {
                if (!ss.Contains(x1))
                {
                    ll.Add(x1);
                    ss.Add(x1);
                }
            }
            return ll.ToArray();
        }

        public static double Variance(this IEnumerable<double> source)
        {
            if (source == null)
            {
                throw new ArgumentNullException(NameOf.GetCallerName() + ": source");
            }

            long n = 0;
            double mean = 0;
            double M2 = 0;

            checked
            {
                foreach (double x in source)
                {
                    n++;
                    double delta = x - mean;
                    mean += delta / n;
                    M2 += delta * (x - mean);
                }
            }

            if (source.Count() < 2)
            {
                if (source.Count() == 1)
                {
                    return double.NaN;
                }

                //TODO : Exception
                throw new InvalidOperationException(NameOf.GetCallerName() + ": source contains no elements.");
            }

            return (double)(M2 / (n - 1));
        }

        public static int WhichMax(this double[] x)
        {
            double max = double.NegativeInfinity;
            int ndx = -1;
            for (int i = 0; i < x.Length; i++)
            {
                if (x[i] > max)
                {
                    max = x[i];
                    ndx = i;
                }
            }
            return ndx;
        }

        public static int WhichMin(this double[] x)
        {
            double min = double.PositiveInfinity;
            int ndx = -1;
            for (int i = 0; i < x.Length; i++)
            {
                if (x[i] < min)
                {
                    min = x[i];
                    ndx = i;
                }
            }
            return ndx;
        }


        //Opérations commutatives

        private static Func<double, double, double> Addition = delegate (double x, double y) { return x + y; };
        private static Func<double, double, double> Multiplication = delegate (double x, double y) { return x * y; };
        //Opérations non commutatives
        private static Func<double, double, double> Substraction = delegate (double x, double y) { return x - y; };
        private static Func<double, double, double> Division = delegate (double x, double y) { return x / y; };
        //Square
        private static Func<double, double> Square = delegate (double x) { return x * x; };

        private static double[] DoBinaryOperation(double[] x, double[] y, Func<double, double, double> operation)
        {
            if ((x == null) || (y == null))
                return new double[0];
            if ((x?.Length == 0) || (y?.Length == 0))
                return new double[0];
            int lx = x.Length, ly = y.Length;
            if (lx != ly)
            {
                bool multiple = lx > ly ? lx % ly == 0 : ly % lx == 0;
                if (!multiple)
                    return new double[0];
            }

            double[] lngst;
            double[] shrst;
            if (lx >= ly)
            {
                lngst = new double[lx];
                Array.Copy(x, lngst, lx);
                shrst = y;
            }
            else
            {
                lngst = new double[ly];
                Array.Copy(y, lngst, ly);
                shrst = x;
            }
            int iln, ish;
            for (iln = 0, ish = 0; iln < lngst.Length; iln++, ish = (ish + 1) % shrst.Length)
            {
                lngst[iln] = operation(lngst[iln], shrst[ish]);
            }
            return lngst;
        }


        private static double[] DoNonCommutativeBinaryOperation(double[] x, double[] y, Func<double, double, double> operation)
        {
            if ((x == null) || (y == null))
                return new double[0];
            if ((x?.Length == 0) || (y?.Length == 0))
                return new double[0];
            int lx = x.Length, ly = y.Length;
            if (lx < ly)
                return new double[0];

            if ((lx != ly) && (lx % ly != 0))
                return new double[0];

            double[] lngst;
            double[] shrst;
            lngst = new double[lx];
            Array.Copy(x, lngst, lx);
            shrst = y;
            int iln, ish;

            for (iln = 0, ish = 0; iln < lngst.Length; iln++, ish = (ish + 1) % shrst.Length)
            {
                lngst[iln] = operation(lngst[iln], shrst[ish]);
            }
            return lngst;
        }

        // traitement des nombres flottants

        public static int BiasedExponent(double x)
        {
            return (int)((BitConverter.DoubleToInt64Bits(x) & _exponentMask) >> 52);
        }

        public static double scalbln(double number, long exponent)
        {
            const long DBL_EXP_MASK = 0x7ff0000000000000L;
            const int DBL_MANT_BITS = 52;
            const long DBL_MANT_MASK = 0x000fffffffffffffL;
            const long DBL_SGN_MASK = -1 - 0x7fffffffffffffffL;
            const long DBL_EXP_CLR_MASK = DBL_SGN_MASK | DBL_MANT_MASK;

            long bits = System.BitConverter.DoubleToInt64Bits(number);
            int exp = (int)((bits & DBL_EXP_MASK) >> DBL_MANT_BITS);
            // Check for infinity or NaN.
            if (exp == 0x7ff)
                return number;
            // Check for 0 or subnormal.
            if (exp == 0)
            {
                // Check for 0.
                if ((bits & DBL_MANT_MASK) == 0)
                    return number;
                // Subnormal, scale number so that it is in [1, 2).
                number *= System.BitConverter.Int64BitsToDouble(0x4350000000000000L); // 2^54
                bits = System.BitConverter.DoubleToInt64Bits(number);
                exp = (int)((bits & DBL_EXP_MASK) >> DBL_MANT_BITS) - 54;
            }
            // Check for underflow.
            if (exponent < -50000)
                return copysign(0D, number);
            // Check for overflow.
            if (exponent > 50000 || (long)exp + exponent > 0x7feL)
                return copysign(System.Double.PositiveInfinity, number);
            exp += (int)exponent;
            // Check for normal.
            if (exp > 0)
                return System.BitConverter.Int64BitsToDouble((bits & DBL_EXP_CLR_MASK) | ((long)exp << DBL_MANT_BITS));
            // Check for underflow.
            if (exp <= -54)
                return copysign(0D, number);
            // Subnormal.
            exp += 54;
            number = System.BitConverter.Int64BitsToDouble((bits & DBL_EXP_CLR_MASK) | ((long)exp << DBL_MANT_BITS));
            return number * System.BitConverter.Int64BitsToDouble(0x3c90000000000000L); // 2^-54
        }

        private static int signbit(double number)
        {
            //public const long DBL_SGN_MASK = -1 - 0x7fffffffffffffffL;
            const long DBL_SGN_MASK = -1 - 0x7fffffffffffffffL;
            if (System.Double.IsNaN(number))
                return ((System.BitConverter.DoubleToInt64Bits(number) & DBL_SGN_MASK) != 0) ? 0 : 1;
            else
                return ((System.BitConverter.DoubleToInt64Bits(number) & DBL_SGN_MASK) != 0) ? 1 : 0;
        }

        private static double copysign(double number1, double number2)
        {
            const long DBL_SGN_MASK = -1 - 0x7fffffffffffffffL;
            const long DBL_SGN_CLR_MASK = 0x7fffffffffffffffL;
            // If number1 is NaN, we have to store in it the opposite of the sign bit.
            long sign = (signbit(number2) == 1 ? DBL_SGN_MASK : 0L) ^ (System.Double.IsNaN(number1) ? DBL_SGN_MASK : 0L);
            return System.BitConverter.Int64BitsToDouble((System.BitConverter.DoubleToInt64Bits(number1) & DBL_SGN_CLR_MASK) | sign);
        }

        internal static double Scalb(double x, int n)
        {
            // Treat special cases first.
            if (x == 0 || double.IsInfinity(x) || double.IsNaN(x))
                return x;
            if (n == 0)
                return x;

            int e = BiasedExponent(x);
            long mantissa = Mantissa(x);
            long sign = ((x > 0) ? 0 : _signMask);

            // Is x denormalized?
            if (e == 0)
            {
                if (n < 0)
                {
                    // n negative means we have to shift the mantissa -n bits to the right.
                    mantissa >>= -n;
                    return BitConverter.Int64BitsToDouble(sign | mantissa);
                }
                else
                {
                    // n positive means we need to shift to the left until we get a normalized number...
                    // if we get there, that is.
                    while (mantissa <= _mantissaMask && n > 0)
                    {
                        mantissa <<= 1;
                        n--;
                    }
                    if (mantissa > _mantissaMask)
                        n++;
                    // The value of n is now the biased exponent.

                    // Does the result overflow?
                    if (n > 2 * _bias)
                        return (x > 0) ? double.PositiveInfinity : double.NegativeInfinity;

                    // n is the biased exponent of the result.
                    return BitConverter.Int64BitsToDouble(sign | ((long)n << 52) | (mantissa & _mantissaMask));
                }
            }

            // If we get here, we know x is normalized.
            // Do scaling. e becomes the biased exponent of the result.
            e = e + n;

            // Check for 0 or denormalized.
            if (e < 0)
            {
                mantissa = ((1L << 52) + mantissa) >> (1 - e);
                return BitConverter.Int64BitsToDouble(sign | mantissa);
            }

            // Check for overflow.
            if (e > 2 * _bias)
                return (x > 0) ? double.PositiveInfinity : double.NegativeInfinity;

            // If we're here, the result is normalized.
            long bits = sign | ((long)e << 52) | mantissa;
            return BitConverter.Int64BitsToDouble(bits);
        }

        private const long _exponentMask = 0x7FF0000000000000;
        private const long _signMask = -1 - 0x7FFFFFFFFFFFFFFF;
        private const long _mantissaMask = 0xFFFFFFFFFFFFF;
        private const int _bias = 1023;

        // Complexes

        internal static double[] Real(this Complex[] x)
        {
            if (x == null)
            {
                return null;
            }

            if (x.Length == 0)
            {
                return new double[0];
            }

            double[] rep = new double[x.Length];
            for (int i = 0; i < x.Length; i++)
            {
                rep[i] = x[i].Real;
            }
            return rep;
        }

        internal static double[] Imaginary(this Complex[] x)
        {
            if (x == null)
            {
                return null;
            }

            if (x.Length == 0)
            {
                return new double[0];
            }

            double[] rep = new double[x.Length];
            for (int i = 0; i < x.Length; i++)
            {
                rep[i] = x[i].Imaginary;
            }
            return rep;
        }

        internal static string ShowR(this Complex x)
        {
            return string.Format("{0:R} + {1:R}i", x.Real, x.Imaginary);
        }

        internal static string ShowR(this Complex[] x)
        {
            if (x == null)
            {
                return "null";
            }

            if (x.Length == 0)
            {
                return "as.complex()";
            }

            string[] zz = x.Select(z => ShowR(z)).ToArray();
            return "c(" + string.Join(",", zz) + ")";
        }

    }// class
}