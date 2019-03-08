namespace Zygotine.Util
{
    using System.Collections.Generic;
    using System.Numerics;

    internal class DescendingComparer : Comparer<double>
    {
        private DescendingComparer() { }

        public override int Compare(double x, double y)
        {
            return y.CompareTo(x);
        }

        internal static readonly Comparer<double> Desc = new DescendingComparer();
    }


    public class ComplexComparer : IComparer<Complex>
    {
        private ComplexComparer() { }

        public int Compare(Complex x, Complex y)
        {
            switch (x.Real.CompareTo(y.Real))
            {
                case -1: return -1;
                case 1: return 1;
                default: return x.Imaginary.CompareTo(y.Imaginary);
            }
        }

        public static readonly IComparer<Complex> Desc = new ComplexComparer();
    }


    public class DoubleComparer : IComparer<double>
    {
        private DoubleComparer() { }

        public int Compare(double x, double y)
        {
            return -x.CompareTo(y);
        }

        public static readonly IComparer<double> Desc = new DoubleComparer();
    }

}
