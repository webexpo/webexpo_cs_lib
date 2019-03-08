namespace Zygotine.WebExpo
{
    using System;
    using Zygotine.Util;

    internal class Bounds
    {
        //#'bounds' is a list with the following elements/dimensions:
        //# x, h, hp, range:   numeric (vectors of length 2)
        //# higher.side:       numeric (1)
        //# include.soln:      logical (1)

        internal double[] H { get; set; }
        internal int HigherSide { get; set; }
        internal double[] Hp { get; set; }
        internal bool IncludeSoln { get; set; }
        internal double[] Range { get; set; }
        internal double[] X { get; set; }

        private Bounds()
        {
        }

        internal Bounds(double[] x, double[] h, double[] hp, double[] range, int higherSide, bool includeSoln)
        {
            this.X = x;
            this.H = h;
            this.Hp = hp;
            this.Range = range;
            this.HigherSide = higherSide;
            this.IncludeSoln = includeSoln;
        }


        internal Bounds Copy()
        {
            Bounds b = new Bounds();
            b.H = Tools.Copy(this.H);
            b.HigherSide = this.HigherSide;
            b.Hp = Tools.Copy(this.Hp);
            b.IncludeSoln = this.IncludeSoln;
            b.Range = Tools.Copy(this.Range);
            b.X = Tools.Copy(this.X);
            return b;
        }

        public override string ToString()
        {
            string s1 = string.Format("X={0}, H={1}, Hp={2}", X.Show(), H.Show(), Hp.Show());
            string s2 = string.Format(", Range={0}", Tools.Show(Range));
            string s3 = string.Format(", HigherSide={0}, IncludeSoln={1}", HigherSide, IncludeSoln);
            return s1 + s2 + s3;
        }
    }
}
