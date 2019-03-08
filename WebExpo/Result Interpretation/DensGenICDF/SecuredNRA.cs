namespace Zygotine.WebExpo
{
    using System;

    internal interface ISecuredNR
    {
        double[] Range { get; set; }
        Func<double, SGNFnAParam, double> H { get; set; }
        Func<double, SGNFnAParam, double> HPrime { get; set; }
        Func<double, SGNFnAParam, double> HSecond { get; set; }
    }

    internal class SecuredNRA : ISecuredNR
    {
        public double[] Range { get; set; }

        public Func<double, SGNFnAParam, double> H { get; set; } = null;
        public Func<double, SGNFnAParam, double> HPrime { get; set; } = null;
        public Func<double, SGNFnAParam, double> HSecond { get; set; } = null;

        internal SecuredNRA(Func<double, SGNFnAParam, double> h, Func<double, SGNFnAParam, double> hp)
        {
            this.H = h;
            this.HPrime = hp;
        }

        internal SecuredNRA(Func<double, SGNFnAParam, double> h, Func<double, SGNFnAParam, double> hp, Func<double, SGNFnAParam, double> hs)
        {
            this.H = h;
            this.HPrime = hp;
            this.HSecond = hs;
        }

        internal SecuredNRA(Func<double, SGNFnAParam, double> h, Func<double, SGNFnAParam, double> hp, double[] range) : this(h, hp)
        {
            this.Range = range;
        }

    }
}
