namespace Zygotine.WebExpo
{
    internal class MuTruncatedGen
    {
        public double Mu { get; private set; }

        private MuTruncatedGen()
        {
        }

        private MuTruncatedGen(GenObject o, double[] range, double muMean, double sigma)
        {
            SGNFnAParam A = o.A.Clone();
            A.MuMean = muMean;
            A.S = sigma;
            A.S2 = sigma * sigma;
            Icdf icdf = new Icdf(o, A, range);
            this.Mu = icdf.Bidon();
        }

        public static MuTruncatedGen GetInstance(GenObject o, double[] range, double muMean, double sigma)
        {
            return new MuTruncatedGen(o, range, muMean, sigma);
        }
    }
}
