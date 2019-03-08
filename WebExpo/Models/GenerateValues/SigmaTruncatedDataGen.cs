namespace Zygotine.WebExpo
{
    using Zygotine.Util;
    internal class SigmaTruncatedDataGen
    {
        public double Sigma { get; private set; }

        private SigmaTruncatedDataGen()
        {
        }

        public static SigmaTruncatedDataGen GetInstance(GenObject o, double[] range, double b, double mu, double currentSigma)
        {
            SigmaTruncatedDataGen rep = new SigmaTruncatedDataGen();
            SGNFnAParam localA = o.A.Clone();
            localA.B = b;
            localA.Mu = mu;
            double[] start = o.Start(localA);
            if(start.Length == 0)
            {
                start = Tools.Combine(currentSigma);
            }

            Icdf icdf = new Icdf(o, localA, range);
            rep.Sigma = icdf.Bidon(start, inestLowerLim: range[0] == 0);
            return rep;
        }

    }
}
