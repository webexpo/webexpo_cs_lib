namespace Zygotine.WebExpo
{
    using System.Linq;
    using Zygotine.Util;

    internal class OutLogoutMoments
    {
        public double Sum { get; private set; }
        public double Sum2 { get; private set; }

        private OutLogoutMoments()
        {
        }

        public static OutLogoutMoments Get(bool anyME, bool logNormalDistrn, DataSummary data, YGen genY, TrueValuesGen genTV)
        {
            OutLogoutMoments olm = new OutLogoutMoments();
            if (anyME)
            {
                double[] y;
                if (logNormalDistrn)
                {
                    y = Tools.Combine(genTV.LogY, genTV.LogGT, genTV.LogLT, genTV.LogI);
                }
                else
                {
                    y = Tools.Combine(genTV.Y, genTV.GT, genTV.LT, genTV.I);
                }

                olm.Sum = y.Sum();
                olm.Sum2 = y.Sqr().Sum();
                return olm;
            }
            else if (logNormalDistrn)
            {
                olm.Sum = data.LogUncensoredSum;
                olm.Sum2 = data.LogUncensoredSum2;
                olm.Sum += genY.LogGT.Sum() + genY.LogLT.Sum() + genY.LogI.Sum();
                olm.Sum2 += genY.LogGT.Sqr().Sum() + genY.LogLT.Sqr().Sum() + genY.LogI.Sqr().Sum();
            }
            else
            {
                olm.Sum = data.UncensoredSum;
                olm.Sum2 = data.UncensoredSum2;
                olm.Sum += genY.GT.Sum() + genY.LT.Sum() + genY.I.Sum();
                olm.Sum2 += genY.GT.Sqr().Sum() + genY.LT.Sqr().Sum() + genY.I.Sqr().Sum();
            }

            return olm;
        }
    }
}
