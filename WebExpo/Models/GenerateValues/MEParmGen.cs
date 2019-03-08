

namespace Zygotine.WebExpo
{
    using System.Linq;
    using Zygotine.Util;

    internal class MEParmGen
    {

        public double Parm { get; private set; }

        private MEParmGen()
        {
        }

        internal static MEParmGen GetInstance(GenObject o, MeasurementError me, DataSummary data, YGen genY, TrueValuesGen genTV)
        {

            MEParmGen instance = new MEParmGen();
            double b;
            double[] tmpY, tmpT;

            if (me.ThroughCV)
            {

                if (o.LogNormalDistrn)
                {
                    tmpY = Tools.Combine(data.LogY, genY.LogGT, genY.LogLT, genY.LogI);
                    tmpT = Tools.Combine(genTV.LogY, genTV.LogGT, genTV.LogLT, genTV.LogI);
                    b = tmpY.Substract(tmpT).Exp().Substract(1.0).Sqr().Sum();
                    b /= 2.0;
                }
                else
                {
                    tmpY = Tools.Combine(data.Y, genY.GT, genY.LT, genY.I);
                    tmpT = Tools.Combine(genTV.Y, genTV.GT, genTV.LT, genTV.I);
                    b = tmpY.Divide(tmpT).Substract(1.0).Sqr().Reverse().Sum();
                    b /= 2.0;
                }

                SGNFnAParam localA = o.A.Clone();
                localA.B = b;
                localA.Range = me.GetRange();
                double[] range = me.GetRange();
                Icdf icdf = new Icdf(o, localA, range);
                instance.Parm = icdf.Bidon(o.Start(localA), inestLowerLim: range[0] == 0);
            }
            else
            {
                tmpY = Tools.Combine(data.Y, genY.GT, genY.LT, genY.I);
                tmpT = Tools.Combine(genTV.Y, genTV.GT, genTV.LT, genTV.I);
                b = tmpY.Substract(tmpT).Sqr().Sum();
                b /= 2.0;
                if (o.LogNormalDistrn)
                {
                    SGNFnAParam localA = o.A.Clone();
                    localA.B = b;
                    localA.Range = me.GetRange();
                    localA.Truevalues = Tools.Copy(tmpT);
                    //me.parm <- dens.gen.icdf(o, A, range=me$range, inestimable.lower.limit=me$range[1]==0)
                    double[] range = me.GetRange();
                    Icdf icdf = new Icdf(o, localA, range);
                    instance.Parm = icdf.Bidon(inestLowerLim: range[0] == 0.0);
                }
                else
                {
                    instance.Parm = WebExpoFunctions3.SqrtInvertedGammaGen(data.N, b, me.GetRange(), o);
                }
            }

            return instance;
        }

    }
}
