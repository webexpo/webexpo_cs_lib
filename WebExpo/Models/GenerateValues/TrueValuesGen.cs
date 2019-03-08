
namespace Zygotine.WebExpo
{
    using System;
    using System.Text;
    using Zygotine.Statistics.Distribution;
    using Zygotine.Util;

    internal class TrueValuesGen
    {
        internal double[] LogY = new double[0];
        internal double[] Y = new double[0];
        internal double[] LogGT { get; set; } = new double[0];
        internal double[] GT { get; set; } = new double[0];
        internal double[] LogLT { get; set; } = new double[0];
        internal double[] LT { get; set; } = new double[0];
        internal double[] LogI { get; set; } = new double[0];
        internal double[] I = new double[0];

        private TrueValuesGen()
        {
        }


        public string Show()
        {
            string fmt = "Y={0}, GT={1}, LT={2}, I={3}";
            string fmt2 = "{0}, **LOG: {1}";
            StringBuilder sb = new StringBuilder();
            sb.AppendFormat(fmt2, string.Format(fmt, Y.ShowR(), GT.ShowR(), LT.ShowR(), I.ShowR()), string.Format(fmt, LogY.ShowR(), LogGT.ShowR(), LogLT.ShowR(), LogLT.ShowR()));
            return sb.ToString();
        }

        private TrueValuesGen(YGen genY, DataSummary data, double mu, double sigma, MeasurementError me, bool logNormalDistrn = true, GenObject o = null)
        {
            if (me.ThroughSD && !logNormalDistrn)
            {
                double meSD = me.Parm;
                double tauStar = 1 / Math.Pow(sigma, 2) + 1 / Math.Pow(meSD, 2);
                double sdStar = 1 / Math.Sqrt(tauStar);
                if (data.YLength > 0)
                {
                    double[] tmpMean = (data.Y.Divide(Math.Pow(meSD, 2)).Add(mu / Math.Pow(sigma, 2))).Divide(tauStar);
                    this.Y = NormalDistribution.RNorm(data.YLength, tmpMean, Tools.Rep(sdStar, tmpMean.Length));
                }

                if (data.GTLength > 0)
                {
                    double[] tmpMean = (genY.GT.Divide(Math.Pow(meSD, 2)).Add(mu / Math.Pow(sigma, 2))).Divide(tauStar);
                    this.Y = NormalDistribution.RNorm(data.GTLength, tmpMean, Tools.Rep(sdStar, tmpMean.Length));
                }

                if (data.LTLength > 0)
                {
                    double[] tmpMean = (genY.LT.Divide(Math.Pow(meSD, 2)).Add(mu / Math.Pow(sigma, 2))).Divide(tauStar);
                    this.Y = NormalDistribution.RNorm(data.GTLength, tmpMean, Tools.Rep(sdStar, tmpMean.Length));
                }

                if (data.IntervalLength > 0)
                {
                    double[] tmpMean = (genY.I.Divide(Math.Pow(meSD, 2)).Add(mu / Math.Pow(sigma, 2))).Divide(tauStar);
                    this.Y = NormalDistribution.RNorm(data.GTLength, tmpMean, Tools.Rep(sdStar, tmpMean.Length));
                }
            }
            else
            {
                o.A.Sigma2 = sigma * sigma;
                this.Y = new double[data.YLength];
                this.GT = new double[data.GTLength];
                this.LT = new double[data.LTLength];
                this.I = new double[data.IntervalLength];

                for (int j = 0; j < data.YLength; j++)
                {
                    this.Y[j] = TrueValueGen(o, me, data.Y[j], mu, logY: logNormalDistrn ? data.LogY[j] : Tools.ND);
                }

                for (int j = 0; j < data.GTLength; j++)
                {
                    this.GT[j] = TrueValueGen(o, me, genY.GT[j], mu, logY: logNormalDistrn ? genY.LogGT[j] : Tools.ND);
                }

                for (int j = 0; j < data.LTLength; j++)
                {
                    this.LT[j] = TrueValueGen(o, me, genY.LT[j], mu, logY: logNormalDistrn ? genY.LogLT[j] : Tools.ND);
                }

                for (int j = 0; j < data.IntervalLength; j++)
                {
                    this.I[j] = TrueValueGen(o, me, genY.I[j], mu, logY: logNormalDistrn ? genY.LogI[j] : Tools.ND);
                }

                if (logNormalDistrn)
                {
                    this.LogY = this.Y.Log();
                    this.LogGT = this.GT.Log();
                    this.LogLT = this.LT.Log();
                    this.LogI = this.I.Log();
                }
            }
        } //# end of truevalues.gen   

        internal static TrueValuesGen GetInstance(YGen genY, DataSummary data, double mu, double sigma, MeasurementError me, bool logNormalDistrn = true, GenObject o = null)
        {
            return new TrueValuesGen(genY, data, mu, sigma, me, logNormalDistrn, o);
        }

        private double TrueValueGen(GenObject genO, MeasurementError me, double y, double mu, double logY)
        {
            SGNFnAParam localA = genO.A.Clone();
            localA.Mu = mu;
            if (genO.LogNormalDistrn)
            {
                localA.LogY = logY;
            }

            localA.Y = y;
            localA.Y2 = y * y;
            if (genO.ThroughSD)
            {
                localA.Ksi = me.Parm;
                localA.Ksi2 = me.Parm * me.Parm;
            }
            else
            {
                localA.CV2 = me.Parm * me.Parm;
            }

            double[] start = genO.Start(localA);
            //start.Any(zz => !zz.IsFinite());
            Icdf icdf = new Icdf(genO, localA, range: genO.Range);
            double x = icdf.Bidon(start, inestLowerLim: !genO.LogNormalDistrn);
            if (genO.LogNormalDistrn)
            {
                x = Math.Exp(x);// # new_0.11
            }

            return x;
        } //# end of truevalue.gen
    }
}
