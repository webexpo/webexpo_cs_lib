namespace Zygotine.WebExpo
{
    using System.Text;
    using Zygotine.Util;

    internal class YGen
    {
        internal double[] LogGT { get; set; } = new double[0];
        internal double[] GT { get; set; } = new double[0];
        internal double[] LogLT { get; set; } = new double[0];
        internal double[] LT { get; set; } = new double[0];
        internal double[] LogI { get; set; } = new double[0];
        internal double[] I = new double[0];


        private YGen()
        {
        }

        public override string ToString()
        {
            string fmt = "GT={0}, LT={1}, I={2}";
            string fmt2 = "{0}, **LOG: {1}";
            StringBuilder sb = new StringBuilder();
            sb.AppendFormat(fmt2, string.Format(fmt, GT.ShowR(), LT.ShowR(), I.ShowR()), string.Format(fmt, LogGT.ShowR(), LogLT.ShowR(), LogLT.ShowR()));
            return sb.ToString();
        }

        private readonly static YGen emptyInstance = new YGen();
        public static YGen EmptyInstance { get { return emptyInstance; } }
        private YGen(MeasurementError me, TrueValuesGen genTV, DataSummary data, double mu, double sigma, bool logNormalDistrn = true)
        {
            if (me.Any)
            {
                // Sample y (censored) values | true values
                if (data.AnyGT)
                {
                    double[] tmpMean = genTV.GT;
                    double[] tmpSD;
                    if (me.ThroughCV)
                    {
                        tmpSD = tmpMean.Multiply(me.Parm);
                    }
                    else
                    {
                        tmpSD = Tools.Rep(me.Parm, tmpMean.Length);
                    }

                    this.GT = RNorm4CensoredMeasures.RNormCensored(tmpMean, tmpSD, lower: data.GT);
                    if (logNormalDistrn)
                    {
                        this.LogGT = this.GT.Log();
                    }
                }

                if (data.AnyLT)
                {
                    double[] tmpMean = genTV.LT;
                    double[] tmpSD;
                    if (me.ThroughCV)
                    {
                        tmpSD = tmpMean.Multiply(me.Parm);
                    }
                    else
                    {
                        tmpSD = Tools.Rep(me.Parm, tmpMean.Length);
                    }

                    this.I = RNorm4CensoredMeasures.RNormCensored(tmpMean, tmpSD, upper: data.LT, negativeDisallowed: logNormalDistrn || me.ThroughCV);
                    if (logNormalDistrn)
                    {
                        this.LogLT = this.LT.Log();
                    }
                }

                if (data.AnyIntervalCensored)
                {
                    double[] tmpMean = genTV.I;
                    double[] tmpSD;
                    if (me.ThroughCV)
                    {
                        tmpSD = tmpMean.Multiply(me.Parm);
                    }
                    else
                    {
                        tmpSD = Tools.Rep(me.Parm, tmpMean.Length);
                    }

                    this.I = RNorm4CensoredMeasures.RNormCensored(tmpMean, tmpSD, lower: data.IntervalGT, upper: data.IntervalLT);
                    if (logNormalDistrn)
                    {
                        this.LogI = this.I.Log();
                    }
                }
            }
            else
            {
                if (logNormalDistrn)
                {
                    if (data.AnyGT)
                    {
                        this.LogGT = RNorm4CensoredMeasures.RNormCensored(mu, sigma, lower: data.LogGT);
                        this.GT = this.LogGT.Exp();
                    }

                    if (data.AnyLT)
                    {
                        this.LogLT = RNorm4CensoredMeasures.RNormCensored(mu, sigma, upper: data.LogLT);
                        this.LT = this.LogLT.Exp();
                    }

                    if (data.AnyIntervalCensored)
                    {
                        this.LogI = RNorm4CensoredMeasures.RNormCensored(mu, sigma, lower: data.LogIntervalGT, upper: data.LogIntervalLT);
                        this.I = this.LogI.Exp();
                    }
                }
                else
                {
                    if (data.AnyGT)
                    {
                        this.GT = RNorm4CensoredMeasures.RNormCensored(mu, sigma, lower: data.GT);
                    }

                    if (data.AnyLT)
                    {
                        this.LT = RNorm4CensoredMeasures.RNormCensored(mu, sigma, upper: data.LT);
                    }

                    if (data.AnyIntervalCensored)
                    {
                        this.I = RNorm4CensoredMeasures.RNormCensored(mu, sigma, lower: data.IntervalGT, upper: data.IntervalLT);
                    }
                }
            }
        }// internal YGen (constructor)

        public static YGen GetEmptyObject()
        {
            return new YGen();
        }

        public static YGen GetInstance(MeasurementError me, TrueValuesGen genTV, DataSummary data, double mu, double sigma, bool logNormalDistrn = true)
        {

            if (me.Any)
            {
                YGen newYGen = new YGen();
                // Sample y (censored) values | true values
                if (data.AnyGT)
                {
                    double[] tmpMean = genTV.GT;
                    double[] tmpSD;
                    if (me.ThroughCV)
                    {
                        tmpSD = tmpMean.Multiply(me.Parm);
                    }
                    else
                    {
                        tmpSD = Tools.Rep(me.Parm, tmpMean.Length);
                    }

                    newYGen.GT = RNorm4CensoredMeasures.RNormCensored(tmpMean, tmpSD, lower: data.GT);
                    if (logNormalDistrn)
                    {
                        newYGen.LogGT = newYGen.GT.Log();
                    }
                }

                if (data.AnyLT)
                {
                    double[] tmpMean = genTV.LT;
                    double[] tmpSD;
                    if (me.ThroughCV)
                    {
                        tmpSD = tmpMean.Multiply(me.Parm);
                    }
                    else
                    {
                        tmpSD = Tools.Rep(me.Parm, tmpMean.Length);
                    }

                    newYGen.LT = RNorm4CensoredMeasures.RNormCensored(tmpMean, tmpSD, upper: data.LT, negativeDisallowed: logNormalDistrn || me.ThroughCV);
                    if (logNormalDistrn)
                    {
                        newYGen.LogLT = newYGen.LT.Log();
                    }
                }

                if (data.AnyIntervalCensored)
                {
                    double[] tmpMean = genTV.I;
                    double[] tmpSD;
                    if (me.ThroughCV)
                    {
                        tmpSD = tmpMean.Multiply(me.Parm);
                    }
                    else
                    {
                        tmpSD = Tools.Rep(me.Parm, tmpMean.Length);
                    }

                    newYGen.I = RNorm4CensoredMeasures.RNormCensored(tmpMean, tmpSD, lower: data.IntervalGT, upper: data.IntervalLT);
                    if (logNormalDistrn)
                    {
                        newYGen.LogI = newYGen.I.Log();
                    }
                }
                return newYGen;
            }
            else
            {
                return YGen.Inits(data, mu, sigma, meThroughCV: data.METhroughCV, logNormalDistrn: logNormalDistrn);
            }
        }

        internal static YGen Inits( DataSummary data, double mu, double sigma, bool meThroughCV, bool logNormalDistrn)
        {
            YGen newYGen = new YGen();

            if (logNormalDistrn)
            {
                if (data.AnyGT)
                {
                    newYGen.LogGT = RNorm4CensoredMeasures.RNormCensored(mu, sigma, lower: data.LogGT);
                    newYGen.GT = newYGen.LogGT.Exp();
                }

                if (data.AnyLT)
                {
                    newYGen.LogLT = RNorm4CensoredMeasures.RNormCensored(mu, sigma, upper: data.LogLT, negativeDisallowed: meThroughCV && !logNormalDistrn);
                    newYGen.LT = newYGen.LogLT.Exp();
                }

                if (data.AnyIntervalCensored)
                {
                    newYGen.LogI = RNorm4CensoredMeasures.RNormCensored(mu, sigma, lower: data.LogIntervalGT, upper: data.LogIntervalLT);
                    newYGen.I = newYGen.LogI.Exp();
                }
            }
            else
            {
                if (data.AnyGT)
                {
                    newYGen.GT = RNorm4CensoredMeasures.RNormCensored(mu, sigma, lower: data.GT);
                }

                if (data.AnyLT)
                {
                    newYGen.LT = RNorm4CensoredMeasures.RNormCensored(mu, sigma, upper: data.LT, negativeDisallowed: meThroughCV && !logNormalDistrn); 
                }

                if (data.AnyIntervalCensored)
                {
                    newYGen.I = RNorm4CensoredMeasures.RNormCensored(mu, sigma, lower: data.IntervalGT, upper: data.IntervalLT);
                }
            }

            return newYGen;
        }
    }
}