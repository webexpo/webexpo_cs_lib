namespace Zygotine.WebExpo
{
    using System;
    using System.Collections.Generic;
    using System.Linq;
    using Zygotine.Util;

    public abstract class AbstractBWModelParameters : AbstractModelParameters
    {
        public bool UseUniformPriorOnSds { get; set; }
        public double MuOverallLower { get; set; } = Tools.ND;
        public double MuOverallUpper { get; set; } = Tools.ND;
        public double InitMuOverall { get; set; } = Tools.ND;
        public double InitSigmaWithin { get; set; } = Tools.ND;

        internal AbstractBWModelParameters(bool useUniformPriorOnSds, bool logNormalDistrn) : base(logNormalDistrn)
        {
            this.UseUniformPriorOnSds = useUniformPriorOnSds;
        }

        internal void SetMuOverallLower(double value)
        {
            this.MuOverallLower = value;
        }

        internal void SetMuOverallUpper(double value)
        {
            this.MuOverallUpper = value;
        }

        internal void SetInitMuOverall(double value)
        {
            this.InitMuOverall = value;
        }

        internal void SetInitSigmaWithin(double value)
        {
            this.InitSigmaWithin = value;
        }

        internal InitialValues DefaultInits(DataSummary data, bool includeCensoredData)
        {
            double muOverall = Tools.ND;
            double sigmaWithin = Tools.ND;
            bool logNormalDstrn = this.LogNormalDstrn;
            bool useUniformPriorOnSds = this.UseUniformPriorOnSds;
            

            List<double> workerMeans = new List<double>();
            List<double> workerSds = new List<double>();
            List<double> normalizedData = new List<double>();
            int[] MeasureCountByWorker = new int[data.NWorkers];
            Dictionary<int, List<double>> measureValueByWorker = new Dictionary<int, List<double>>();

            if (data.YLength > 0)
            {
                //Une liste pour chaque travailleurs. wo == workerOrdinal
                for (int wo = 0; wo < data.NWorkers; wo++)
                {
                    measureValueByWorker[wo] = new List<double>();
                }

                Measure[] lm;
                //On ajoute les mesures pour chacun des travailleurs
                lm = data.UncensoredMeasures;
                for (int i = 0; i < lm.Length; i++)
                {
                    int wo = lm[i].Worker.Ordinal;
                    List<double> lst = measureValueByWorker[wo];
                    double value = logNormalDstrn ? Math.Log(lm[i].A) : lm[i].A;
                    lst.Add(value);
                    normalizedData.Add(value);
                }
            }

            if (includeCensoredData)
            {
                Measure[] lm;
                lm = data.GTMeasures;
                for (int i = 0; i < lm.Length; i++)
                {
                    int wo = lm[i].Worker.Ordinal;
                    List<double> lst = measureValueByWorker[wo];
                    double value = logNormalDstrn ? Math.Log(lm[i].A) : lm[i].A;
                    lst.Add(value);
                    normalizedData.Add(value);
                }

                lm = data.LTMeasures;
                for (int i = 0; i < lm.Length; i++)
                {
                    int wo = lm[i].Worker.Ordinal;
                    List<double> lst = measureValueByWorker[wo];
                    double value = logNormalDstrn ? Math.Log(lm[i].A) : lm[i].A;
                    lst.Add(value);
                    normalizedData.Add(value);
                }

                lm = data.IntervalMeasures;
                for (int i = 0; i < lm.Length; i++)
                {
                    int wo = lm[i].Worker.Ordinal;
                    List<double> lst = measureValueByWorker[wo];
                    double midP = lm[i].A + (lm[i].B - lm[i].A) / 2;
                    double value = logNormalDstrn ? Math.Log(midP) : midP;
                    lst.Add(value);
                    normalizedData.Add(value);
                }
            }

            if (normalizedData.Count > 0)
            {
                for (int o = 0; o < data.NWorkers; o++)
                {
                    List<double> src = measureValueByWorker[o];
                    if (src.Count == 0)
                    {
                        continue;
                    }

                    if (src.Count == 1)
                    {
                        workerMeans.Add(src[0]);
                        continue; // pas de sd
                    }

                    int n = 0;
                    double mean = 0;
                    double M2 = 0;

                    foreach (double x in src)
                    {
                        n++;
                        double delta = x - mean;
                        mean += delta / n;
                        M2 += delta * (x - mean);
                    }

                    workerMeans.Add(mean);
                    if (n > 1)
                    {
                        workerSds.Add(Math.Sqrt(M2 / (n - 1)));
                    }
                }

                //means.Length > 0
                muOverall = workerMeans.Average();

                //# modif_0.13
                if (workerSds.Count == 0)
                {
                    if (workerMeans.Count > 1)
                    {
                        sigmaWithin = Math.Sqrt(workerMeans.Variance());
                    }
                }
                else
                {
                    sigmaWithin = workerSds.Mean();
                }
            }
            else
            {
                muOverall = 0.0;
                normalizedData.Clear();
            }

            // # new_0.13
            if (Tools.IsND(sigmaWithin) && !includeCensoredData)
            {
                //# call this fct again, but this time including censored data
                InitialValues tmp = DefaultInits(data, includeCensoredData: true);
                sigmaWithin = tmp.SigmaWithin;
            }

            // # new_0.13
            if (Tools.IsND(sigmaWithin))
            {
                sigmaWithin = GetSigmaWithin();
            }

            return new InitialValues(muOverall, sigmaWithin, normalizedData);
        } //# end of Default.inits.local

        protected abstract double GetSigmaWithin();
    }
}
