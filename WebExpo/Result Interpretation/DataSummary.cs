namespace Zygotine.WebExpo
{
    using System;
    using System.Collections.Generic;
    using System.Linq;
    using Zygotine.Util;

    public class DataSummary
    {
        internal class RawData
        {
            internal RawData()
            {
            }
            internal double[] Y;
            internal double[] LT;
            internal double[] GT;
            internal Tuple<double[], double[]> IntervalCensored;
            internal double UncensoredSum = double.NaN;
            internal double UncensoredSum2 = double.NaN;
        }

        private bool logNormalDistrn;
        private RawData logOfRaw = null;
        private RawData raw = new RawData();
        public int N { get; private set; }

        //booléen
        public bool MEAny { get; private set; }
        public bool METhroughCV { get; private set; }
        public bool METhroughSD { get; private set; }
        public bool AnyCensored { get; private set; }
        public bool AnyLT { get; private set; }
        public bool AnyGT { get; private set; }
        public bool AnyIntervalCensored { get; private set; }


        //DATA
        public double UncensoredSum { get { return raw.UncensoredSum; } }
        public double UncensoredSum2 { get { return raw.UncensoredSum2; } }
        public double[] Y { get { return raw.Y; } }
        public double[] LT { get { return raw.LT; } }
        public double[] GT { get { return raw.GT; } }
        public double[] IntervalGT { get { return raw.IntervalCensored.Item1; } }
        public double[] IntervalLT { get { return raw.IntervalCensored.Item2; } }

        public int YLength { get { return raw.Y.Length; } }
        public int LTLength { get { return raw.LT.Length; } }
        public int GTLength { get { return raw.GT.Length; } }
        public int IntervalLength { get { return raw.IntervalCensored.Item1.Length; } }

        public double[] LogY { get { return logOfRaw.Y; } }
        public double LogUncensoredSum { get { return logOfRaw.UncensoredSum; } }
        public double LogUncensoredSum2 { get { return logOfRaw.UncensoredSum2; } }
        public double[] LogLT { get { return logOfRaw.LT; } }
        public double[] LogGT { get { return logOfRaw.GT; } }
        public double[] LogIntervalGT { get { return logOfRaw.IntervalCensored.Item1; } }
        public double[] LogIntervalLT { get { return logOfRaw.IntervalCensored.Item2; } }

        public MeasureList MeasureList { get; private set; }
        public SortedList<string, Worker> WorkersByTag { get; private set; }
        public Worker[] WorkersByOrdinal { get; private set; }
        public int[] WorkerCounts { get; private set; }
        public int NWorkers { get { return this.WorkersByOrdinal.Length; } }

        public Measure[] UncensoredMeasures { get; private set; }
        public Measure[] GTMeasures { get; private set; }
        public Measure[] LTMeasures { get; private set; }
        public Measure[] IntervalMeasures { get; private set; }

        internal Dictionary<Worker, List<Measure>> MeasuresByWorker { get { return MeasureList.MeasuresByWorker; } }

        internal DataSummary(MeasureList ml, bool logNormalDistrn)
        {
            this.MeasureList = ml;
            this.UncensoredMeasures = MeasureList.measuresByCensoringType[Measure.MeasureType.Uncensored].ToArray();
            this.LTMeasures = MeasureList.measuresByCensoringType[Measure.MeasureType.LessThan].ToArray();
            this.GTMeasures = MeasureList.measuresByCensoringType[Measure.MeasureType.GreaterThan].ToArray();
            this.IntervalMeasures = MeasureList.measuresByCensoringType[Measure.MeasureType.Interval].ToArray();

            this.logNormalDistrn = logNormalDistrn;
            this.N = ml.N;
            this.MEAny = ml.MEAny;
            this.METhroughCV = ml.METhroughCV;
            this.METhroughSD = ml.METhroughSD;

            this.raw = new RawData();
            this.raw.Y = ml.measuresByCensoringType[Measure.MeasureType.Uncensored].Select(x => x.A).ToArray();
            this.raw.UncensoredSum = this.raw.Y.Sum();
            this.raw.UncensoredSum2 = this.raw.Y.Sqr().Sum();
            this.raw.LT = ml.measuresByCensoringType[Measure.MeasureType.LessThan].Select(x => x.A).ToArray();
            this.raw.GT = ml.measuresByCensoringType[Measure.MeasureType.GreaterThan].Select(x => x.A).ToArray();
            this.raw.IntervalCensored = Tuple.Create(
                    ml.measuresByCensoringType[Measure.MeasureType.Interval].Select(x => x.A).ToArray(),
                    ml.measuresByCensoringType[Measure.MeasureType.Interval].Select(x => x.B).ToArray());
            this.AnyGT = GTLength > 0;
            this.AnyLT = LTLength > 0;
            this.AnyIntervalCensored = IntervalLength > 0;
            this.AnyCensored = this.AnyGT || this.AnyLT || this.AnyIntervalCensored;

            if (ml.ME.ThroughCV || logNormalDistrn)
            {
                Func<double, bool> test;
                if (logNormalDistrn)
                {
                    test = delegate (double x)
                    {
                        return x <= 0;
                    };
                }
                else
                {
                    test = delegate (double x)
                    {
                        return x < 0;
                    };
                }

                string firstErrSrc = "";
                while (true)
                {
                    if (this.raw.Y.Any(xx => test(xx)))
                    {
                        firstErrSrc = "uncensored measure";
                        break;
                    }

                    if (this.raw.GT.Any(xx => test(xx)))
                    {
                        firstErrSrc = "'greater than' measure";
                        break;
                    }

                    if (this.raw.LT.Any(xx => test(xx)))
                    {
                        firstErrSrc = "'less than' measure";
                        break;
                    }

                    if (this.raw.IntervalCensored.Item1.Any(xx => test(xx)))
                    {
                        firstErrSrc = "'interval' measure";
                        break;
                    }

                    break;
                }

                if (firstErrSrc != string.Empty)
                {
                    //TODO:  Il faudra vérifier si cette exception peut ou non se produire.
                    string fmt = "{1} and negative data points (as found in {0}) are irreconciliable.";
                    firstErrSrc = string.Format(
                                    fmt,
                                    firstErrSrc,
                                    logNormalDistrn ? "A log-normally distributed outcome" : "Measurement Error specified through a Coefficient of Variation");
                    
                    throw new WEException(firstErrSrc);
                }
            }

            if (logNormalDistrn)
            {
                this.logOfRaw = new RawData();
                this.logOfRaw.Y = this.raw.Y.Log();
                this.logOfRaw.UncensoredSum = this.logOfRaw.Y.Sum();
                this.logOfRaw.UncensoredSum2 = this.logOfRaw.Y.Sqr().Sum();
                this.logOfRaw.GT = this.raw.GT.Log();
                this.logOfRaw.LT = this.raw.LT.Log();
                this.logOfRaw.IntervalCensored = Tuple.Create(
                    this.raw.IntervalCensored.Item1.Log(),
                    this.raw.IntervalCensored.Item2.Log());
            }

            this.WorkersByTag = MeasureList.WorkersByTag;
            for (int iw = 0; iw < this.WorkersByTag.Values.Count; iw++)
            {
                this.WorkersByTag.Values[iw].Ordinal = iw;
            }

            this.WorkerCounts = new int[this.MeasuresByWorker.Keys.Count];
            foreach(KeyValuePair<Worker, List<Measure>> kv in this.MeasuresByWorker)
            {
                this.WorkerCounts[kv.Key.Ordinal] = kv.Value.Count;
            }

            this.WorkersByOrdinal = this.WorkersByTag.Values.ToArray();
        }
    }
}
