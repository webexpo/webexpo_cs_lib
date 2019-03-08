namespace Zygotine.WebExpo
{
    using System;
    using System.Collections.Generic;
    using System.Linq;
    using Zygotine.Statistics.Distribution;
    using Zygotine.Util;

    public class BetweenWorkerModel : Model
    {
        public double LogSigmaBetweenMu { get; private set; }
        public double LogSigmaBetweenPrec { get; private set; }
        public double LogSigmaWithinMu { get; private set; }
        public double LogSigmaWithinPrec { get; private set; }
        public double MuOverallLower { get; private set; }
        public double MuOverallUpper { get; private set; }
        public double[] MuOverallRange { get; private set; }
        public bool UseUniformPriorOnSds { get; private set; }
        public double[] SigmaBetweenRange { get; private set; }
        public double[] SigmaWithinRange { get; private set; }
        public double[] NormalizedY { get; private set; }
        public double InitMuOverall = Tools.ND;
        public double InitSigmaWithin = Tools.ND;

        static BetweenWorkerModel()
        {
            MeasurementErrorSupport = new MESupport(me_CV: false, meSD: false);
        }

        public BetweenWorkerModel(MeasureList measures, AbstractBWModelParameters specificParams, McmcParameters mcmcParams)
            : base(measures, specificParams.LogNormalDstrn, mcmcParams)
        {
            this.Result = new ModelResult(this, "muOverall", "sigmaWithin", "sigmaBetween");
            if (!measures.IsWorkerInfoComplete)
            {
                this.Messages.AddError("Incomplete worker information.", this.ClassName);
            }

            if (this.Data.MEAny)
            {
                if (this.Data.METhroughSD && !MeasurementErrorSupport.MESD)
                {
                    this.Messages.AddError("No support for METhrougSD.", this.ClassName);
                }
                else if (this.Data.METhroughCV && !MeasurementErrorSupport.me_CV)
                {
                    this.Messages.AddError("No support for METhrougCV.", this.ClassName);
                }
            }

            if (this.Messages.Level == Message.Severity.Error)
            {
                return;
            }

            this.UseUniformPriorOnSds = specificParams.UseUniformPriorOnSds;
            if (this.UseUniformPriorOnSds)
            {
                this.SigmaBetweenRange = Tools.Copy(((BWModelParameters_UniformPriorOnSDs)specificParams).SigmaBetweenRange);
                this.SigmaWithinRange = ((BWModelParameters_UniformPriorOnSDs)specificParams).SigmaWithinRange;
            }
            else
            {
                this.LogSigmaBetweenMu = ((BWModelParameters)specificParams).LogSigmaBetweenMu;
                this.LogSigmaBetweenPrec = ((BWModelParameters)specificParams).LogSigmaBetweenPrec;
                this.LogSigmaWithinMu = ((BWModelParameters)specificParams).LogSigmaWithinMu;
                this.LogSigmaWithinPrec = ((BWModelParameters)specificParams).LogSigmaWithinPrec;
            }

            //4 paramètres communs aux 2 formes que prend le modèle.
            this.MuOverallLower = specificParams.MuOverallLower;
            this.MuOverallUpper = specificParams.MuOverallUpper;
            this.InitMuOverall = specificParams.InitMuOverall;
            this.InitSigmaWithin = specificParams.InitSigmaWithin;

            this.MuOverallRange = Tools.Combine(this.MuOverallLower, this.MuOverallUpper);
            this.Result.AddWorkerMuChains(Data.WorkersByTag.Values);

            if (!Data.AnyCensored && UseUniformPriorOnSds)
            {
                //
                //# If there is no variation observed within subjects, 
                //# then we require a non-null lower bound for sigma.within
                //
                if (SigmaWithinRange[0] == 0)
                {
                    bool problem = (Data.NWorkers == Data.N); //il devrait donc y avoir plus d'une mesure pour au moins un travailleur
                    if (!problem)
                    {
                        problem = true;
                        foreach (KeyValuePair<Worker, List<Measure>> kv in Data.MeasuresByWorker)
                        {
                            if (kv.Value.Count <= 1)
                            {
                                continue;
                            }
                            //at least two measures
                            double vari = kv.Value.Select(x => x.A).Variance();
                            if (double.IsNaN(vari) || (vari == 0))
                            {
                                continue;
                            }
                            else
                            {
                                problem = false;
                                break;
                            }
                        } //foreach
                    }

                    if (problem)
                    {
                        this.Messages.AddError("Lower bound for SigmaWithinRange must be > 0 due to observed null within-subjects variance.", this.ClassName);
                        return;
                    }
                }
            }

            if (specificParams.GetType() == typeof(BWModelParameters))// && (Tools.IsND(this.InitMuOverall) || Tools.IsND(this.InitSigmaWithin)))
            {
                InitialValues initValues = specificParams.DefaultInits(this.Data, includeCensoredData: false);

                if (Tools.IsND(this.InitMuOverall))
                {
                    this.InitMuOverall = initValues.Mu;
                }

                if (Tools.IsND(this.InitSigmaWithin))
                {
                    this.InitSigmaWithin = initValues.SigmaWithin;
                }

                this.NormalizedY = initValues.NormalizedY.ToArray();
            }
        }

        internal override void Run()
        {
            ChainNamePair b_s;
            b_s = Mcmc.GetChainNames("muOverall");
            double[] muOverallBurnin = Result.Chains.GetChain(b_s.Burnin);
            double[] muOverallSample = Result.Chains.GetChain(b_s.Sample);
            b_s = Mcmc.GetChainNames("sigmaWithin");
            double[] sigmaWithinBurnin = Result.Chains.GetChain(b_s.Burnin);
            double[] sigmaWithinSample = Result.Chains.GetChain(b_s.Sample);
            b_s = Mcmc.GetChainNames("sigmaBetween");
            double[] sigmaBetweenBurnin = Result.Chains.GetChain(b_s.Burnin);
            double[] sigmaBetweenSample = Result.Chains.GetChain(b_s.Sample);

            double[][] workerBurnin = new double[Data.NWorkers][];
            double[][] workerSample = new double[Data.NWorkers][];

            int iTag = 0;
            foreach (string tag in Data.WorkersByTag.Keys)
            {
                b_s = Mcmc.GetWorkerChainNames(tag);
                workerBurnin[iTag] = Result.Chains.GetChain(b_s.Burnin);
                workerSample[iTag] = Result.Chains.GetChain(b_s.Sample);
                iTag++;
            }

            int iter = -1, savedIter = 0;
            try
            {
                double logSigmaWithinSD = 1 / Math.Sqrt(LogSigmaWithinPrec);
                double logSigmaBetweenSD = 1 / Math.Sqrt(LogSigmaBetweenPrec);

                //# Prepare dens.gen.icdf objects
                if (this.ME.Any)
                {
                    //o.tv < -truevalue.gen.object(me, outcome.is.logNormally.distributed)
                }

                if (this.ME.ThroughCV && !this.OutcomeIsLogNormallyDistributed)
                {
                    //o.mu.overall < -mu.truncatedData.gen.local.object(data$N, data$worker$count)
                    //o.mu.worker < -mu.worker.truncatedData.gen.object(data$worker$count)
                }

                GenObject oSB = null, oSW = null;
                if (this.UseUniformPriorOnSds)
                {
                    if (Data.NWorkers <= 1)
                    {
                        oSB = new Sigma_woLM_GenObject(Data.NWorkers);
                    }
                    else
                    {
                        oSB = null;
                    }
                }
                else
                {
                    oSB = new Sigma_LM_GenObject(Data.NWorkers, this.LogSigmaBetweenMu, logSigmaBetweenSD);
                }

                if (this.ME.ThroughCV && !this.OutcomeIsLogNormallyDistributed)
                {
                    //if (use.uniform.prior.on.sds)
                    //{
                    //    o.sw < -sigma.within.truncatedData.gen.object(data$N, data$worker$count, T, range = sigma.within.range)
                    //}
                    //else
                    //{
                    //    o.sw < -sigma.within.truncatedData.gen.object(data$N, data$worker$count, F, lnorm.mu = log.sigma.within.mu, lnorm.sd = log.sigma.within.sd)
                    //}
                }
                else
                {
                    if (this.UseUniformPriorOnSds)
                    {
                        if (Data.N <= 1)
                        {
                            oSW = new Sigma_woLM_GenObject(Data.N);
                        }
                        else
                        {
                            oSW = null;
                        }
                    }
                    else
                    {
                        oSW = new Sigma_LM_GenObject(this.Data.N, lNormMu: this.LogSigmaWithinMu, lNormSd: logSigmaWithinSD);
                    }
                }

                if (this.ME.Any && !this.ME.Known)
                {
                    //o.me < -me.gen.object(me, outcome.is.logNormally.distributed, data$N)
                }

                double muOverall = this.InitMuOverall;
                double sigmaWithin = InitSigmaWithin;

                //# Initialize measured values for subjects with censored values [new_0.10]
                YGen genY = YGen.GetEmptyObject();
                WorkerData workerData = new WorkerData(this);

                this.ResultParams = new Object[1];
                this.ResultParams[0] = this.Data.WorkersByTag.Keys;

                if (this.Data.AnyCensored)
                {
                    genY = YGen.Inits(data: this.Data, mu: muOverall, sigma: sigmaWithin, meThroughCV: false, logNormalDistrn: OutcomeIsLogNormallyDistributed);
                    workerData.UpdateGeneratedValues(genY);
                }


                double[] muWorker = workerData.MuWorker;
                workerData.AdjustMuOverall(this.MuOverallLower, this.MuOverallUpper);
                muOverall = workerData.MuOverall;
                muWorker = muWorker.Substract(muOverall); // # center mu.worker

                double[] predicted = workerData.GetPredictedMeans(muOverall);
                sigmaWithin = workerData.GetSigmaWithin();


                int nIterations = NBurnin + NIter * NThin;

                for (iter = 0; iter < nIterations; iter++)
                {
                    if (this.ME.Any)
                    {
                        //true.values < -truevalues.gen.local(gen.y, data, me, mu.worker, o = o.tv)
                    }

                    //# Sample y values for censored observations
                    if (this.Data.AnyCensored)
                    {
                        //function(true.values, data, me, mu.worker, mu=mu.overall, sigma=sigma.within, logNormal.distrn=outcome.is.logNormally.distributed)
                        genY = workerData.GenYLocal(muWorker, muOverall, sigmaWithin);
                        workerData.UpdateGeneratedValues(genY);
                    }

                    double[] yWorkerAvg = workerData.MuWorker;
                    double yAvg = workerData.Average;

                    //# Sample from f(sigma.within | other parms)
                    double[] residuals = workerData.GetGenValues().Substract(muWorker.Extract(workerData.WorkerIds)).Substract(muOverall);
                    double b = residuals.Sqr().ToArray().Sum() / 2.0;
                    SGNFnAParam localA = null;

                    if (this.ME.ThroughCV && !this.OutcomeIsLogNormallyDistributed)
                    {
                        //A < -c(o.sw$A, list(b = b, mu = mu.overall, muk = mu.worker))
                        //sigma.within < -dens.gen.icdf(o.sw, A, range = o.sw$range, start = sigma.within, inestimable.lower.limit = o.sw$inestimable.lower.limit)
                    }
                    else
                    {
                        if (this.UseUniformPriorOnSds)
                        {
                            sigmaWithin = WebExpoFunctions3.SqrtInvertedGammaGen(Data.N, b, SigmaWithinRange.Copy(), oSW);
                        }
                        else
                        {
                            localA = oSW.A.Clone();
                            localA.B = b;
                            localA.Mu = muOverall;
                            localA.Muk = muWorker;
                            Icdf icdf = new Icdf(oSW, localA, range: Tools.Combine(0, double.PositiveInfinity));
                            sigmaWithin = icdf.Bidon(start: oSW.Start(localA), inestLowerLim: true);
                        }
                    }

                    //# Sample from f(sigma.between | other parms)
                    b = muWorker.Sqr().Sum() / 2.0;

                    double sigmaBetween = 0;
                    if (UseUniformPriorOnSds)
                    {
                        sigmaBetween = WebExpoFunctions3.SqrtInvertedGammaGen(Data.NWorkers, b, this.SigmaBetweenRange.Copy(), oSB);
                    }
                    else
                    {
                        localA = oSB.A.Clone();
                        localA.B = b;
                        Icdf icdf = new Icdf(oSB, localA, range: Tools.Combine(0, double.PositiveInfinity));
                        sigmaBetween = icdf.Bidon(start: oSB.Start(localA), inestLowerLim: true);
                    }

                    //# Sample from f(mu.overall | other parms)
                    //# modif_0.10
                    double tmpMean = yAvg - (muWorker.Multiply(this.Data.WorkerCounts).Sum() / this.Data.N);


                    if (this.ME.ThroughCV && !OutcomeIsLogNormallyDistributed)
                    {
                        //muOverall = mu.truncatedData.gen.local(o.mu.overall, tmp.mean, sigma.within, mu.worker, mu.overall.range, current.value = mu.overall)
                    }
                    else
                    {
                        double tmpSD = sigmaWithin / Math.Sqrt(this.Data.N);
                        muOverall = RNorm4CensoredMeasures.RNormCensored(tmpMean, tmpSD, lower: this.MuOverallLower, upper: this.MuOverallUpper);
                    }

                    //# Sample from f(mu.worker's | other parms)
                    //# modif_0.10

                    double[] sigma2A = Tools.Rep(Math.Pow(sigmaWithin, 2.0), this.Data.NWorkers).Divide(this.Data.WorkerCounts);   // # vector of length '# of workers' 
                    double sigma2B = Math.Pow(sigmaBetween, 2);                                     // # scalar
                    double[] muA = yWorkerAvg.Substract(muOverall);                                 // # vector of length '# of workers'
                    double[] mukStar = muA.Multiply(sigma2B).Divide(sigma2A.Add(sigma2B));          // # vector of length '# of workers'
                    double[] s2kStar = sigma2A.Multiply(sigma2B).Divide(sigma2A.Add(sigma2B));      // # vector of length '# of workers'

                    if (this.ME.ThroughCV && !OutcomeIsLogNormallyDistributed)
                    {
                        //muWorker = mu.worker.truncatedData.gen(o.mu.worker, muk.star, s2k.star, mu.overall, sigma.within, mu.worker)
                    }
                    else
                    {
                        muWorker = NormalDistribution.RNorm(this.Data.NWorkers, mu: mukStar, sigma: s2kStar.Sqrt());
                    }

                    if (iter < NBurnin)
                    {
                        if (MonitorBurnin)
                        {
                            muOverallBurnin[iter] = muOverall;
                            sigmaBetweenBurnin[iter] = sigmaBetween;
                            sigmaWithinBurnin[iter] = sigmaWithin;
                            SaveWorkerChainData(iter, muWorker, workerBurnin);
                        }
                    }
                    else if ((iter - NBurnin) % NThin == 0)
                    {
                        muOverallSample[savedIter] = muOverall;
                        sigmaBetweenSample[savedIter] = sigmaBetween;
                        sigmaWithinSample[savedIter] = sigmaWithin;
                        SaveWorkerChainData(savedIter, muWorker, workerSample);
                        savedIter++;
                    }
                } //for ...
            }
            catch (Exception ex)
            {
                this.Result.Messages.AddError(WEException.GetStandardMessage(ex, iter, Result.PRNGSeed), this.ClassName);
                return;
            }
        }

        private void SaveWorkerChainData(int position, double[] dataForWorker, double[][] target)
        {
            for (int iworker = 0; iworker < dataForWorker.Length; iworker++)
            {
                target[iworker][position] = dataForWorker[iworker];
            }
        }

        private struct Position
        {
            public int Start;
            public int Length;

            public Position(int begin, int len)
            {
                this.Start = begin;
                this.Length = len;
            }
        }

        private class WorkerData
        {
            BetweenWorkerModel parentModel;
            private DataSummary Data;
            private WorkerDigest[] workerDigests;
            private DfRecord[] df;

            Dictionary<Measure.MeasureType, List<DfRecord>> DFRecordByMeasureType = new Dictionary<Measure.MeasureType, List<DfRecord>>();
            internal double MuOverall { get; set; } = 0;
            internal double[] MuWorker { get; private set; }
            internal int[] MeasureCountByWorker { get; private set; }
            internal double Average { get; private set; }
            private Position[] Positions { get; set; } = new Position[4];
            internal int[] WorkerIds { get; set; }
            private bool logNormalDistrn;

            internal WorkerData(BetweenWorkerModel parent)
            {
                this.parentModel = parent;
                this.logNormalDistrn = this.parentModel.OutcomeIsLogNormallyDistributed;
                Data = parent.Data;
                workerDigests = new WorkerDigest[Data.NWorkers];
                int i = 0;
                this.MeasureCountByWorker = new int[Data.NWorkers];

                foreach (Worker w in Data.WorkersByOrdinal)
                {
                    MeasureCountByWorker[i] = Data.MeasuresByWorker[w].Count;
                    workerDigests[i] = new WorkerDigest(ordinal: i, mean: 0, measureCount: Data.MeasuresByWorker[w].Count);
                    i++;
                }

                DFRecordByMeasureType.Add(Measure.MeasureType.Uncensored, new List<DfRecord>());
                DFRecordByMeasureType.Add(Measure.MeasureType.GreaterThan, new List<DfRecord>());
                DFRecordByMeasureType.Add(Measure.MeasureType.LessThan, new List<DfRecord>());
                DFRecordByMeasureType.Add(Measure.MeasureType.Interval, new List<DfRecord>());

                df = new DfRecord[Data.N];
                int j = 0;
                Positions[(int)Measure.MeasureType.Uncensored] = new Position(j, Data.UncensoredMeasures.Length);
                for (i = 0; i < Data.UncensoredMeasures.Length; i++)
                {
                    int workerOrdinal = Data.UncensoredMeasures[i].Worker.Ordinal;
                    Measure m = Data.UncensoredMeasures[i];
                    m.WorkerDataOrdinal = j;
                    df[j] = new DfRecord(
                        measureOrdinal: j,
                        genValue: this.logNormalDistrn ? Math.Log(m.A) : m.A,
                        workerOrdinal: workerOrdinal,
                        workerDigest: workerDigests[workerOrdinal]);
                    DFRecordByMeasureType[Measure.MeasureType.Uncensored].Add(df[j]);
                    j++;
                }

                Positions[(int)Measure.MeasureType.GreaterThan] = new Position(j, Data.GT.Length);
                for (i = 0; i < Data.GTMeasures.Length; i++)
                {
                    int workerOrdinal = Data.GTMeasures[i].Worker.Ordinal;
                    Measure m = Data.GTMeasures[i];
                    m.WorkerDataOrdinal = j;
                    df[j] = new DfRecord(
                        measureOrdinal: j,
                        genValue: Tools.ND,
                        workerOrdinal: workerOrdinal,
                        workerDigest: workerDigests[workerOrdinal]);
                    DFRecordByMeasureType[Measure.MeasureType.GreaterThan].Add(df[j]);
                    j++;
                }

                Positions[(int)Measure.MeasureType.LessThan] = new Position(j, Data.LT.Length);
                for (i = 0; i < Data.LTMeasures.Length; i++)
                {
                    int workerOrdinal = Data.LTMeasures[i].Worker.Ordinal;
                    Measure m = Data.LTMeasures[i];
                    m.WorkerDataOrdinal = j;
                    df[j] = new DfRecord(
                        measureOrdinal: j,
                        genValue: Tools.ND,
                        workerOrdinal: workerOrdinal,
                        workerDigest: workerDigests[workerOrdinal]);
                    DFRecordByMeasureType[Measure.MeasureType.LessThan].Add(df[j]);
                    j++;
                }

                Positions[(int)Measure.MeasureType.Interval] = new Position(j, Data.IntervalGT.Length);
                for (i = 0; i < Data.IntervalMeasures.Length; i++)
                {
                    int workerOrdinal = Data.IntervalMeasures[i].Worker.Ordinal;
                    Measure m = Data.IntervalMeasures[i];
                    m.WorkerDataOrdinal = j;
                    df[j] = new DfRecord(
                        measureOrdinal: j,
                        genValue: Tools.ND,
                        workerOrdinal: workerOrdinal,
                        workerDigest: workerDigests[workerOrdinal]);
                    DFRecordByMeasureType[Measure.MeasureType.Interval].Add(df[j]);
                    j++;
                }

                this.WorkerIds = df.Select(rec => rec.WorkerOrdinal).ToArray();
                UpdatePublicProperties();
            }

            internal void UpdateGeneratedValues(YGen genY)
            {
                int begin, len, p;

                p = (int)Measure.MeasureType.GreaterThan;
                begin = Positions[p].Start;
                len = Positions[p].Length;
                for (int i = 0; i < len; i++)
                {
                    df[begin + i].GenValue = this.logNormalDistrn ? genY.LogGT[i] : genY.GT[i];
                }

                p = (int)Measure.MeasureType.LessThan;
                begin = Positions[p].Start;
                len = Positions[p].Length;
                for (int i = 0; i < len; i++)
                {
                    df[begin + i].GenValue = this.logNormalDistrn ? genY.LogLT[i] : genY.LT[i];
                }

                p = (int)Measure.MeasureType.Interval;
                begin = Positions[p].Start;
                len = Positions[p].Length;
                for (int i = 0; i < len; i++)
                {
                    df[begin + i].GenValue = this.logNormalDistrn ? genY.LogI[i] : genY.I[i];
                }

                UpdatePublicProperties();
            }

            internal void AdjustMuOverall(double muOverallLower, double muOverallUpper)
            {
                if (this.MuOverall < muOverallLower)
                {
                    this.MuOverall = muOverallLower; // # new_0.10
                }
                else if (this.MuOverall > muOverallUpper)
                {
                    this.MuOverall = muOverallUpper; // # new_0.10
                }
            }

            internal double[] GetResiduals(double[] muWorkerUncentered)
            {

                double[] residuals = new double[df.Length];

                for (int i = 0; i < df.Length; i++)
                    //residuals < -df$y - mu.overall - mu.worker[df$id]
                    residuals[i] = df[i].GenValue - muWorkerUncentered[df[i].WD.Ordinal];
                return residuals;
            }

            internal YGen GenYLocal(double[] muWorker, double mu, double sigma)
            {
                YGen rep = YGen.GetEmptyObject();
                if (this.parentModel.ME.Any)
                {
                    //TODO non implanté
                    return null;
                }
                else
                {
                    if (this.logNormalDistrn)
                    {
                        if (Data.AnyGT)
                        {
                            int[] workerIds = DFRecordByMeasureType[Measure.MeasureType.GreaterThan].Select(rec => rec.WorkerOrdinal).ToArray();
                            double[] means = muWorker.Extract(workerIds).Add(mu);
                            int len = this.Data.LogGT.Length;
                            rep.GT = new double[len];
                            rep.LogGT = new double[len];
                            for (int i = 0; i < len; i++)
                            {
                                rep.LogGT[i] = RNorm4CensoredMeasures.RNormCensored(means[i], sigma, lower: this.Data.LogGT[i]);
                            }

                            rep.GT = rep.LogGT.Exp();
                        }

                        if (Data.AnyLT)
                        {
                            int[] workerIds = DFRecordByMeasureType[Measure.MeasureType.LessThan].Select(rec => rec.WorkerOrdinal).ToArray();
                            double[] means = muWorker.Extract(workerIds).Add(mu);
                            int len = this.Data.LogLT.Length;
                            rep.LT = new double[len];
                            rep.LogLT = new double[len];
                            for (int i = 0; i < len; i++)
                            {
                                rep.LogLT[i] = RNorm4CensoredMeasures.RNormCensored(means[i], sigma, upper: this.Data.LogLT[i]);
                            }

                            rep.LT = rep.LogLT.Exp();
                        }

                        if (Data.AnyIntervalCensored)
                        {
                            int[] workerIds = DFRecordByMeasureType[Measure.MeasureType.Interval].Select(rec => rec.WorkerOrdinal).ToArray();
                            double[] means = muWorker.Extract(workerIds).Add(mu);
                            int len = this.Data.LogIntervalGT.Length;
                            rep.I = new double[len];
                            rep.LogI = new double[len];
                            for (int i = 0; i < len; i++)
                            {
                                rep.LogI[i] = RNorm4CensoredMeasures.RNormCensored(means[i], sigma, lower: this.Data.LogIntervalGT[i], upper: this.Data.LogIntervalLT[i]);
                            }

                            rep.I = rep.LogI.Exp();
                        }
                    }
                    else
                    {
                        if (Data.AnyGT)
                        {
                            int[] workerIds = DFRecordByMeasureType[Measure.MeasureType.GreaterThan].Select(rec => rec.WorkerOrdinal).ToArray();
                            double[] means = muWorker.Extract(workerIds).Add(mu);
                            int len = this.Data.GT.Length;
                            rep.GT = new double[len];
                            for (int i = 0; i < len; i++)
                            {
                                rep.GT[i] = RNorm4CensoredMeasures.RNormCensored(means[i], sigma, lower: this.Data.GT[i]);
                            }
                        }

                        if (Data.AnyLT)
                        {
                            int[] workerIds = DFRecordByMeasureType[Measure.MeasureType.LessThan].Select(rec => rec.WorkerOrdinal).ToArray();
                            double[] means = muWorker.Extract(workerIds).Add(mu);
                            int len = this.Data.LT.Length;
                            rep.LT = new double[len];
                            for (int i = 0; i < len; i++)
                            {
                                rep.LT[i] = RNorm4CensoredMeasures.RNormCensored(means[i], sigma, upper: this.Data.LT[i]);
                            }
                        }

                        if (Data.AnyIntervalCensored)
                        {
                            int[] workerIds = DFRecordByMeasureType[Measure.MeasureType.Interval].Select(rec => rec.WorkerOrdinal).ToArray();
                            double[] means = muWorker.Extract(workerIds).Add(mu);
                            int len = this.Data.IntervalGT.Length;
                            rep.I = new double[len];
                            rep.LogI = new double[len];
                            for (int i = 0; i < len; i++)
                            {
                                rep.I[i] = RNorm4CensoredMeasures.RNormCensored(means[i], sigma, lower: this.Data.IntervalGT[i], upper: this.Data.IntervalLT[i]);
                            }
                        }
                    }
                }
                return rep;
            }//# end of y.gen.local

            private void UpdatePublicProperties()
            {   // This will set MuWorkerCentered, MuWorkerUncentered and MuOverall
                CalculateMeans();
            }

            private void CalculateMeans()
            {
                foreach (WorkerDigest wd in workerDigests)
                {
                    wd.Clear();
                }

                foreach (DfRecord dfRec in df)
                {
                    dfRec.WD.AddValue(dfRec);
                }

                this.Average = df.Select(rec => rec.GenValue).Average();

                MuOverall = 0;
                foreach (WorkerDigest wd in workerDigests)
                {
                    wd.CalculateMean();
                    MuOverall += wd.Mean;
                }

                MuOverall = MuOverall / workerDigests.Length;
                MuWorker = workerDigests.Select(dg => dg.Mean).ToArray();
            }

            internal double GetSigmaWithin()
            {
                return Math.Sqrt(df.Select(xx => xx.GenValue).ToArray().Substract(GetPredictedMeans(this.MuOverall)).Variance());
            }

            internal double[] GetPredictedMeans(double muOverall)
            {
                double[] muWorker = this.MuWorker.Substract(muOverall); // # center mu.worker
                double[] predicted = muWorker.Extract(this.WorkerIds).Add(muOverall);
                return predicted;
            }

            internal double[] GetGenValues()
            {
                return df.Select(rec => rec.GenValue).ToArray();
            }

            private class WorkerDigest
            {
                private double Accumulator = 0;
                public int Ordinal { get; private set; }
                public double Mean { get; private set; }
                public int MeasureCount { get; private set; }

                internal WorkerDigest(int ordinal, double mean, int measureCount)
                {
                    Ordinal = ordinal;
                    Mean = mean;
                    MeasureCount = measureCount;
                }

                internal void Clear()
                {
                    Mean = 0;
                    Accumulator = 0;
                }

                internal void AddValue(DfRecord rec)
                {
                    Accumulator += rec.GenValue;
                }

                internal void CalculateMean()
                {
                    Mean = Accumulator / MeasureCount;
                }
            }

            private class DfRecord
            {
                public int Ordinal { get; private set; }
                public double GenValue { get; internal set; }
                public double LogGenValue { get; internal set; }
                public int WorkerOrdinal { get; private set; }
                public WorkerDigest WD { get; private set; }

                internal DfRecord(int measureOrdinal, double genValue, int workerOrdinal, WorkerDigest workerDigest)
                {
                    Ordinal = measureOrdinal;
                    GenValue = genValue;
                    WorkerOrdinal = workerOrdinal;
                    WD = workerDigest;
                }
            }
        }
    }// end of class
}
