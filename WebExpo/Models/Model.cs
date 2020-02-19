namespace Zygotine.WebExpo
{
    using System;
    using System.Collections.Generic;
    using System.Text.RegularExpressions;
    using System.Linq;

    public abstract class Model
    {
        public ModelResult Result { get; internal set; } = null;
        public Object[] ResultParams { get; set; } = null;
        internal DataSummary Data { get; set; }
        internal int NIter { get; set; }
        internal int NBurnin { get; set; }
        internal int NThin { get; set; }
        internal bool MonitorBurnin { get; set; }
        internal double OEL { get; set; }
        public bool OutcomeIsLogNormallyDistributed { get; private set; } = true;
        public static Version Version { get; private set; } = new Version(0, 3);
        internal MeasurementError ME { get; private set; }
        public MeasureList Measures { get; private set; }
        public Messages Messages { get; private set; } = new Messages();
        public uint PRNGSeed { get; private set; }
        protected string ClassName;

        internal class OelStandardizer
        {
            Model M;

            internal OelStandardizer(Model m)
            {
                M = m;
            }

            internal void StandardizeObservations()
            {
                if (!M.Measures.ObsStandardized)
                {
                    foreach (Measure m in M.Measures.measuresList)
                    {
                        if (!double.IsNaN(m.A))
                        {
                            m.A = m.A / M.Measures.OEL;
                        }
                        if (!double.IsNaN(m.B))
                        {
                            m.B = m.B / M.Measures.OEL;
                        }
                    }

                    M.Measures.ObsStandardized = true;
                }
            }

            internal void StandarizePastData()
            {
                ((SEGInformedVarModel)M).PastData.Mean -= Math.Log(M.OEL);
            }

            internal void AdjustMuChains()
            {
                Regex muChainRegex = new Regex(@"^mu.*Sample$");
                string[] chainIdsStandarize = M.Result.GetChainNames().Where(cn => muChainRegex.Match(cn).Success).ToArray<string>();
                foreach (string chainId in chainIdsStandarize)
                {
                    double[] c = M.Result.GetChainByName(chainId);
                    bool workerChain = chainId != "muSample" && chainId != "muOverallSample";
                    double[] muOverallChain = workerChain ? M.Result.GetChainByName("muOverallSample") : null;
                    for (int i = 0; i < c.Length; i++)
                    {
                        double delta = workerChain? muOverallChain[i] : Math.Log(M.OEL);
                        c[i] += delta;
                    }
                }
            }
        }

        internal OelStandardizer OelStdz;

        public Model(MeasureList measures, bool outcomeIsLogNormallyDistributed, McmcParameters mcmcParams)
        {
            OelStdz = new OelStandardizer(this);
            Measures = measures;

            this.ClassName = this.GetType().Name;
            if (mcmcParams == null)
            {
                mcmcParams = new McmcParameters(); // default values
                this.Messages.AddWarning("Mcmc parameters were undefined, using default values", this.ClassName);
            }

            this.Messages.Add(measures.Error);
            this.Messages.Add(mcmcParams.Error);
            this.NIter = mcmcParams.NIter;
            this.NBurnin = mcmcParams.NBurnin;
            this.NThin = mcmcParams.NThin;
            this.MonitorBurnin = mcmcParams.MonitorBurnin;
            this.OutcomeIsLogNormallyDistributed = outcomeIsLogNormallyDistributed;
            if (this.OutcomeIsLogNormallyDistributed)
            {
                OelStdz.StandardizeObservations();
            }
            this.Data = new DataSummary(measures, this.OutcomeIsLogNormallyDistributed);
            this.ME = measures.ME;
            this.OEL = measures.OEL;
        }
        protected static Model.MESupport MeasurementErrorSupport = new MESupport(me_CV: false, meSD: false);

        public Dictionary<string, double[]> GetQuantiles(double[] probs)
        {
            return this.Result.Chains.GetQuantiles(probs);
        }

        public class MESupport
        {
            public bool MEAny { get { return me_CV || MESD; } }
            public bool me_CV { get; private set; } = false;
            public bool MESD { get; private set; } = false;

            public MESupport(bool me_CV, bool meSD)
            {
                this.me_CV = me_CV;
                this.MESD = meSD;
            }
        }

        public bool Valid { get {
                return this.Messages.Level != Message.Severity.Error;
            } }

        internal abstract void Run();

        public void Compute()
        {
            uint prngSeed = Zygotine.Statistics.RNG.GetNewSeed();
            this.Compute(prngSeed);
        }

        public void Compute(uint prngSeed)
        {
            if (!this.Valid)
            {
                this.Result.Messages.AddError("Tentative de faire rouler un modèle qui n'est pas valide.", this.ClassName);
                return;
            }

            Zygotine.Statistics.RNG.SetSeed(prngSeed);
            this.PRNGSeed = prngSeed;
            this.Result.PRNGSeed = prngSeed;
            if (this.OutcomeIsLogNormallyDistributed)
            {
                if (this.GetType() == typeof(SEGInformedVarModel))
                {
                    OelStdz.StandarizePastData();
                }
            }
            this.Run();
            if (this.OutcomeIsLogNormallyDistributed) {
                OelStdz.AdjustMuChains();
            }
        }
    }
}
