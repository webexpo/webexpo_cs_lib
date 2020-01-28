namespace Zygotine.WebExpo
{
    using System;
    using System.Collections.Generic;

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
        public Model(MeasureList measures, bool outcomeIsLogNormallyDistributed, McmcParameters mcmcParams)
        {
            Measures = measures;

            this.ClassName = this.GetType().Name;
            if (mcmcParams == null)
            {
                mcmcParams = new McmcParameters(); // default values
                this.Messages.AddWarning("Mcmc parameters were undefined, using default values", this.ClassName);
            }

            if ( outcomeIsLogNormallyDistributed )
            {
                measures.StandardizeObservations(measures.OEL);
            }

            this.Messages.Add(measures.Error);
            this.Messages.Add(mcmcParams.Error);
            this.NIter = mcmcParams.NIter;
            this.NBurnin = mcmcParams.NBurnin;
            this.NThin = mcmcParams.NThin;
            this.MonitorBurnin = mcmcParams.MonitorBurnin;
            this.OutcomeIsLogNormallyDistributed = outcomeIsLogNormallyDistributed;
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
            if(!this.Valid)
            {
                this.Result.Messages.AddError("Tentative de faire rouler un modèle qui n'est pas valide.", this.ClassName);
                return;
            }

            Zygotine.Statistics.RNG.SetSeed(prngSeed);
            this.PRNGSeed = prngSeed;
            this.Result.PRNGSeed = prngSeed;
            this.Run();
            if ( this.OutcomeIsLogNormallyDistributed )
            {
                this.Result.StandardizeChains(this.OEL);
            }
        }

    }
}
