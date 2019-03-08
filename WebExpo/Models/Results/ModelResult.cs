namespace Zygotine.WebExpo
{
    using System.Collections.Generic;
    using System.Linq;

    public class ModelResult
    {
        public McmcParameters McmcParameters { get; internal set; }
        public string MeasureList { get; internal set; }
        public string PastData = PastDataSummary.EmptyObject.ToString();
        public Messages Messages { get; internal set; } = new Messages();
        public bool Valid { get { return this.Messages.Level != Message.Severity.Error; } }
        public string ModelName { get; private set; }
        public string[] WorkerIds { get; private set; }
        internal Mcmc Chains { get; set; }
        internal Model ParentModel { get; set; }
        public uint PRNGSeed { get; internal set; }

        public ModelResult(Model model, params string[] mcmcPrefixes)
        {
            if((model == null) || (!model.Valid))
            {
                throw new WEException("In ModelResult constructor, parameter model must refer to a valid Model.");
            }

            this.ParentModel = model;
            this.ModelName = model.GetType().Name;
            this.Chains = new Mcmc(this, mcmcPrefixes);
            this.MeasureList = model.Data.MeasureList.ToString();
            this.WorkerIds = new string[model.Data.MeasureList.WorkerTags.Length];
            model.Data.MeasureList.WorkerTags.CopyTo(this.WorkerIds, 0);
            this.PastData = null;
            this.McmcParameters = new McmcParameters(nIter: ParentModel.NIter, nBurnin: ParentModel.NBurnin, nThin: ParentModel.NThin, monitorBurnin: ParentModel.MonitorBurnin);
        }

        public string[] GetChainNames()
        {
            return Chains.ChainDictionary.Keys.ToArray();
        }

        public double[] GetChainByName(string chainName)
        {
            double[] rep = Chains.GetChainByName(chainName);
            return rep;
        }

        public void ToCsvFile(string fileName)
        {
            Chains.ToCsvFile(fileName);
        }

        public Dictionary<string, double[]> GetQuantiles(double[] probs)
        {

            return this.Chains.GetQuantiles(probs);
        }

        internal void AddWorkerMuChains(IList<Worker> workers)
        {
            Chains.AddWorkerMuChain(workers);
        }

    }
}
