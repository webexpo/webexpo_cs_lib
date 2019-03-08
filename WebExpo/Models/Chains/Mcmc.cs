namespace Zygotine.WebExpo
{
    using System.Collections.Generic;
    using System.Linq;

    internal class Mcmc
    {
        private int burninSize, sampleSize;
        private ChainPairCreator cpc;
        internal Dictionary<string, Chain> ChainDictionary { get; private set; } = null;
        private Statistics.Quantile quantileDef;
        public string QuantileHeader { get; private set; }

        private Mcmc(ModelResult modelResult)
        {
            this.burninSize = modelResult.ParentModel.MonitorBurnin ? modelResult.ParentModel.NBurnin : 0;
            this.sampleSize = modelResult.ParentModel.NIter;
            this.cpc = new ChainPairCreator(burninSize, sampleSize);
            this.ChainDictionary = new Dictionary<string, Chain>();
            this.quantileDef = new Statistics.Quantile();
            this.QuantileHeader = "," + string.Join(",", quantileDef.Probabilities.Select(p => string.Format("{0:P1}", p)));
            //chainNo = chainNumber;
        }

        /**
         * Le paramètre modelResult ne doit pas être null, il permet d'obtenir les valeurs de MonitorBurnin, NBurnin, NIter.
         */
        internal Mcmc(ModelResult modelResult, params string[] chainNamePrefixes) : this(modelResult)
        {
            foreach (string prefix in chainNamePrefixes)
            {
                ChainPairCreator.ChainPair pair = this.cpc.GetPair(prefix);
                ChainDictionary.Add(pair.Burnin.Label, pair.Burnin);
                ChainDictionary.Add(pair.Sample.Label, pair.Sample);
            }
        }

        internal Dictionary<string, double[]> GetQuantiles(double[] probs)
        {
            Zygotine.Statistics.Quantile q = new Zygotine.Statistics.Quantile(probs);
            Dictionary<string, double[]> rep = new Dictionary<string, double[]>();


            foreach (KeyValuePair<string, Chain> kv in ChainDictionary)
            {
                if (kv.Value.Length > 0)
                {
                    rep.Add(kv.Key,  q.Compute(kv.Value.Values));
                } else
                {
                    rep.Add(kv.Key, new double[0]);
                }
            }
            return rep;
        }

        internal double[] GetChainByName(string chainName)
        {
            return
                    (ChainDictionary.ContainsKey(chainName))
                    ? ChainDictionary[chainName].Values
                    : null;
        }

        internal void Add(string chainNamePrefix)
        {
            ChainPairCreator.ChainPair pair = cpc.GetPair(chainNamePrefix);
            ChainDictionary.Add(pair.Burnin.Label, pair.Burnin);
            ChainDictionary.Add(pair.Sample.Label, pair.Sample);
        }

        internal void AddWorkerMuChain(IList<Worker> worker)
        {
            ChainPairCreator cpc = new ChainPairCreator(burninSize, sampleSize);
            foreach (Worker w in worker)
            {
                ChainPairCreator.ChainPair cp = cpc.GetPair(GetWorkerChainPrefix(w.Tag));
                ChainDictionary.Add(cp.Burnin.Label, cp.Burnin);
                ChainDictionary.Add(cp.Sample.Label, cp.Sample);
            }
        }

        internal double[] GetChain(string chainName)
        {
            if (ChainDictionary.ContainsKey(chainName))
                return ChainDictionary[chainName].Values;
            else
                return null;
        }

        public virtual void ToCsvFile(string filePath)
        {
            string[] oneLine = new string[ChainDictionary.Count / 2];
            Chain[] sampleChains = new Chain[ChainDictionary.Count / 2];
            int columnIndex = 0;

            foreach (Chain chain in ChainDictionary.Values)
            {
                if (chain.IsSample)
                {
                    sampleChains[columnIndex++] = chain;
                }
            }

            columnIndex = 0;
            foreach (Chain chain in sampleChains)
            {
                oneLine[columnIndex++] = string.Format(@"""{0}""", chain.Header);
            }

            using (System.IO.StreamWriter sw = new System.IO.StreamWriter(filePath))
            {
                sw.WriteLine(string.Join(",", oneLine)); //headers
                for (int j = 0; j < sampleSize; j++)
                {
                    columnIndex = 0;
                    foreach (Chain chain in sampleChains)
                    {
                        oneLine[columnIndex++] = string.Format("{0:G16}", chain.Values[j]);
                    }

                    sw.WriteLine(string.Join(",", oneLine));
                }
            }
        }

        public static ChainNamePair GetChainNames(string prefix)
        {
            return ChainPairCreator.GetChainNames(prefix);
        }

        public static string GetWorkerChainPrefix(string workerTag)
        {
            return "mu_" + workerTag;
        }

        public static ChainNamePair GetWorkerChainNames(string workerTag)
        {
            string prefix = "mu_" + workerTag;
            return ChainPairCreator.GetChainNames(prefix);
        }

    } // Mcmc


}
