namespace Zygotine.WebExpo
{
    using System;
    using System.Linq;
    using Zygotine.Statistics.Distribution;
    using Zygotine.Util;

    public class SEGInformedVarModel : Model
    {

        //    LogDot3 = -1.20397280432594,
        //    Log2Dot5 = 0.916290731874155;

        public static readonly Version version = new Version(0, 3);
        public double LogSigmaMu { get; private set; }
        public double LogSigmaPrec { get; private set; }
        public double InitSigma { get; private set; }
        public double MuLower { get; private set; }
        public double MuUpper { get; private set; }
        public double InitMu { get; private set; }
        public PastDataSummary PastData { get; private set; }

        public SEGInformedVarModel(MeasureList measures, SEGInformedVarModelParameters specificParams, McmcParameters mcmcParams = null, PastDataSummary pastDataSummary = null)
            : base(measures, specificParams.LogNormalDstrn, mcmcParams)
        {
            this.PastData = pastDataSummary == null ? PastDataSummary.EmptyObject : pastDataSummary;

            this.Messages.Add(this.PastData.Messages);
            this.LogSigmaMu = specificParams.LogSigmaMu;
            this.LogSigmaPrec = specificParams.LogSigmaPrec;
            this.InitSigma = specificParams.InitSigma;
            this.MuLower = specificParams.MuLower;
            this.MuUpper = specificParams.MuUpper;
            this.InitMu = specificParams.InitMu;
            this.Result = new ModelResult(this, "mu", "sd");

            if (this.ME.ThroughCV)
            {
                this.Result.Chains.Add("cv");
            }
        }// end constructor

        /*
         * La méthode étant internal, elle peut être invoquée d'un programme externe à la librairie.
         * La seule méthode qui peut invoquer Run, c'est la méthode Compute de Model qui ne le fera que si le modèle est
         * jugé valide.
         */
        internal override void Run()
        {
            SGNFnAParam localA = null;
            GenObject oTV = null;
            GenObject oMu = null;
            GenObject oSigma = null;
            GenObject oME = null;
            YGen genY = YGen.EmptyInstance;
            TrueValuesGen genTV = null;
            double[] burninMu;
            double[] burninSigma;
            double[] burninCV = null;
            double[] sampleMu;
            double[] sampleSigma;
            double[] sampleCV = null;
            double mu;
            double sigma;
            int iter = -1, savedIter;
            double muCondMean;
            double yBar;
            double muCondSD;
            double[] pLim = new double[2];
            double p;
            double[] muLim = new double[] { this.MuLower, this.MuUpper };
            double logSigmaSD;

            try
            {
                logSigmaSD = 1 / Math.Sqrt(this.LogSigmaPrec);
                if (ME.Any)
                {
                    if (ME.ThroughCV)
                    {
                        if (OutcomeIsLogNormallyDistributed)
                        {
                            oTV = new TrueValue_CV_LogN_GenObject();
                        }
                        else
                        {
                            oTV = new TrueValue_CV_Norm_GenObject();
                        }
                    }
                    else
                    {
                        //oTV = new TrueValue_SD_GenObject();
                    }
                }

                //# modif_0.12
                int combinedN = this.Data.N + (this.PastData.Defined ? PastData.N : 0);
                if (ME.ThroughCV && !OutcomeIsLogNormallyDistributed)
                {
                    oMu = new MuTruncatedData_GenObject(combinedN); //# modif_0.12 
                    oSigma = GenObject.GetSigmaTruncatedDataLNormGenObject(combinedN, this.LogSigmaMu, logSigmaSD); //# modif_0.12 
                }
                else
                {
                    oSigma = GenObject.GetSigmaGenObject(combinedN, this.LogSigmaMu, logSigmaSD); //# modif_0.12 
                }

                localA = oSigma.A.Clone();
                if (ME.Any && !ME.Known)
                {
                    oME = GenObject.GetMeGenObject(this.ME, this.OutcomeIsLogNormallyDistributed, this.Data.N);
                }

                int nIterations = NBurnin + NIter * NThin;
                //les tableaux pour les chaines
                sampleMu = Result.Chains.GetChain("muSample");
                sampleSigma = Result.Chains.GetChain("sdSample");
                burninMu = Result.Chains.GetChain("muBurnin");
                burninSigma = Result.Chains.GetChain("sdBurnin");
                if (ME.ThroughCV)
                {
                    sampleCV = Result.Chains.GetChain("cvSample");
                    burninCV = Result.Chains.GetChain("cvBurnin");
                }

                bool inestimableLowerLimit = false;

                //Initial values for mu and sigma
                mu = InitMu;
                sigma = InitSigma;
                savedIter = 0; // pour les échantillons
                if (this.Data.AnyCensored)
                {
                    genY = YGen.Inits(this.Data, mu, sigma, meThroughCV: this.ME.ThroughCV, logNormalDistrn: OutcomeIsLogNormallyDistributed);
                }

                if (ME.Any)
                {
                    ME.Parm = ME.InitialValue;
                }


                //Boucle principale
                for (iter = 0; iter < nIterations; iter++)
                {
                    if (ME.Any)
                    {
                        genTV = TrueValuesGen.GetInstance(genY, this.Data, mu, sigma, this.ME, logNormalDistrn: OutcomeIsLogNormallyDistributed, o: oTV);
                    }

                    if (this.Data.AnyCensored)
                    {
                        //y.gen(true.values, data, sigma, me, outcome.is.logNormally.distributed, mu=mu)
                        //On ne tient pas compte de true.values, ni de me ...
                        genY = YGen.GetInstance(this.ME, genTV, this.Data, mu, sigma, OutcomeIsLogNormallyDistributed);
                    }

                    OutLogoutMoments moments = OutLogoutMoments.Get(this.ME.Any, this.OutcomeIsLogNormallyDistributed, this.Data, genY, genTV);
                    double sigmaBeta = (moments.Sum2 - 2 * mu * moments.Sum + this.Data.N * mu * mu) / 2.0;
                    if (PastData.Defined)
                    {
                        sigmaBeta = sigmaBeta + PastData.N / 2.0 * Math.Pow(PastData.Mean - mu, 2) + PastData.NS2 / 2.0;
                    }

                    double[] start = new double[0];
                    if (this.ME.ThroughCV && !OutcomeIsLogNormallyDistributed)
                    {
                        //ici
                        //        A <- c(o.sigma$A, list(b=sigma.beta, mu=mu))
                        localA = oSigma.A.Clone();
                        localA.B = sigmaBeta;
                        localA.Mu = mu;
                        start = Tools.Combine(sigma);
                        inestimableLowerLimit = false;
                    }
                    else
                    {
                        localA.B = sigmaBeta;
                        start = oSigma.Start(localA);
                        inestimableLowerLimit = true;
                    }

                    Icdf icdf = new Icdf(oSigma, localA, Tools.Combine(0, double.PositiveInfinity));
                    sigma = icdf.Bidon(start, inestimableLowerLimit);
                    yBar = moments.Sum / this.Data.N;
                    muCondMean = this.PastData.Defined ? (moments.Sum + PastData.N * PastData.Mean) / combinedN : yBar; // # new_0.12
                    if (this.ME.ThroughCV && !this.OutcomeIsLogNormallyDistributed)
                    {
                        mu = MuTruncatedGen.GetInstance(oMu, muLim, muCondMean, sigma).Mu;
                    }
                    else
                    {
                        muCondSD = sigma / Math.Sqrt(combinedN);
                        pLim = NormalDistribution.PNorm(muLim.Substract(muCondMean).Divide(muCondSD));
                        p = UniformDistribution.RUnif(1, pLim[0], pLim[1])[0];
                        mu = NormalDistribution.QNorm(p, mu: muCondMean, sigma: muCondSD);
                    }

                    //# Sample Measurement Error from its posterior density
                    if (this.ME.Any && !this.ME.Known)
                    {
                        this.ME.Parm = MEParmGen.GetInstance(oME, this.ME, this.Data, genY, genTV).Parm;
                    }

                    if (iter < NBurnin)
                    {
                        if (MonitorBurnin)
                        {
                            burninMu[iter] = mu;
                            burninSigma[iter] = sigma;
                            if (this.ME.Any && !this.ME.Known)
                            {
                                burninCV[iter] = ME.Parm;
                            }
                        }
                    }
                    else if ((iter - NBurnin) % NThin == 0)
                    {
                        sampleMu[savedIter] = mu;
                        sampleSigma[savedIter] = sigma;
                        if (this.ME.Any && !this.ME.Known)
                        {
                            sampleCV[savedIter] = ME.Parm;
                        }

                        savedIter++;
                    }

                }// for( int iter = 1 ...

            }
            catch (Exception ex)
            {
                this.Result.Messages.AddError(WEException.GetStandardMessage(ex, iter, Result.PRNGSeed), this.ClassName);
                return;
            }

        } //end Run
    }
}










