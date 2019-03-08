namespace Zygotine.WebExpo
{
    using System;
    using Zygotine.Util;

    public class SEGUninformativeModel : Model
    {
        //    LogDot3 = -1.20397280432594,
        //    Log2Dot5 = 0.916290731874155;
        public static readonly Version version = new Version(0, 3);
        public double IntervalLower { get; private set; }
        public double IntervalUpper { get; private set; }
        public double MuLower { get; private set; }
        public double MuUpper { get; private set; }
        public double[] SDRange { get; private set; }
        public double InitMu { get; private set; }
        public double InitSD { get; private set; }


        public SEGUninformativeModel(MeasureList measures, UninformativeModelParameters specificParams, McmcParameters mcmcParams = null)
            : base(measures, specificParams.LogNormalDstrn, /*new Range(specificParams.MuLower, specificParams.MuUpper, specificParams.InitMu),*/ mcmcParams)
        {
            this.MuLower = specificParams.MuLower;
            this.MuUpper = specificParams.MuUpper;
            this.SDRange = specificParams.SDRange;
            this.InitMu = specificParams.InitMu;
            this.InitSD = specificParams.InitSD;
            this.Result = new ModelResult(this, "mu", "sd");
            bool testSD = double.IsNaN(this.InitSD) || (this.InitSD == 0);

            if (Tools.IsND(this.InitMu) || testSD)
            {
                InitialValues initVals = WebExpoFunctions3.DefaultInits(this.Data, this.OutcomeIsLogNormallyDistributed, Tools.Combine(this.MuLower, this.MuUpper), this.SDRange, includeCensoredData: false);
                if (Tools.IsND(this.InitMu))
                {
                    this.InitMu = initVals.Mu;
                }

                if (testSD)
                {
                    this.InitSD = initVals.SigmaWithin;
                }
            }

            if (this.ME.ThroughCV)
            {
                this.Result.Chains.Add("cv");
            }
        }// end constructor

        internal override void Run()
        {
            YGen genY = YGen.EmptyInstance; // Pourra probablement prendre plusieurs formes.
            TrueValuesGen genTV = null;
            GenObject oSigma = EmptyGenObject.Instance;
            GenObject oMu = EmptyGenObject.Instance;
            GenObject oME = null;
            GenObject oTV = null;
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
            double muCondSD;

            try
            {
                //# Prepare dens.gen.icdf objects
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

                if (ME.ThroughCV && !OutcomeIsLogNormallyDistributed)
                {
                    oMu = new MuTruncatedData_GenObject(this.Data.N); //# modif_0.12 
                    oSigma = GenObject.GetSigmaTruncatedDataGenObject(this.Data.N); //# modif_0.12 
                }
                else
                {
                    oSigma = EmptyGenObject.Instance;
                }

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

                //Initial values for mu and sigma
                mu = this.InitMu;
                sigma = this.InitSD;

                //# Initialize measured values for subjects with censored values [modif_0.10]
                if (this.Data.AnyCensored)
                {
                    genY = YGen.Inits(data: this.Data, mu: mu, sigma: sigma, meThroughCV: true, logNormalDistrn: OutcomeIsLogNormallyDistributed);
                }

                if (ME.Any)
                {
                    ME.Parm = ME.InitialValue;
                }

                //# Start MCMC
                savedIter = 0; // pour les échantillons
                for (iter = 0; iter < nIterations; iter++)
                {
                    //# Sample true values (in presence of measurement error) [new_0.10]
                    if (ME.Any)
                    {
                        genTV = TrueValuesGen.GetInstance(genY, this.Data, mu, sigma, this.ME, OutcomeIsLogNormallyDistributed, o: oTV);
                    }

                    //# Sample y values for censored observations
                    if (this.Data.AnyCensored)
                    {
                        //y.gen(true.values, data, sigma, me, outcome.is.logNormally.distributed, mu=mu)
                        //On ne tient pas compte de true.values, ni de me ...
                        genY = YGen.GetInstance(this.ME, genTV, this.Data, mu, sigma, OutcomeIsLogNormallyDistributed);
                    }

                    //# Compute data points sum and square sum
                    OutLogoutMoments moments = OutLogoutMoments.Get(this.ME.Any, this.OutcomeIsLogNormallyDistributed, this.Data, genY, genTV);

                    //# Sample from f(sigma | mu)
                    //# modif_0.10
                    double sigmaBeta = (moments.Sum2 - 2 * mu * moments.Sum + this.Data.N * mu * mu) / 2.0;
                    if (sigmaBeta < 1e-6)
                    {
                        sigmaBeta = 0; // # protection against numeric imprecision
                    }

                    if (this.ME.ThroughCV && !OutcomeIsLogNormallyDistributed)
                    {
                        sigma = SigmaTruncatedDataGen.GetInstance(oSigma, SDRange, sigmaBeta, mu, sigma).Sigma;
                    }
                    else
                    {
                        sigma = WebExpoFunctions3.SqrtInvertedGammaGen(Data.N, sigmaBeta, this.SDRange, oSigma);
                    }

                    muCondMean = moments.Sum / this.Data.N;

                    if (this.ME.ThroughCV && !this.OutcomeIsLogNormallyDistributed)
                    {
                        mu = MuTruncatedGen.GetInstance(oMu, Tools.Combine(MuLower, MuUpper), muCondMean, sigma).Mu;
                    }
                    else
                    {
                        muCondSD = sigma / Math.Sqrt(this.Data.N);
                        mu = RNorm4CensoredMeasures.RNormCensored(muCondMean, muCondSD, lower: this.MuLower, upper: this.MuUpper);
                    }

                    //# Sample Measurement Error from its posterior density
                    if (this.ME.Any && !ME.Known)
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
            } catch(Exception ex)
            {
                this.Result.Messages.AddError(WEException.GetStandardMessage(ex, iter, Result.PRNGSeed), this.ClassName);
                return;
            }
        }//compute
    }
}










