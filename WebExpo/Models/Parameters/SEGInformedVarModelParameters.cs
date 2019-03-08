namespace Zygotine.WebExpo
{
    using System;

    public class SEGInformedVarModelParameters : AbstractModelParameters
    {
        public double MuLower { get; set; }
        public double MuUpper { get; set; }
        public double LogSigmaMu { get; set; }
        public double LogSigmaPrec { get; set; }
        public double InitMu { get; set; }
        public double InitSigma { get; set; }

        private SEGInformedVarModelParameters()
        {
        }

        public SEGInformedVarModelParameters(
            double muLower,
            double muUpper,
            double logSigmaMu,
            double logSigmaPrec,
            double initMu, /* log(0.3) */
            double initSigma,/*log(2.5)*/
            bool logNormalDstrn = true) : base(logNormalDstrn)
        {
            this.MuLower = muLower;
            this.MuUpper = muUpper;
            this.LogSigmaMu = logSigmaMu;
            this.LogSigmaPrec = logSigmaPrec;
            this.InitMu = initMu;
            this.InitSigma = initSigma;
        }

        public static SEGInformedVarModelParameters GetDefaults(bool logNormalDstrn)
        {
            SEGInformedVarModelParameters p = new SEGInformedVarModelParameters();
            p.LogNormalDstrn = logNormalDstrn;
            // modifié le 31 mai 2017 (voir jl : Paramètres par défaut pour les fonctions de McGill May2017.docx)
            if (p.LogNormalDstrn)
            {
                p.MuLower = -20;
                p.MuUpper = 20;
                p.LogSigmaMu = -0.1744;
                p.LogSigmaPrec = 2.5523;
                p.InitMu = -1.2039728043259361; /* log(0.3) */
                p.InitSigma = 0.91629073187415511; /*log(2.5)*/
            }
            else
            {
                p.MuLower = 40;
                p.MuUpper = 125;
                p.LogSigmaMu = 1.0986122886681098; //#(GM = 3) // corrigé le 31 mai log(3)
                p.LogSigmaPrec = 1.191059; //(GSD=2.5)
                p.InitMu = 85;
                p.InitSigma = 3;
            }
            return p;
        }
    }
}
