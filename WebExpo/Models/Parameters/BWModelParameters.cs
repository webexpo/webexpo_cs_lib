namespace Zygotine.WebExpo
{
    using System;
    using Zygotine.Util;

    public class BWModelParameters : AbstractBWModelParameters
    {
        // LogSigmaBetweenMu, LogSigmaBetweenPrec, LogSigmaWithinMu, LogSigmaWithinPrec
        // are relevant only when UseUniformPriorOnSds = false
        public double LogSigmaBetweenMu { get; set; } = Tools.ND;
        public double LogSigmaBetweenPrec { get; set; } = Tools.ND;
        public double LogSigmaWithinMu { get; set; } = Tools.ND;
        public double LogSigmaWithinPrec { get; set; } = Tools.ND;

        protected BWModelParameters(bool logNormalDstrn) : base(false, logNormalDstrn)
        {
        }

        public BWModelParameters(bool logNormalDstrn, double muOverallLower, double muOverallUpper, double initMuOverall, double initSigmaWithin) : this(logNormalDstrn)
        {
            this.MuOverallLower = muOverallLower;
            this.MuOverallUpper = muOverallUpper;
            this.InitSigmaWithin = initSigmaWithin;
            this.InitMuOverall = initMuOverall;
        }

        public void SetLogSigmaBetweenMu(double value)
        {
            this.LogSigmaBetweenMu = value;
        }

        public void SetLogSigmaBetweenPrec(double value)
        {
            this.LogSigmaBetweenPrec = value;
        }

        public void SetLogSigmaWithinMu(double value)
        {
            this.LogSigmaWithinMu = value;
        }

        public void SetLogSigmaWithinPrec(double value)
        {
            this.LogSigmaWithinPrec = value;
        }

        public static BWModelParameters GetDefaults(bool logNormalDstrn)
        {
            BWModelParameters parameters = new BWModelParameters(logNormalDstrn);

            // si les valeurs de InitMuOverall et/ou de InitSigmaWithin sont laissées indéfinies
            // la méthode DefaultInits pourra alors être appelé. Sinon l'usager doit les définir explicitement.
            // L'appel à DefaultInits se fera dans le constructeur du modèle, et elle ne sera appelé 
            // que si InitMuOverall ou InitSigmaWithin n'est pas défini.

            if (logNormalDstrn)
            {
                parameters.MuOverallLower = -20;
                parameters.MuOverallUpper = 20;

                parameters.LogSigmaBetweenMu = -0.8786;
                parameters.LogSigmaBetweenPrec = 1.634;
                parameters.LogSigmaWithinMu = -0.4106;
                parameters.LogSigmaWithinPrec = 1.9002;

                parameters.InitMuOverall  = Math.Log(0.3);
                parameters.InitSigmaWithin = Math.Log(2.5);
            }
            else
            {
                parameters.MuOverallLower = 40;
                parameters.MuOverallUpper = 125;

                parameters.LogSigmaBetweenMu = 1.098612;// #(GM = 3) 
                parameters.LogSigmaBetweenPrec = 1.191059; //# (GSD=2.5)
                parameters.LogSigmaWithinMu = 1.098612; // #(GM = 3) 
                parameters.LogSigmaWithinPrec = 1.191059; // # (GSD=2.5)

                parameters.InitMuOverall  = 85;
                parameters.InitSigmaWithin = 3;
            }

            return parameters;
        }

        protected override double GetSigmaWithin()
        {
            return Math.Exp(this.LogSigmaWithinMu);
        }
    }
}
