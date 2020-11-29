namespace Zygotine.WebExpo
{
    using System;
    using System.Collections.Generic;
    using System.Reflection;
    using Zygotine.Util;
    public class BWModelParameters_UniformPriorOnSDs : AbstractBWModelParameters
    {
        //SigmaBetweenRange and SigmaWithinRange 
        // are relevant only when UseUniformPriorOnSds = true
        public double[] SigmaBetweenRange { get; set; } = new double[0];
        public double[] SigmaWithinRange { get; set; } = new double[0];

        protected BWModelParameters_UniformPriorOnSDs(bool logNormalDstrn) : base(true, logNormalDstrn)
        {
        }

        public BWModelParameters_UniformPriorOnSDs(bool logNormalDstrn, double muOverallLower, double muOverallUpper, double initMuOverall, double initSigmaWithin, double[] sigmaBetweenRange, double[] sigmaWithinRange) : this(logNormalDstrn)
        {
            this.MuOverallLower = muOverallLower;
            this.MuOverallUpper = muOverallUpper;
            this.InitSigmaWithin = initSigmaWithin;
            this.InitMuOverall = initMuOverall;
            this.SigmaBetweenRange = sigmaBetweenRange;
            this.SigmaWithinRange = sigmaWithinRange;
        }

        public static BWModelParameters_UniformPriorOnSDs GetDefaults(bool logNormalDstrn)
        {
            BWModelParameters_UniformPriorOnSDs parameters = new BWModelParameters_UniformPriorOnSDs(logNormalDstrn);

            // les valeurs de InitMuOverall et de InitSigmaWithin doivent être définis
            // modifié le 31 mai 2017 (voir jl : Paramètres par défaut pour les fonctions de McGill May2017.docx)
            if (logNormalDstrn)
            {
                parameters.MuOverallLower = -100;
                parameters.MuOverallUpper = 100;
                parameters.InitMuOverall = 0.3;
                parameters.InitSigmaWithin = 1.0;

                parameters.SigmaBetweenRange = Tools.Combine(0, 2.3);
                parameters.SigmaWithinRange = Tools.Combine(0.095, 2.3);
            }
            else
            {
                parameters.MuOverallLower = 40;
                parameters.MuOverallUpper = 125;
                parameters.InitMuOverall = 85;
                parameters.InitSigmaWithin = 5;

                parameters.SigmaBetweenRange = Tools.Combine(0, 20);
                parameters.SigmaWithinRange = Tools.Combine(0.1, 20);
            }

            return parameters;
        }

        protected override double GetSigmaWithin()
        {
            return Tools.Mean(this.SigmaWithinRange);
        }

        public IEnumerator<PropertyInfo> GetEnumerator()
        {
            foreach (var property in typeof(BWModelParameters_UniformPriorOnSDs).GetProperties())
            {
                yield return property;
            }
        }
    }
}
