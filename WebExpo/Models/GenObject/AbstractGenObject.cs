namespace Zygotine.WebExpo
{
    using System.Collections.Generic;

    internal abstract class GenObject
    {
        internal SGNFnAParam A { get; set; }
        internal bool LogNormalDistrn { get; set; }
        internal double MathLowerLimit { get; set; } = 1E-6;
        internal bool PotentiallyBimodal { get; set; }
        internal double[] Range { get; set; } = new double[0];
        internal bool ThroughCV { get; set; } = false;
        internal bool ThroughSD { get; set; } = false;
        internal abstract double F(double x, SGNFnAParam A);
        internal abstract double LogF(double x, SGNFnAParam A);
        internal abstract double LogFPrime(double x, SGNFnAParam A);
        internal abstract double LogFSecond(double x, SGNFnAParam A);
        internal abstract List<double> LogFInvRemote(double x, SGNFnAParam A);
        internal abstract double[] Start(SGNFnAParam A);
        internal virtual string AAsString(SGNFnAParam A)
        {
            return "";
        }

        /*
         * Pour Sigma_Gen Objects
         */

        internal static GenObject GetSigmaGenObject(int N, double lNormMu, double lNormSd)
        {
            return new Sigma_LM_GenObject(N, lNormMu, lNormSd);
        }

        internal static GenObject GetSigmaTruncatedDataGenObject(int N)
        {
            return new SigmaTruncatedData_GenObject(N);
        }

        internal static GenObject GetSigmaTruncatedDataLNormGenObject(int N, double lNormMu, double lNormSd)
        {
            return new SigmaTruncatedDataLNorm_GenObject(N, lNormMu, lNormSd);
        }

        private static GenObject GetMuTruncatedDataGenObject(int n)
        {
            return new MuTruncatedData_GenObject(n);
        }

        internal static GenObject GetSigmaGenObject(int N)
        {
            // on a pas les infos sur lNormMu ni sur lNormSd ...
            return new Sigma_woLM_GenObject(N);
        }

        /*
         * Pour ME_Gen Objects
         */

        internal static GenObject GetMeGenObject(MeasurementError me, bool logNDstrn, int N)
        {
            if (me.ThroughCV)
            {
                return new ME_CV_GenObject(N, logNDstrn);
            }
            else
            {
                // me.ThroughSD est vrai
                return null;
            }
        }
    }

}
