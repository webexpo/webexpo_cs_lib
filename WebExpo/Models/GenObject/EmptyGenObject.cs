namespace Zygotine.WebExpo
{
    using System.Collections.Generic;
    using Zygotine.Util;

    internal class EmptyGenObject : GenObject
    {
        private readonly static EmptyGenObject emptyInstance = new EmptyGenObject();

        public static EmptyGenObject Instance { get { return emptyInstance; } }

        internal override double F(double sigma, SGNFnAParam A)
        {
            return Tools.ND;
        }

        internal override double LogF(double sigma, SGNFnAParam A)
        {
            return Tools.ND;
        }

        internal override double LogFPrime(double sigma, SGNFnAParam A)
        {
            return Tools.ND;
        }

        internal override double LogFSecond(double sigma, SGNFnAParam A)
        {
            return Tools.ND;
        }

        internal override List<double> LogFInvRemote(double target, SGNFnAParam A)
        {
            return new List<double>();
        } //# end of log.f.inv.remote

        internal override double[] Start(SGNFnAParam A)
        {
            return new double[0];
        }

        private EmptyGenObject()
        {
        }
    }

}
