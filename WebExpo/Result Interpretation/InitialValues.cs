namespace Zygotine.WebExpo
{
    using System.Collections.Generic;

    public class InitialValues
    {
        public double Mu { get; private set; }
        public double SigmaWithin { get; private set; }
        public IEnumerable<double> NormalizedY { get; private set; } = new double[0];

        internal InitialValues(double mu, double sigma, IEnumerable<double> normalizedY = null)
        {
            this.Mu = mu;
            this.SigmaWithin = sigma;
            if (normalizedY != null)
            {
                this.NormalizedY = normalizedY;
            }
        }
    }
}
