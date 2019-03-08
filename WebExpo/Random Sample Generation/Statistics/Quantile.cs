namespace Zygotine.Statistics
{
    using System;
    using System.Collections.Generic;
    using System.Linq;
    using Zygotine.Util;

    /* 
     * R file: quantile.R
     * subset of quantile.defaultl function
     */
    public class Quantile
    {
        private const double eps = 2.22045E-14;
        private static double[] defaultProbs = new double[] { .025, .05, .1, .2, .3, .4, .5, .6, .7, .8, .9, .95, .975 };
        double[] probs = null;

        public double[] Probabilities
        {
            get
            {
                if (probs == null) return probs;
                return (double[])probs.Clone();
            }
        }

        public Quantile() : this(Quantile.defaultProbs)
        {
        }

        public Quantile(double prob) : this(new double[1] { prob })
        {
        }

        public Quantile(double[] probs)
        {
            bool probsOk = probs != null;

            if (probsOk)
            {
                this.probs = probs.Where(x => x.IsFinite()).ToArray();
                probsOk = probsOk && (this.probs.Length > 0);
            }


            if (probsOk)
            {
                foreach (double d in probs)
                {
                    if ((d < -eps) || (d > 1 + eps))
                    {
                        probsOk = false;
                        break;
                    }
                }
            }
            if (!probsOk)
                throw new ArgumentException("Invalid parameter 'probs'");

            Array.Sort(this.probs);

        }

        public double[] Compute(double[] x)
        {
            // attention aucune prise en compte des NaN (incluant les NA)
            double[] qs;
            int n = x.Length;

            int np = probs.Length;

            if ((n > 0) && (np > 0)) // np est nécessairement > 0, voir la validation faite dans le constructeur.
            {
                double[] index = probs.Multiply(n - 1);

                int[] lo = index.Apply(Math.Floor).ToInt32();
                int[] hi = index.Apply(Math.Ceiling).ToInt32();

                double[] y = new double[n];
                Array.Copy(x, y, n);
                Array.Sort(y);
                qs = y.SubArray(lo);
                List<int> iTmp = new List<int>();
                for (int k = 0; k < index.Length; k++)
                {
                    if (index[k] > lo[k])
                        iTmp.Add(k);
                }

                foreach (int k in iTmp)
                {
                    double h = (index[k] - lo[k]);
                    qs[k] = (1 - h) * qs[k] + h * y[hi[k]];
                }
            }
            else
            {
                qs = Tools.GetConstantArray(np, Tools.NA);
            }
            return qs;
        }
    }
}
