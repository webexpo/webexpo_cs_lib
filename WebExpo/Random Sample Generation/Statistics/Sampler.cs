namespace Zygotine.Statistics
{
    using System;
    using Zygotine.Statistics.Distribution;
    using Zygotine.Util;

    public class SamplerWithoutReplacement
    {

        private double[] weights;
        private int sampleSize;
        public SamplerWithoutReplacement(double[] weights, int sampleSize)
        {
            this.weights = weights;
            this.sampleSize = sampleSize;
        }

        /* R file: random.c
         * R fucntion name: ProbSampleNoReplace
         * On tire depuis 0..weights.Length-1 une suite d'entiers distincts dont la cardinalité est 'sampleSize' 
         */
        public int[] GetSample()
        {
            int[] perm = new int[weights.Length];
            int[] ans = new int[sampleSize];
            double rT, mass, totalMass;
            int i, j, k, n1;
            double[] normalizedWeights = new double[weights.Length];
            Array.Copy(weights, normalizedWeights, weights.Length);

            FixUpProb(normalizedWeights, ans.Length, false);
            for (i = 0; i < weights.Length; i++)
            {
                perm[i] = i;
            }

            Array.Sort(normalizedWeights, perm, Zygotine.Util.DescendingComparer.Desc);
            totalMass = 1.0;
            for (i = 0, n1 = weights.Length - 1; i < ans.Length; i++, n1--)
            {
                rT = totalMass * UniformDistribution.RUnif();
                mass = 0.0;
                for (j = 0; j < n1; j++)
                {
                    mass += normalizedWeights[j];
                    if (rT <= mass)
                    {
                        break;
                    }
                }
                ans[i] = perm[j];
                totalMass -= normalizedWeights[j];
                for (k = j; k < n1; k++)
                {
                    normalizedWeights[k] = normalizedWeights[k + 1];
                    perm[k] = perm[k + 1];
                }
            }
            return ans;
        }

        /* 
         * R file: random.c
         */
        private void FixUpProb(double[] p, int require_k, bool replace)
        {
            double sum;
            int i, npos;
            npos = 0;
            sum = 0.0;
            for (i = 0; i < p.Length; i++)
            {
                if (!p[i].IsFinite())
                {
                    // error(_("NA in probability vector"));
                }
                if (p[i] < 0)
                {
                    // error(_("non-positive probability"));
                }

                if (p[i] > 0)
                {
                    npos++;
                    sum += p[i];
                }
            }

            if (npos == 0 || (!replace && require_k > npos))
            {
                // error(_("too few positive probabilities"));
            }
            for (i = 0; i < p.Length; i++)
                p[i] /= sum;
        }
    }
}
