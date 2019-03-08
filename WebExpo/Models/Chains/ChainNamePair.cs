
namespace Zygotine.WebExpo
{
    public class ChainNamePair
    {
        public string Burnin { get; private set; }
        public string Sample { get; private set; }

        internal ChainNamePair(string burnin, string sample)
        {
            Burnin = burnin;
            Sample = sample;
        }
    }
}
