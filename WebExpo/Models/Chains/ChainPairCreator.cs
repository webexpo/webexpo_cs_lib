
namespace Zygotine.WebExpo
{
    internal class ChainPairCreator
    {
        private const bool IsSample = true;


        private int _burninLength;
        private int _sampleLength;

        internal ChainPairCreator(int burnin, int sample)
        {
            _burninLength = burnin;
            _sampleLength = sample;
        }

        public ChainPair GetPair(string prefix)
        {
            return new ChainPair(this, prefix);
        }

        internal class ChainPair
        {
            public Chain Burnin { get; internal set; }
            public Chain Sample { get; internal set; }

            internal ChainPair(ChainPairCreator parent, string prefix)
            {
                ChainNamePair deux = GetChainNames(prefix);
                Burnin = new Chain(deux.Burnin, prefix, parent._burninLength, ! IsSample);
                Sample = new Chain(deux.Sample, prefix, parent._sampleLength, IsSample);
            }
        }

        public static ChainNamePair GetChainNames(string prefix)
        {
            return new ChainNamePair(prefix + "Burnin", prefix + "Sample");
        }
    }
}
