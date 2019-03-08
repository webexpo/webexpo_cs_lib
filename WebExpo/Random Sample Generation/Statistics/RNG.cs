namespace Zygotine.Statistics
{
    public static class RNG
    {
        public static uint LastSeedUsed { get; private set; }
        public static Numerics.MersenneTwister MTPrng { get; private set; }
        static RNG()
        {
            MTPrng = new Numerics.MersenneTwister();
            SetSeed();
        }

        public static uint GetNewSeed()
        {
            System.Random rnd = new System.Random();
            return (uint)rnd.Next() + (uint)rnd.Next();
        }

        public static uint SetSeed()
        {
            uint seed = GetNewSeed();
            return RNG.SetSeed(seed);
        }

        public static uint SetSeed(uint s)
        {
            MTPrng.init_genrand(s);
            RNG.LastSeedUsed = s;
            return s;
        }
    }
}
