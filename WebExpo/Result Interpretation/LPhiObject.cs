namespace Zygotine.WebExpo
{
    internal class LPhiObject
    {
        public double[] R { get; internal set; }
        public double[] R2 { get; internal set; }

        public LPhiObject(int i)
        {
            R = new double[i];
            R2 = new double[i];
        }
    }
}
