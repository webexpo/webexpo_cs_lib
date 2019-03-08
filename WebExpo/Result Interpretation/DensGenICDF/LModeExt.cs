namespace Zygotine.WebExpo
{
    using Zygotine.Util;

    public class LModeExt : LMode
    {
        internal double[] Remote { get; private set; }
        internal double[] Range { get; private set; }

        private LModeExt()
        {
        }

        public LModeExt(LMode mode, double remoteLeft, double remoteRight) : base(mode.X, mode.H, mode.Found)
        {
            this.Remote = Tools.Combine(remoteLeft, remoteRight);
        }

        public override string ToString()
        {
            return string.Format("{0}, Remote={1}", base.ToString(), Tools.Show(this.Remote));
        }
    }
}
