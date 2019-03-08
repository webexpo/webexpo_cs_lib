
namespace Zygotine.WebExpo
{
    using Zygotine.Util;

    public class LMode
    {
        internal double X { get; set; } = Tools.NA;
        internal double H { get; set; } = Tools.NA;
        internal bool Found { get; set; } = false;

        internal LMode()
        {
        }

        internal LMode(double x, double h, bool found = false)
        {
            X = x;
            H = h;
            Found = found;
        }

        public override string ToString()
        {
            return string.Format("X={0}, H={1}, Found={2}", this.X, this.H, this.Found);
        }
    }
}
