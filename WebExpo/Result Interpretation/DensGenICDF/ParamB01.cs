
namespace Zygotine.WebExpo
{
    public class ParamB01
    {
        internal double X { get; private set; }
        internal double H { get; private set; }
        internal double Hp { get; private set; }
        internal bool Cntn { get; private set; }


        internal ParamB01(double x, double h, double hp)
        {
            this.X = x;
            this.H = h;
            this.Hp = hp;
        }

        internal ParamB01(double x, double h, double hp, bool cntn) : this(x,h,hp)
        {
            this.Cntn = cntn;
        }

    }
}
