namespace Zygotine.WebExpo
{
    using Zygotine.Util;

    public class UninformativeModelParameters : AbstractModelParameters
    {

        public double MuLower { get; set; }
        public double MuUpper { get; set; }
        public double[] SDRange { get; set; }
        public double InitMu { get; set; }
        public double InitSD { get; set; }

        private UninformativeModelParameters()
        {
        }

        public static UninformativeModelParameters GetDefaults(bool logNormalDstrn)
        {
            UninformativeModelParameters p = new UninformativeModelParameters();
            p.LogNormalDstrn = logNormalDstrn;
            if (p.LogNormalDstrn)
            {
                p.MuLower = -20;
                p.MuUpper = 20;
                p.SDRange = Tools.Combine(.095, 2.3);
                p.InitMu = -1.2039728043259361; // log(0.3);
                p.InitSD = 0.91629073187415511; // log(2.5)
            }
            else
            {
                p.MuLower = 40;
                p.MuUpper = 125;
                p.SDRange = Tools.Combine(.1, 20);
                p.InitMu = 85;
                p.InitSD = 3;
            }

            return p;
        }

        public override string ToString()
        {
            return string.Format("MuLower={0:R}, MuUpper=={1:R}, SDRange={2}, InitMu={3:R}, InitSD={4:R}", this.MuLower, this.MuUpper, this.SDRange.Show(), this.InitMu, this.InitSD);
        }

    }
}
