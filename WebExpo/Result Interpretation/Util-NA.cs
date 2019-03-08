namespace Zygotine.Util
{
    using System;

    public static partial class Tools
    {
        public static readonly double NA;
        public static readonly double ND;

        static Tools()
        {
            byte[] b = BitConverter.GetBytes(double.NaN);
            // on a dans notre environnement 0,0,0,0,0,0,248,255

            b[6] = 252;
            NA = BitConverter.ToDouble(b, 0);
            // on aura le bit du signe à 1, 
            // les 11 bits de l'exposant à 1;
            // et les premier et second bits de la mantisse à 1. NaN n'a que le premier bit levé.
            b[6] = 254;
            ND = BitConverter.ToDouble(b, 0);
            // on aura le bit du signe à 1, 
            // les 11 bits de l'exposant à 1;
            // et les 3 premiers à 1. NaN n'a que le premier bit levé.
        }

        
    }
}
