namespace Zygotine
{
    using System;
    using System.Linq;

    public class WEException : Exception
    {
        private const string defaultMsg = "Unspecified WebExpo exception!";

        public WEException() : base(defaultMsg)
        {
        }

        public WEException(string message) : base(message)
        {
        }


        private static string StackTraceFmt = "Class: {0}, {1}: {2}, File: {3} (line number={4})";

        public static string GetStackTrace(Exception e)
        {
            return FormatStackTrace(new System.Diagnostics.StackTrace(e, true));
        }

        private static string FormatStackTrace(System.Diagnostics.StackTrace st)
        {
            string[] q = st.GetFrames().Select(frame =>

                      string.Format(
                            StackTraceFmt,
                            frame.GetMethod().DeclaringType.Name, // Classe
                            frame.GetMethod().MemberType, // Méthode ou Constructeur
                            frame.GetMethod().Name, //Méthode
                            System.IO.Path.GetFileName(frame.GetFileName()), // Nom du fichier
                            frame.GetFileLineNumber())
                  ).ToArray(); // Numéro de la ligne
            return string.Join("\r\n", q);
        }

        public static string GetStandardMessage(Exception e, int iter, uint prngSeed)
        {
            return string.Format("Exception occured at iter={0} using rng seed of {1}: {2}", iter, prngSeed, e.Message) + "\r\nStack trace\r\n" + WEException.GetStackTrace(e);

        }
    }
}
