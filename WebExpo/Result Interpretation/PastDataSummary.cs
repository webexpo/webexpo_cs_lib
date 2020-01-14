namespace Zygotine.WebExpo
{
    using System.Text.RegularExpressions;
    using System;

    public class PastDataSummary
    {
        public bool Defined { get; set; }
        public int N { get; set; }
        public double SD { get; set; }
        public double Mean { get; set; }
        private readonly static PastDataSummary emptyObject;
        public double NS2 { get; private set; }
        public string ClassName { get; private set; }
        internal Messages Messages { get; private set; } = new Messages();
        public static PastDataSummary EmptyObject { get { return emptyObject; } }

        static PastDataSummary()
        {
            emptyObject = new PastDataSummary();
        }

        public PastDataSummary()
        {
            ClassName = this.GetType().Name;
            this.SD = double.NaN;
            this.Mean = double.NaN;
            this.N = 0;
            this.Defined = false;
        }

        public PastDataSummary(string def) : this()
        {
            def = Regex.Replace(def, @" ", "");
            string[] parameters = def.Split(fieldSeparators);
            if (parameters.Length < 3)
            {
                this.Messages.AddError(string.Format("Invalid data. 3 parameters are expected (mean, sd, N).", def), this.ClassName);
            }
            else
            {
                try
                {
                    this.Mean = double.Parse(parameters[0]);
                    this.SD = double.Parse(parameters[1]);
                    this.N = int.Parse(parameters[2]);
                    this.Defined = true;
                    this.SetParameters();
                }
                catch
                {
                    this.Messages.AddError(string.Format("Error while parsing {0}.", def), this.ClassName);
                }
            }
        }

        private void SetParameters()
        {
            this.NS2 = (this.N - 1) * (this.SD * this.SD);
        }

        public PastDataSummary(double mean, double sd, int n) : this()
        {
            if (double.IsNaN(mean) || double.IsNaN(sd) || n < 1)
            {
                //TODO exception
            }

            this.Defined = true;
            this.N = n;
            this.SD = sd;
            this.Mean = mean;
            SetParameters();
            
        } //# end of past.data.summary

        public override string ToString()
        {

            return string.Format("{0:R}; {1:R}; {2}",
                this.Mean,
                this.SD,
                this.N);
        }

        public string ToString(bool verbose = true)
        {
            if (!verbose)
            {
                return this.ToString();
            }
            else
            {
                return string.Format("mean ={0:R}, sd={1:R}, N={2}, ns2={3:R}, used={4:R}",
                    this.Mean,
                    this.SD,
                    this.N,
                    this.NS2,
                    this.Defined);
            }
        }

        internal string ForRLanguage(string listName)
        {
            if (this == PastDataSummary.EmptyObject)
            {
                return string.Format("{0}$pastData <- list(mean=double(0), sd=double(0), n=integer(0))", listName);
            }
            else
            {
                return string.Format("{0}$pastData <- list(mean={1}, sd={2}, n={3})", listName, this.Mean, this.SD, this.N);
            }
        }

        private static char[] fieldSeparators = new char[] { '\t', ';' };

    }
}
