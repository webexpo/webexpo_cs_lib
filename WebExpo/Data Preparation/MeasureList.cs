namespace Zygotine.WebExpo
{
    using System;
    using System.Collections.Generic;
    using System.Globalization;
    using System.Linq;
    using System.Text;
    using System.Text.RegularExpressions;

    public class MeasureList
    {
        public double OEL { get; set; }
        private static readonly Regex rgxME = new Regex(@"(sd|cv)\(.*\)", RegexOptions.IgnoreCase | RegexOptions.Compiled);

        #region Constructor
        private MeasureList()
        {
            this.ClassName = this.GetType().Name;
            foreach (Measure.MeasureType ctype in Enum.GetValues(typeof(Measure.MeasureType)))
                measuresByCensoringType.Add(ctype, new List<Measure>());
        }

        public MeasureList(string s, double oel = 99.9) : this()
        {
            this.OEL = oel;

            string[] records = s.Replace(" ", "").Split(recordSeparators, StringSplitOptions.RemoveEmptyEntries);
            if (records.Length == 0)
            {
                this.Error.AddError(string.Format("No measure found in string '{0}'!", s), this.ClassName);
                return;
            }

            for (int i = 0; i < records.Length; i++)
            {
                records[i] = records[i].Replace(@" ", "");
            }

            bool hasME = rgxME.IsMatch(records[0]);
            if (hasME)
            {
                this.AddME(records[0]);
                if (this.ME.Equals(MeasurementError.None))
                {
                    // Si l'insertion ne s'est pas faite, alors il y a erreur de syntaxe.
                    this.Error.AddError(string.Format("Could not build measurement error from string '{0}'!", records[0]), this.ClassName);
                    return;
                }
            }

            for (int i = hasME ? 1 : 0; i < records.Length; i++)
            {
                Measure m = ParseMeasure(records[i]);
                if (m != null)
                {
                    this.Add(m);
                }
                else
                {
                    if (this.Error.Level == Message.Severity.Error)
                    {
                        //Une erreur s'est produite, on arrète la construction de la liste.
                        return; // en js on poursuit, pour tenter de documenter le plus d'erreurs possibles.
                    }

                 }
            }
        }

	// Convenience constructor added
        public MeasureList(Array measures, double oel) : this(ConvertMeasuresToString(measures), oel)
        {
        }

	// Convenience constructor added
        public MeasureList(Array measures, double oel, double[] measErrRange) : this(ConvertMeasErrToString(measErrRange) + ConvertMeasuresToString(measures), oel)
        {

        }


        private static String ConvertMeasuresToString(Array measures)
        {
            Object[] meas = new object[measures.Length];
            measures.CopyTo(meas, 0);
            string str = String.Join("|", meas);
            return str;
        }

        private static String ConvertMeasErrToString(double[] measErrVarCoeffRange)
        {
            String str = "";
            if (measErrVarCoeffRange != null && measErrVarCoeffRange.Length == 2 )
            {
                str = "cv(" + ShowDouble(measErrVarCoeffRange[0] / 100) + "~" + ShowDouble(measErrVarCoeffRange[1] / 100) + ")|";
            }
            return str;
        }

        private static string ShowDouble(double d)
        {
            return ToSignificantDigits(d, 3);
        }

        private static string ToSignificantDigits(double value, int significant_digits)
        {
            // Use G format to get significant digits.
            // Then convert to double and use F format.
            string format1 = "{0:G" + significant_digits.ToString() + "}";
            string result = Convert.ToDouble(
                String.Format(format1, value)).ToString("F99");

            // Rmove trailing 0s.
            result = result.TrimEnd('0');

            // Rmove the decimal point and leading 0s,
            // leaving just the digits.
            string test = result.Replace(",", ".").Replace(".", "").TrimStart('0');

            // See if we have enough significant digits.
            if (significant_digits > test.Length)
            {
                // Add trailing 0s.
                result += new string('0', significant_digits - test.Length);
            }
            else
            {
                // See if we should remove the trailing decimal point.
                if (/*(significant_digits < test.Length) &&*/ result.EndsWith(","))
                    result = result.Substring(0, result.Length - 1);
            }

            return result;
        }

        public void RemoveMeasurementError()
        {
            if (this.MEAny)
            {
                this.ME = MeasurementError.None;
                this.ME.Any = false;
                this.METhroughCV = false;
                this.METhroughSD = false;
            }
        }

        #endregion

        #region Add methods

        private void AddME(string me)
        {
            string x = Regex.Replace(me, @" ", "");
            if (x.Length < 5)
            {
                return;
            }

            string meType = x.Substring(0, 2).ToLower();
            if (meType != "sd" && meType != "cv")
            {
                return;
            }

            if (x[2] != '(')
            {
                return;
            }

            if (x.Last() != ')')
            {
                return;
            }

            string str = x.Substring(3, x.Length - 4);
            if (str == string.Empty)
            {
                return;
            }

            string[] numbers = str.Split(new char[] { '~' }, StringSplitOptions.RemoveEmptyEntries);
            if ((numbers.Length == 0) || (numbers.Length > 2))
            {
                return;
            }

            double a, b;
            bool test;
            test = ParseDouble(numbers[0], out a);
            if (!test)
            {
                return;
            }

            b = a;
            if (numbers.Length == 2)
            {
                test = ParseDouble(numbers[1], out b);
                if (!test)
                {
                    return;
                }
            }

            if (a > b)
            {
                Error.AddError(string.Format("({0}, {1}) do not represent a valid range for the measurement error.", a, b), this.ClassName);
                return;
            }

            this.ME = (meType == "sd") ? MeasurementError.NewSDError(a, b) : MeasurementError.NewCVError(a, b);
            return;
        }

        private void Add(Measure m)
        {
            m.Ordinal = this.measuresList.Count;
            measuresList.Add(m);
            if (m.Censoring != Measure.MeasureType.Uncensored)
            {
                HasCensoredMeasure = true;
            }

            //Update collections
            measuresByCensoringType[m.Censoring].Add(m);

            if (m.Worker != null)
            {
                if (WorkersByTag.ContainsKey(m.Worker.Tag))
                {
                    m.Worker = WorkersByTag[m.Worker.Tag]; // will be the same object, the same worker : tag must be unique.
                    measuresByWorker[m.Worker].Add(m);
                }
                else
                {
                    WorkersByTag.Add(m.Worker.Tag, m.Worker);
                    List<Measure> lm = new List<Measure>();
                    lm.Add(m);
                    measuresByWorker.Add(m.Worker, lm);
                }
            }
            else
            {
                IsWorkerInfoComplete = false;
            }
        }
        #endregion

        #region ToString
        public override string ToString()
        {
            return ToString(RecordSeparator.Pipe, FieldSeparator.Semicolon);
        }

        public string ToString(RecordSeparator rs, FieldSeparator fs)
        {

            if (this.Count == 0)
            {
                return "";
            }

            string rsStr = GetRecordSeparator(rs);
            string fsStr = GetFieldSeparator(fs);
            int length = this.Count;
            if (this.MEAny)
            {
                length++;
            }


            string[] rep = new string[length];
            int position = 0;

            // L'erreur de mesure
            if (this.MEAny)
            {
                rep[position] = string.Format("{0}({1}, {2})",
                    ME.ThroughCV ? "CV" : "SD",
                    this.ME.Minimum,
                    this.ME.Maximum);
                position++;
            }

            StringBuilder sb = new StringBuilder(measuresList[0].ToString());
            for (int i = 0; i < measuresList.Count; i++)
            {
                rep[position] = measuresList[i].ToString();
                position++;
            }

            return string.Join(rsStr, rep);
        }

        private Tuple<string, string, string> ForRLanguageIntervals(string listName, List<Measure> intervals)
        {
            //RData$uncensored < -c(0.406, 0.597, 0.718312913)
            //RData$rightCensored < -c(0.2, 0.2, 0.2, 0.2)
            //RData$leftCensored < -c(0.6)
            //RData$intervalCensored$left < -c(RData$intervalCensored$right < -0.4, 0.4)
            //c(0.6, 0.6)

            StringBuilder sbA = new StringBuilder();
            StringBuilder sbB = new StringBuilder();
            StringBuilder sbW = new StringBuilder();
            if (intervals.Count == 0)
            {
                sbA.AppendFormat("{0}${1}$left$values <- numeric(0)\n", listName, Measure.MeasureType.Interval);
                sbB.AppendFormat("{0}${1}$right$values <- numeric(0)\n", listName, Measure.MeasureType.Interval);
                sbW.AppendFormat("{0}${1}$worker <- character(0)\n", listName, Measure.MeasureType.Interval);
            }
            else
            {
                sbA.AppendFormat("{0}${1}$left$values <- c(", listName, Measure.MeasureType.Interval);
                sbB.AppendFormat("{0}${1}$right$values <- c(", listName, Measure.MeasureType.Interval);
                sbW.AppendFormat("{0}${1}$worker <- c(", listName, Measure.MeasureType.Interval);

                Measure currentMeasure = intervals[0];

                sbA.AppendFormat("{0:R}", currentMeasure.A);
                sbB.AppendFormat("{0:R}", currentMeasure.B);

                if (currentMeasure.Worker == null)
                    sbW.Append("NA");
                else
                    sbW.AppendFormat("\"{0}\"", currentMeasure.Worker.Tag);


                for (int i = 1; i < intervals.Count; i++)
                {
                    currentMeasure = intervals[i];
                    sbA.AppendFormat(", {0}", currentMeasure.A);
                    sbB.AppendFormat(", {0}", currentMeasure.B);
                    if (currentMeasure.Worker == null)
                        sbW.Append(", NA");
                    else
                        sbW.AppendFormat(", \"{0}\"", currentMeasure.Worker.Tag);
                }

                sbA.Append(")\n");
                sbB.Append(")\n");
                sbW.Append(")\n");
            }

            return new Tuple<string, string, string>(sbA.ToString(), sbB.ToString(), sbW.ToString());
        }

        internal string ForRLanguage(string listName)
        {
            const string ME_CV_RANGE = "me_CVRange";
            const string ME_SD_RANGE = "MESDRange";
            const string NUMERIC0 = "numeric(0)";

            StringBuilder sb = new StringBuilder();
            StringBuilder sbA;
            StringBuilder sbW;

            sb.AppendFormat("{0} <- list()\n", listName);

            // L'erreur de mesure
            if (!this.MEAny)
            {
                sb.AppendFormat("{0}${1} <- {2}\n", listName, ME_CV_RANGE, NUMERIC0);
                sb.AppendFormat("{0}${1} <- {2}\n", listName, ME_SD_RANGE, NUMERIC0);
            }
            else
            {
                sb.AppendFormat("{0}${1} <- c({2:R}, {3:R})\n",
                    listName,
                    ME.ThroughCV ? ME_CV_RANGE : ME_SD_RANGE,
                    this.ME.Minimum,
                    this.ME.Maximum);

                sb.AppendFormat("{0}${1} <- {2}\n",
                    listName,
                    ME.ThroughCV ? "MESDRange" : "me_CVRange",
                    NUMERIC0);
            }

            sb.AppendFormat("{0}$WorkerInfoComplete <- {1}", listName, IsWorkerInfoComplete ? "T" : "F").AppendLine();
            foreach (Measure.MeasureType ct in Enum.GetValues(typeof(Measure.MeasureType)))
            {
                if (ct == Measure.MeasureType.Interval)
                {
                    Tuple<string, string, string> t = ForRLanguageIntervals(listName, measuresByCensoringType[ct]);
                    sb.Append(t.Item1).Append(t.Item2).Append(t.Item3);
                    continue;
                }
                sbA = new StringBuilder();
                sbW = new StringBuilder();
                sbA.AppendFormat("{0}${1}$values <-", listName, ct.ToString());
                sbW.AppendFormat("{0}${1}$worker <-", listName, ct.ToString());

                if (measuresByCensoringType[ct].Count == 0)
                {
                    sbA.AppendFormat(" {0}\n", NUMERIC0);
                    sbW.Append(" character(0)\n");
                }
                else
                {
                    sbA.Append("c(");
                    sbW.Append("c(");
                    Measure currentMeasure = measuresByCensoringType[ct][0];
                    sbA.AppendFormat("{0:R}", currentMeasure.A);
                    if (currentMeasure.Worker == null)
                        sbW.Append("NA");
                    else
                        sbW.AppendFormat("\"{0}\"", currentMeasure.Worker.Tag);

                    for (int i = 1; i < measuresByCensoringType[ct].Count; i++)
                    {
                        currentMeasure = measuresByCensoringType[ct][i];
                        sbA.AppendFormat(", {0:R}", currentMeasure.A);

                        if (currentMeasure.Worker == null)
                            sbW.Append(", NA");
                        else
                            sbW.AppendFormat(", \"{0}\"", currentMeasure.Worker.Tag);
                    }

                    sbA.Append(")\n");
                    sbW.AppendFormat(")\n");
                }
                sb.Append(sbA).Append(sbW);
            }
            string sbStr = sb.ToString();
            return sbStr;
        }
        #endregion

        #region Indexer

        public Measure this[int index]
        {
            get { return measuresList[index]; }
        }
        #endregion

        #region fields
        // Some collections
        internal Dictionary<Measure.MeasureType, List<Measure>> measuresByCensoringType = new Dictionary<Measure.MeasureType, List<Measure>>();
        internal SortedList<string, Worker> WorkersByTag = new SortedList<string, Worker>();
        private Dictionary<Worker, List<Measure>> measuresByWorker = new Dictionary<Worker, List<Measure>>();

        // measuresList donne les mesures dans l'ordre où elles ont été entrées.
        private List<Measure> measuresList = new List<Measure>();
        #endregion

        public void StandardizeObservations(double oel)
        {
            foreach (Measure m in measuresList)
            {
                if (!double.IsNaN(m.A))
                {
                    m.A = m.A / oel;
                }
                if (!double.IsNaN(m.B))
                {
                    m.B = m.B / oel;
                }
            }
        }

        #region properties
        public bool Valid { get { return Count > 5; } }
        public bool WorkerValid
        {
            get
            {
                if (!this.IsWorkerInfoComplete)
                {
                    return false;
                }
                if (this.WorkerCount < 2)
                {
                    return false;
                }

                int n = 0;
                foreach (KeyValuePair<Worker, List<Measure>> kv in this.measuresByWorker)
                {
                    if (kv.Value.Count < 4)
                    {
                        return false;
                    }
                    n += kv.Value.Count;
                }

                return n == this.Count;
            }
        }

        public int Count { get { return measuresList.Count; } }
        public bool HasCensoredMeasure { get; private set; } = false;
        public int N { get { return Count; } }
        public int WorkerCount { get { return WorkersByTag.Count; } }
        public bool IsWorkerInfoComplete { get; private set; } = true;
        public string[] WorkerTags { get { return WorkersByTag.Keys.ToArray(); } }
        public MeasurementError ME { get; private set; } = MeasurementError.None;
        public bool METhroughSD { get { return this.ME.ThroughSD; } private set { this.ME.ThroughSD = value; } }
        public bool METhroughCV { get { return this.ME.ThroughCV; } private set { this.ME.ThroughCV = value; } }
        public bool MEAny { get { return this.ME.Any; } private set { this.ME.Any = value; } }
        public Dictionary<Worker, List<Measure>> MeasuresByWorker { get { return measuresByWorker; } }
        public Messages Error { get; private set; } = new Messages();
        public string ClassName { get; private set; }
        #endregion

        #region fields and records separators
        public enum RecordSeparator { CarriageReturn = 0, Linefeed, CarriageReturnLinefeed, Pipe }
        public enum FieldSeparator { Semicolon, Tab }
        // for output to the user
        private static string[] recordSeparators = new string[] { "\r", "\n", "\r\n", "|" };
        private static string[] fieldSeparators = new string[] { ";", "\t" };
        public static string GetFieldSeparator(FieldSeparator fs) { return fieldSeparators[(int)fs]; }
        public static string GetRecordSeparator(RecordSeparator rs) { return recordSeparators[(int)rs]; }
        #endregion

        #region helper struct
        private struct PreMeasure
        {
            internal Measure.MeasureType measureType;
            internal double a;
            internal double b;
            internal bool hasError;
            internal PreMeasure(double a, double b)
            {
                this.measureType = Measure.MeasureType.Uncensored;
                this.a = a;
                this.b = b;
                this.hasError = false;
            }
        }
        #endregion

        #region private static methods Parsing

        private Measure ParseMeasureValue(string str)
        {
            // ! string.IsNullOrEmpty(str) is always true

            PreMeasure pm = new PreMeasure();
            pm.hasError = true;
            string strOrig = str;
            switch (str[0])
            {
                case '[':
                    if (str.Last() != ']')
                    {
                        break;
                    }

                    str = str.Substring(1, str.Length - 2);
                    if (str == string.Empty)
                    {
                        break;
                    }

                    string[] nlines = str.Split(new char[] { ',', '-' }, StringSplitOptions.RemoveEmptyEntries);
                    if (nlines.Length != 2)
                    {
                        this.Error.AddError(string.Format("Erreur survenue lors de lecture de mesure censurée '{0}'", strOrig), this.ClassName);
                        break;
                    }

                    pm.measureType = Measure.MeasureType.Interval;
                    pm.hasError = !ParseDouble(nlines[0], out pm.a);
                    if (pm.hasError)
                    {
                        break;
                    }

                    pm.hasError = !ParseDouble(nlines[1], out pm.b);
                    pm.hasError = pm.hasError || (pm.a > pm.b);
                    break;

                case '<':
                    pm.measureType = Measure.MeasureType.LessThan;
                    pm.hasError = !ParseDouble(str.Substring(1), out pm.a);
                    break;

                case '>':
                    pm.measureType = Measure.MeasureType.GreaterThan;
                    pm.hasError = !ParseDouble(str.Substring(1), out pm.a);
                    break;

                default:
                    pm.measureType = Measure.MeasureType.Uncensored;
                    pm.hasError = !ParseDouble(str, out pm.a);
                    break;

            }

            if (pm.hasError)
            {
                return null;
            }

            switch (pm.measureType)
            {
                case Measure.MeasureType.GreaterThan:
                    {
                        return new GTMeasure(pm.a);
                    }

                case Measure.MeasureType.LessThan:
                    {
                        return new LTMeasure(pm.a);
                    }

                case Measure.MeasureType.Interval:
                    {
                        return new IntervalMeasure(pm.a, pm.b);
                    }

                default:
                    {
                        return new UncensoredMeasure(pm.a);
                    }
            }
        }

        private Measure ParseMeasure(string s)
        {

            /* 
             * s décrit une mesure, et peut contenir un 'tag' identifiant un travailleur
             * Les 2 éléments d'information par | (piping)
             * L'ordre est le suivant : 
             *
             * 1) La mesure comme telle. 
             * 2) Le travailleur est une chaine dont la longueur est >=1 dont les blancs seront supprimés
             * 
             * La mesure comme telle est obligatoire. 
             * exemples 
             * [3,4]\tFrancois   
             * <4;Francois
             * >4\tFrancois       
             * <4.1;
             * <4                
             */

            Measure m = null;
            const string exceptionMsg = "MeasureList: Parse: Invalid syntax: \"{0}\".";

            string record = Regex.Replace(s, @" ", ""); // remove whitespaces
            if (record == string.Empty)
            {
                return null; // null, the calling program will ignore that 'empty' measure.
            }

            string[] fields = record.Split(fieldSeparators, StringSplitOptions.None);

            if (fields.Length == 0)
            {
                return null; // calling program shall ignore this 'empty' measure.
            }

            if (fields.Length <= 2)
            {
                // measure + worker tag
                if (fields[0].Length > 0)
                {
                    m = ParseMeasureValue(fields[0]);
                    if (m == null)
                    {
                        return null;
                    }
                    if (fields.Length > 1)
                    {
                        m.SetWorker(fields[1]);  // si fields[1] est vide ou n'est composé que d'espaces, m.Worker demeurera null.
                    }
                }
            }
            if (m == null)
            {
                this.Error.AddError(string.Format(exceptionMsg, s), this.ClassName);
                return null;
            }
            return m;
        }

        private bool ParseDouble(string s, out double value)
        {
            bool ok = double.TryParse(s, NumberStyles.Float, CultureInfo.InvariantCulture, out value);
            if (!ok)
            {
                s = s.Contains('.') ? s.Replace(".", ",") : s.Replace(",", ".");
                ok = double.TryParse(s, NumberStyles.Float, CultureInfo.InvariantCulture, out value);
            }
            if (!ok)
            {
                Error.AddError(string.Format("Error while parsing string '{0}'", s), this.ClassName);
            }
            return ok;
        }
        #endregion

        internal static MeasureList NewList(string s, bool jeromeStyle = true, bool discardCV = false)
        {
            if (!jeromeStyle)
            {
                //On fait alors l'hypothèse que c'est une liste 'normale' et on ne tient pas compte de discardCV
                return new MeasureList(s);
            }

            string[] records = s.Replace(" ", "").Split(recordSeparators, StringSplitOptions.RemoveEmptyEntries);
            if (records.Length == 0)
            {
                return new MeasureList(s);
            }

            string fs = GetFieldSeparator(FieldSeparator.Semicolon);
            string rs = GetRecordSeparator(RecordSeparator.Pipe);
            List<string> newList = new List<string>();
            for (int i = 0; i < records.Length; i++)
            {
                string currentMeasureValue = "";
                string currentWorker = "";

                if (records[i] == string.Empty)
                {
                    continue;
                }

                string[] fields = records[i].Split(fieldSeparators, StringSplitOptions.None);
                //3 champs possiblement, mesure, cv et travailleur
                if (fields.Length == 0)
                {
                    continue;
                }

                currentMeasureValue = fields[0];
                if (fields.Length >= 3)
                {
                    currentWorker = fields[2];
                }

                string m = currentMeasureValue;
                if (currentWorker != string.Empty)
                {
                    m += fs + currentWorker;
                }
                newList.Add(m);

            } // for (int i ...

            return new MeasureList(string.Join(rs, newList));
        }
    }
}
