namespace Zygotine.WebExpo
{
    public class Message
    {
        public static int sequence = 0;
        public string Text { get; private set; }
        public Severity Level { get; private set; }
        public string Source { get; private set; }
        public long Sequence { get; private set; }
        public enum Severity { Ok = 0, Warning, Error, Info };

        public Message( string message, Severity lvl, string src)
        {
            this.Text = message;
            this.Level = lvl;
            this.Source = src;
            this.Sequence = Message.sequence++;
        }

        public override string ToString()
        {
            return string.Format("{0} - {1}: {2} From: {3}.", this.Sequence, this.Level.ToString(), this.Text, this.Source);
        }
    }
}
