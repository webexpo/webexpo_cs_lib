namespace Zygotine.WebExpo
{
    using System.Collections.Generic;
    using System.Linq;

    public class Messages
    {
        public List<Message> lst { get; } = new List<Message>();
        public Message.Severity Level { get; private set; } = Message.Severity.Ok;

        internal Messages()
        {
        }

        internal void Add(Messages msgs)
        {
            foreach (Message msg in msgs.lst)
            {
                this.Add(msg);
            }
        }

        public List<Message> GetMessages(bool errorOnly = false)
        {
            List<Message> rep;
            if (errorOnly)
            {
                if (this.Level == Message.Severity.Error)
                {
                    rep = new List<Message>(); // une liste vide                    
                }
                else
                {
                    rep = lst.Where(e => e.Level == Message.Severity.Error).ToList();
                }
            }
            else
            {
                rep = new List<Message>();
                rep.AddRange(this.lst);
            }

            return rep;
        }

        internal void Add(Message msg)
        {
            this.lst.Add(msg);
            if(this.Level < msg.Level)
            {
                this.Level = msg.Level;
            }
        }

        internal void AddError(string msg, string src)
        {
            this.Add(msg, Message.Severity.Error, src);
        }

        internal void AddWarning(string msg, string src)
        {
            this.Add(msg, Message.Severity.Warning, src);
        }

        private void Add(string msg, Message.Severity lvl, string src)
        {
            this.Add(new Message(msg, lvl, src));
        }
    }
}
