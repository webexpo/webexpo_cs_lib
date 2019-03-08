namespace Zygotine.Numerics
{
    using System;
    using System.Text;

    public class Matrix
    {
        private double[][] matrix;

        public int NCols { get; private set; }
        public int NRows { get; private set; }

        public Matrix(int nRows, int nCols)
        {
            this.NRows = nRows;
            this.NCols = nCols;
            this.matrix = new double[NRows][];
            for (int iRow = 0; iRow < NRows; iRow++)
            {
                this.matrix[iRow] = new double[NCols];
            }
        }

        public Matrix(int n) : this(n, n)
        {
        }

        /*
         * data.Length >= nRows x nCols
         */
        public Matrix(int nRows, int nCols, double[] data, bool byCol = true) : this(nRows, nCols)
        {
            int iData = 0;
            if (byCol)
            {
                for (int iCol = 0; iCol < this.NCols; iCol++)
                {
                    for (int iRow = 0; iRow < this.NRows; iRow++)
                    {
                        this.matrix[iRow][iCol] = data[iData];
                        iData++;
                    }
                }
            }
            else
            {
                for (int iRow = 0; iRow < this.NRows; iRow++)
                {
                    for (int iCol = 0; iCol < this.NCols; iCol++)
                    {
                        this.matrix[iRow][iCol] = data[iData];
                        iData++;
                    }
                }

            }
        }

        public Matrix(int nRows, int nCols, double val, bool byCol = true) : this(nRows, nCols)
        {
            for (int iCol = 0; iCol < this.NCols; iCol++)
            {
                for (int iRow = 0; iRow < this.NRows; iRow++)
                {
                    this.matrix[iRow][iCol] = val;
                }
            }
        }

        public double[] ToArray(bool byCol = true)
        {
            double[] data = new double[NRows * NCols];
            int iData = 0;
            if (byCol)
            {
                for (int iCol = 0; iCol < this.NCols; iCol++)
                {
                    for (int iRow = 0; iRow < this.NRows; iRow++)
                    {
                        data[iData] = this.matrix[iRow][iCol];
                        iData++;
                    }
                }
            }
            else
            {
                for (int iRow = 0; iRow < this.NRows; iRow++)
                {
                    for (int iCol = 0; iCol < this.NCols; iCol++)
                    {
                        data[iData] = this.matrix[iRow][iCol];
                        iData++;
                    }
                }
            }

            return data;
        }

        // Access this matrix as a 2D array
        public double this[int iRow, int iCol]
        {
            get { return this.matrix[iRow][iCol]; }
            set { this.matrix[iRow][iCol] = value; }
        }

        public void MultiplyBy(double x)
        {
            for (int iRow = 0; iRow < this.NRows; iRow++)
            {
                for (int iCol = 0; iCol < this.NCols; iCol++)
                {
                    this.matrix[iRow][iCol] *= x;
                }
            }
        }

        public string AsString()
        {
            StringBuilder sb = new StringBuilder();
            for (int i = 0; i < matrix.Length; ++i)
            {
                for (int j = 0; j < matrix[i].Length; ++j)
                {
                    sb.AppendFormat("{0,15:E} ", this.matrix[i][j]);
                }

                sb.AppendLine();
            }
            return sb.ToString();
        }

        public static Matrix CreateColumn(int n)
        {
            return new Matrix(1, n);
        }

        public static Matrix CreateRow(int n)
        {
            return new Matrix(n, 1);
        }


        public static Matrix Add(Matrix m1, Matrix m2)
        {
            if ((m1.NRows != m2.NRows) || (m1.NCols != m2.NCols))
            {
                //TODO Gérer l'exception
                throw new Exception("invalid Matrix dimensions");
            }

            Matrix rep = new Matrix(m1.NRows, m1.NCols);
            for (int iRow = 0; iRow < m1.NRows; iRow++)
            {
                for (int iCol = 0; iCol < m1.NCols; iCol++)
                {
                    rep.matrix[iRow][iCol] = m1.matrix[iRow][iCol] + m2.matrix[iRow][iCol];
                }
            }

            return rep;
        }


        public static Matrix Product(Matrix m1, Matrix m2)
        {
            if (m1.NCols != m2.NRows)
            {
                //TODO Gérer l'exception
                throw new Exception("invalid Matrix dimensions");
            }

            Matrix rep = new Matrix(m1.NRows, m2.NCols);
            for (int iRow = 0; iRow < m1.NRows; iRow++)
            {
                for (int iCol = 0; iCol < m2.NCols; iCol++)
                {
                    rep[iRow, iCol] = 0;
                    for (int i = 0; i < m1.NCols; i++)
                    {
                        rep[iRow, iCol] += m1[iRow, i] * m2[i, iCol];
                    }
                }
            }

            return rep;
        }
    }
}
