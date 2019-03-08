
namespace Zygotine.Numerics
{
    using System;
    using System.Text;
    using Zygotine.Util;

    public class XInverse4Y
    {
        private int nRow; // number of rows 
        private int nCol; // number of columns including Y
        private double[,] z = null;
        private XInverse4Y(double[] byCol, double[] y)
        {
            nRow = y.Length;
            if ((nRow < 2) || (byCol.Length != nRow * nRow))
            {
                //TODO exception à formaliser
                throw new Exception();
            }

            nCol = nRow + 1;
            z = new double[nRow, nCol];
            for (int i = 0; i < byCol.Length; i++)
            {
                int icol = i / nRow;
                int ilig = i % nRow;
                z[ilig, icol] = byCol[i];
            }

            for (int i = 0; i < y.Length; i++)
            {
                z[i, nCol - 1] = y[i];
            }

        }
        private double[] Diagonal()
        {
            double[] rep = new double[nRow];
            for (int i = 0; i < nRow; i++)
            {
                rep[i] = z[i, i];
            }

            return rep;
        }

        private double[] Compute()
        {

            for (int iR = 0; iR < nRow - 1; iR++)
            {
                for (int iS = iR + 1; iS < nRow; iS++)
                {
                    //z[s,] - z[s,r]/z[r,r]*z[r,]
                    double tmp = -z[iS, iR] / z[iR, iR];
                    double[] row = GetRowAt(iS).Add(GetRowAt(iR).Multiply(tmp));
                    SetRowAt(iS, row);
                    //                       z[s,] < -z[s,] - z[s, r] / z[r, r] * z[r,]
                }
            }

            for (int iR = 0; iR < nRow - 1; iR++)
            {
                for (int iS = iR + 1; iS < nRow; iS++)
                {
                    double tmp = -z[iR, iS] / z[iS, iS];
                    double[] row = GetRowAt(iR).Add(GetRowAt(iS).Multiply(tmp));
                    SetRowAt(iR, row);
                    //z[r,] <- z[r,] - z[r,s]/z[s,s]*z[s,]
                }
            }

            double[] y = this.GetColumnAt(nCol - 1);
            double[] diagonal = this.Diagonal();
            double[] logAbsTheta = y.Abs().Log().Substract(diagonal.Abs().Log());
            double[] theta = logAbsTheta.Exp();
            for (int i = 0; i < nRow; i++)
            {
                theta[i] = theta[i] * (Math.Sign(y[i]) * Math.Sign(diagonal[i]));
            }
            return theta;
            //9.0909091 -3.0909091 -1.3333333 -0.6666667
        }

        private string Format(double[] tmp)
        {
            StringBuilder sb = new StringBuilder();
            for (int i = 0; i < tmp.Length; i++)
            {
                sb.AppendFormat("{0,16}", tmp[i]);
            }
            sb.AppendLine();

            return sb.ToString();
        }

        private double[] GetRowAt(int iRow)
        {
            double[] row = new double[nCol];
            for (int iC = 0; iC < nCol; iC++)
            {
                row[iC] = z[iRow, iC];
            }
            return row;
        }

        public override string ToString()
        {
            StringBuilder sb = new StringBuilder();
            for (int iR = 0; iR < nRow; iR++)
            {
                for (int iC = 0; iC < nCol; iC++)
                {
                    sb.AppendFormat("{0,16}", z[iR, iC]);
                }
                sb.AppendLine();
            }

            return sb.ToString();
        }

        private double[,] SetRowAt(int iRow, double[] row)
        {
            for (int iC = 0; iC < nCol; iC++)
            {
                z[iRow, iC] = row[iC];
            }
            return this.z;
        }

        private double[] GetColumnAt(int iCol)
        {
            double[] col = new double[nRow];
            for (int iRow = 0; iRow < nRow; iRow++)
            {
                col[iRow] = z[iRow, iCol];
            }
            return col;
        }

        public static double[] Compute( double [] x, double[] y)
        {
            return new XInverse4Y(x, y).Compute();
        }

        /*
         * Extrapolation m = new Extrapolation(new double[] { 1.1, 2, 3, 4, 0, 2, 3, 4, 0, 0, 3, 4, 0, 0, 0, 4 }, new double[] { 10, 12, 14, 16 });
         * m.Compute();
         */
    }


}
