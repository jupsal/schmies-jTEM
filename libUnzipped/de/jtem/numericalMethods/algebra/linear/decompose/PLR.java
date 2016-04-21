/**
This file is part of a jTEM project.
All jTEM projects are licensed under the FreeBSD license 
or 2-clause BSD license (see http://www.opensource.org/licenses/bsd-license.php). 

Copyright (c) 2002-2009, Technische Universität Berlin, jTEM
All rights reserved.

Redistribution and use in source and binary forms, with or without modification, 
are permitted provided that the following conditions are met:

-	Redistributions of source code must retain the above copyright notice, 
	this list of conditions and the following disclaimer.

-	Redistributions in binary form must reproduce the above copyright notice, 
	this list of conditions and the following disclaimer in the documentation 
	and/or other materials provided with the distribution.
 
THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" 
AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, 
THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE 
ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS 
BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, 
OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT 
OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS 
INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, 
STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING 
IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY 
OF SUCH DAMAGE.
**/

package de.jtem.numericalMethods.algebra.linear.decompose;

/**
 * Computes for a square matrix A the matrices P, L and R complying P�A=L�R.
 * The algorithm used is the <em>Gauss algorithm with pivoting</em>.
 * <dl>
 * <dt>Advantages</dt>
 * <dd>� Fast</dd>
 * <dd>� Few additional resources</dd>
 * <dd>� Gives determinant</dd>
 * <dt>Disadvantages</dt>
 * <dd>� No control on condition of triangular matrices</dd>
 * </dl>
 * <p />
 * This algorithm is implemented for both real and complex valued equation
 * systems.
 * And each method is again available in two version.
 * <p />
 * One version is an optimized version that consumes as few resources as
 * possible.
 * <p />
 * The other version is an easy to use version that takes care of for example
 * not changing any arguments where not expected.
 * But this version can sometimes be very resource consuming!
 * <p />
 * Method versions for P�A=L�R:
 * <blockquote>
 * <table border="2" cellpadding="6">
 * <tr>
 * <th></th><th>optimized</th><th>easy to use</th>
 * </tr>
 * <tr>
 * <th>real valued</th>
 * <td>{@link #decompose(int[],double[][]) decompose(p,A)}</td>
 * <td>{@link #decompose(int[][],double[][],double[][],double[][]) decompose(P,A,L,R)}</td>
 * </tr>
 * <tr>
 * <th>complex valued</th>
 * <td>{@link #decompose(int[],double[][],double[][]) decompose(p,A)}</td>
 * <td>{@link #decompose(int[][],double[][],double[][],double[][],double[][],double[][],double[][]) decompose(P,A,L,R)}</td>
 * </tr>
 * </table>
 * </blockquote>
 * <em>
 * This class does not check any arguments for validity!<br />
 * For identical arguments results are undefined!
 * </em>
 *
 * @author  Samy Khadem-Al-Charieh
 * @see de.jtem.numericalMethods.algebra.linear.decompose decompose
 */
public class PLR
{
    /**
     * Decomposes a real matrix.
     * <p />
     * The left and right triangular matrices are stored merged into the
     * original matrix.
     * The diagonal belongs to the right triangular matrix because the diagonal
     * elements of the left triangular matrix are per definition equal to one.
     * For unmerging them see the methods in
     * {@link Unmerge#triangular(double[][],double[][]) Unmerge.triangular}.
     * <p />
     * The permutation is stored in an index vector and can be used like this:
     * <blockquote>
     * (P�A)<sub>ij</sub> &equiv; <code>A[p[i]][j]</code>
     * </blockquote>
     * <p />
     * @param   p   An index vector for taking the permutation.
     * @param   A   The matrix to decompose and afterwards the merged left and
     *              right triangular matrix.
     * @return  The determinant of the matrix.
     * @see Unmerge#triangular(double[][],double[][]) Unmerge.triangular
     */
    public static double decompose(int[] p, double[][] A)
    {
        int n = A.length;
        for (int i = 0; i < n; i++)
            p[i] = i;
        double det = 1.;
        for (int k = 0; k < n; k++)
        {
            double abs = Math.abs(A[k][k]);
            int piv = k;
            for (int i = k + 1; i < n; i++)
            {
                double t = Math.abs(A[i][k]);
                if (t > abs)
                {
                    piv = i;
                    abs = t;
                }
            }
            if (k != piv)
            {
                det *= -1.;
                int pk = p[k];
                p[k] = p[piv];
                p[piv] = pk;
                double[] Ak = A[k];
                A[k] = A[piv];
                A[piv] = Ak;
            }
            det *= A[k][k];
            for (int i = k + 1; i < n; i++)
            {
                if (abs != 0.)
                    A[i][k] /= A[k][k];
                else
                    A[i][k] = 0.;
                for (int j = k + 1; j < n; j++)
                    A[i][j] -= A[i][k] * A[k][j];
            }
        }
        return det;
    }

    /**
     * Decomposes a real matrix.
     * <p />
     * @param   P   A matrix taking the permutation matrix.
     * @param   A   The matrix to decompose.
     * @param   L   A matrix taking the left triangular matrix.
     * @param   R   A matrix taking the right triangular matrix.
     * @return  The determinant of the matrix.
     * @see Unmerge#triangular(double[][],double[][]) Unmerge.triangular
     */
    public static double decompose(int[][] P, double[][] A, double[][] L, double[][] R)
    {
        int n = A.length;
        for (int i = 0; i < n; i++)
            for (int j = 0; j < n; j++)
                L[i][j] = A[i][j];
        int[] p = new int[n];
        double det = decompose(p, L);
        Unmerge.triangular(L, R);
        for (int i = 0; i < n; i++)
        {
            for (int j = 0; j < n; j++)
                P[i][j] = 0;
            P[i][p[i]] = 1;
        }
        return det;
    }

    /**
     * Decomposes a complex matrix.
     * <p />
     * The left and right triangular matrices are stored merged into the
     * original matrix.
     * The diagonal belongs to the right triangular matrix because the diagonal
     * elements of the left triangular matrix are per definition equal to one.
     * For unmerging them see the methods in
     * {@link Unmerge#triangular(double[][],double[][],double[][],double[][]) Unmerge.triangular}.
     * <p />
     * The permutation is stored in an index vector and can be used like this:
     * <blockquote>
     * (P�A)<sub>ij</sub> &equiv; <code>A[p[i]][j]</code>
     * </blockquote>
     * <p />
     * @param   p       An index vector for taking the permutation.
     * @param   A_re    The real part of the matrix to decompose and afterwards
     *                  the real parts of the merged left and right triangular
     *                  matrix.
     * @param   A_im    The imaginary part of the matrix to decompose and
     *                  afterwards the imaginary parts of the merged left and
     *                  right triangular matrix.
     * @return  The determinant of the matrix.
     * @see Unmerge#triangular(double[][],double[][],double[][],double[][]) Unmerge.triangular
     */
    public static double[] decompose(int[] p, double[][] A_re, double[][] A_im)
    {
        int n = A_re.length;
        for (int i = 0; i < n; i++)
            p[i] = i;
        double[] det = new double[] {1., 0.};
        for (int k = 0; k < n; k++)
        {
            // find the pivot element
            int piv = k;
            double x_re = A_re[piv][k];
            double x_im = A_im[piv][k];
            double abs = x_re * x_re + x_im * x_im;
            for (int i = k + 1; i < n; i++)
            {
                double a_re = A_re[i][k];
                double a_im = A_im[i][k];
                double t = a_re * a_re + a_im * a_im;
                if (t > abs)
                {
                    piv = i;
                    abs = t;
                }
            }
            // exchange current row with pivot row
            if (k != piv)
            {
                det[0] *= -1;
                det[1] *= -1;
                int pk = p[k];
                p[k] = p[piv];
                p[piv] = pk;
                double[] Ak = A_re[k];
                A_re[k] = A_re[piv];
                A_re[piv] = Ak;
                Ak = A_im[k];
                A_im[k] = A_im[piv];
                A_im[piv] = Ak;
            }
            double d0 = det[0] * A_re[k][k] - det[1] * A_im[k][k];
            det[1] = det[0] * A_im[k][k] + det[1] * A_re[k][k];
            det[0] = d0;
            for (int i = k + 1; i < n; i++)
            {
                // Frobenius element: A(i,k)=A(i,k)/A(k,k)
                if (abs != 0.)
                {
                    double a_re = A_re[i][k];
                    double a_im = A_im[i][k];
                    double b_re = A_re[k][k];
                    double b_im = - A_im[k][k];
                    double t = b_re * b_re + b_im * b_im;
                    A_re[i][k] = (a_re * b_re - a_im * b_im) / t;
                    A_im[i][k] = (a_im * b_re + a_re * b_im) / t;
                }
                else
                {
                    A_re[i][k] = 0.;
                    A_im[i][k] = 0.;
                }
                // eliminate
                for (int j = k + 1; j < n; j++)
                {
                    // A(i,j)-=A(i,k)�A(k,j)
                    double a_re = A_re[i][k];
                    double a_im = A_im[i][k];
                    double b_re = A_re[k][j];
                    double b_im = A_im[k][j];
                    A_re[i][j] -= a_re * b_re - a_im * b_im;
                    A_im[i][j] -= a_im * b_re + a_re * b_im;
                }
            }
        }
        return det;
    }

    /**
     * Decomposes a complex matrix.
     * <p />
     * @param   P       A matrix taking the permutation matrix.
     * @param   A_re    The real part of the matrix to decompose.
     * @param   A_im    The imaginary part of the matrix to decompose.
     * @param   L_re    A matrix taking the real part of the left triangular
     *                  matrix.
     * @param   L_im    A matrix taking the imaginary part of the left
     *                  triangular matrix.
     * @param   R_re    A matrix taking the real part of the right triangular
     *                  matrix.
     * @param   R_im    A matrix taking the imaginary part of the right
     *                  triangular matrix.
     * @return  The determinant of the matrix.
     * @see Unmerge#triangular(double[][],double[][]) Unmerge.triangular
     */
    public static double[] decompose(int[][] P, double[][] A_re, double[][] A_im, double[][] L_re, double[][] L_im, double[][] R_re, double[][] R_im)
    {
        int n = A_re.length;
        for (int i = 0; i < n; i++)
            for (int j = 0; j < n; j++)
            {
                L_re[i][j] = A_re[i][j];
                L_im[i][j] = A_im[i][j];
            }
        int[] p = new int[n];
        double[] det = decompose(p, L_re, L_im);
        Unmerge.triangular(L_re, L_im, R_re, R_im);
        for (int i = 0; i < n; i++)
        {
            for (int j = 0; j < n; j++)
                P[i][j] = 0;
            P[i][p[i]] = 1;
        }
        return det;
    }
}
