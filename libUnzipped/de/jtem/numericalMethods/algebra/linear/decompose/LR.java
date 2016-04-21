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
 * Computes for a square matrix A the matrices L and R complying A=L�R.
 * The algorithm used is the <em>Gauss algorithm</em>.
 * <p />
 * If a zero happens to be on the diagonal while computing, this method
 * doesn't work.
 * So if the method is not successful that does <em>not</em> mean this
 * matrix is singular.
 * For a method having not this problem see {@link PLR}.
 * <dl>
 * <dt>Advantages</dt>
 * <dd>� Very fast</dd>
 * <dd>� No additional resources</dd>
 * <dd>� When successful gives determinant</dd>
 * <dt>Disadvantages</dt>
 * <dd>� No control on condition of triangular matrices</dd>
 * <dd>� Does not decompose all regular matrices</dd>
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
 * Method versions for A=L�R:
 * <blockquote>
 * <table border="2" cellpadding="6">
 * <tr>
 * <th></th><th>optimized</th><th>easy to use</th>
 * </tr>
 * <tr>
 * <th>real valued</th>
 * <td>{@link #decompose(double[][]) decompose(A)}</td>
 * <td>{@link #decompose(double[][],double[][],double[][]) decompose(A,L,R)}</td>
 * </tr>
 * <tr>
 * <th>complex valued</th>
 * <td>{@link #decompose(double[][],double[][]) decompose(A)}</td>
 * <td>{@link #decompose(double[][],double[][],double[][],double[][],double[][],double[][]) decompose(A,L,R)}</td>
 * </tr>
 * </table>
 * </blockquote>
 * <em>
 * This class does not check any arguments for validity!<br />
 * For identical arguments results are undefined!
 * </em>
 *
 * @author  Samy Khadem-Al-Charieh
 * @see PLR
 * @see de.jtem.numericalMethods.algebra.linear.decompose decompose
 */
public class LR
{
    /**
     * Decomposes a real matrix.
     * <p>
     * The left and right triangular matrices are stored merged into
     * the original matrix.
     * The diagonal belongs to the right triangular matrix because the diagonal
     * elements of the left triangular matrix are per definition equal to one.
     * For unmerging them see
     * {@link Unmerge#triangular(double[][],double[][]) Unmerge.triangular}.
     * </p>
     * @param   A   The matrix to decompose and afterwards the merged left and
     *              right triangular matrix.
     * @return  The determinant if the matrix could successfully be decomposed.
     *          A value of zero is returned for both either singular matrices
     *          or regular matrices that this algorithm didn't work on.
     * @see PLR
     * @see Unmerge#triangular(double[][],double[][]) Unmerge.triangular
     */
    public static double decompose(double[][] A)
    {
        int n = A.length;
        double det = 1.;
        for (int i = 0; i < n; i++)
        {
            for (int k = i; k < n; k++)
                for (int j = 0; j < i; j++)
                    A[i][k] -= A[i][j] * A[j][k];
            double d = A[i][i];
            det *= d;
            if (d != 0.)
                for (int k = i + 1; k < n; k++)
                {
                    double s = 0.;
                    for (int j = 0; j < i; j++)
                        s += A[k][j] * A[j][i];
                    A[k][i] = (A[k][i] - s) / d;
                }
            else
                for (int k = i + 1; k < n; k++)
                    A[k][i] = 0.;
        }
        return det;
    }

    /**
     * Decomposes a real matrix.
     * <p>
     * @param   A   The matrix to decompose.
     * @param   L   A matrix taking the left triangular matrix.
     * @param   R   A matrix taking the right triangular matrix.
     * @return  The determinant if the matrix could successfully be decomposed.
     *          A value of zero is returned for both either singular matrices
     *          or regular matrices that this algorithm didn't work on.
     * @see PLR
     */
    public static double decompose(double[][] A, double[][] L, double[][] R)
    {
        int n = A.length;
        for (int i = 0; i < n; i++)
            for (int j = 0; j < n; j++)
                L[i][j] = A[i][j];
        double det = decompose(L);
        Unmerge.triangular(L, R);
        return det;
    }

    /**
     * Decomposes a complex matrix.
     * <p>
     * The left and right triangular matrices are stored merged into the
     * original matrix.
     * The diagonal belongs to the right triangular matrix because the diagonal
     * elements of the left triangular matrix are per definition equal to one.
     * For unmerging them see
     * {@link Unmerge#triangular(double[][],double[][],double[][],double[][]) Unmerge.triangular}.
     * </p>
     * @param   A_re    The real part of the matrix to decompose and afterwards
     *                  the real parts of the merged left and right triangular
     *                  matrix.
     * @param   A_im    The imaginary part of the matrix to decompose and
     *                  afterwards the imaginary parts of the merged left and
     *                  right triangular matrix.
     * @return  The determinant if the matrix could successfully be decomposed.
     *          A value of zero is returned for both either singular matrices
     *          or regular matrices that this algorithm didn't work on.
     * @see PLR
     * @see Unmerge#triangular(double[][],double[][],double[][],double[][]) Unmerge.triangular
     */
    public static double[] decompose(double[][] A_re, double[][] A_im)
    {
        int n = A_re.length;
        double[] det = new double[] {1., 0.};
        for (int i = 0; i < n; i++)
        {
            // compute row in upper right triangle
            for (int k = i; k < n; k++)
            {
                double s_re = 0.;
                double s_im = 0.;
                for (int j = 0; j < i; j++)
                {
                    // s+=A(ij)�A(jk)
                    double a_re = A_re[i][j];
                    double a_im = A_im[i][j];
                    double b_re = A_re[j][k];
                    double b_im = A_im[j][k];
                    s_re += a_re * b_re - a_im * b_im;
                    s_im += a_im * b_re + a_re * b_im;
                }
                // A(ik)=A(ik)-s
                A_re[i][k] -= s_re;
                A_im[i][k] -= s_im;
            }
            // compute column in lower left triangle
            double d_re = A_re[i][i];
            double d_im = A_im[i][i];
            double d = det[0] * d_re - det[1] * d_im;
            det[1] = det[0] * d_im + det[1] * d_re;
            det[0] = d;
            d = d_re * d_re + d_im * d_im;
            if (d != 0.)
                for (int k = i + 1; k < n; k++)
                {
                    double s_re = 0.;
                    double s_im = 0.;
                    for (int j = 0; j < i; j++)
                    {
                        // s+=A(kj)�A(ji)
                        double a_re = A_re[k][j];
                        double a_im = A_im[k][j];
                        double b_re = A_re[j][i];
                        double b_im = A_im[j][i];
                        s_re += a_re * b_re - a_im * b_im;
                        s_im += a_im * b_re + a_re * b_im;
                    }
                    // A(ki)=(A(ki)-s)/d
                    s_re = A_re[k][i] - s_re;
                    s_im = A_im[k][i] - s_im;
                    A_re[k][i] = (s_re * d_re + s_im * d_im) / d;
                    A_im[k][i] = (s_im * d_re - s_re * d_im) / d;
                }
            else
                for (int k = i + 1; k < n; k++)
                {
                    A_re[k][i] = 0.;
                    A_im[k][i] = 0.;
                }
        }
        return det;
    }

    /**
     * Decomposes a complex matrix.
     * <p />
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
     * @return  The determinant if the matrix could successfully be decomposed.
     *          A value of zero is returned for both either singular matrices
     *          or regular matrices that this algorithm didn't work on.
     * @see PLR
     */
    public static double[] decompose(double[][] A_re, double[][] A_im, double[][] L_re, double[][] L_im, double[][] R_re, double[][] R_im)
    {
        int n = A_re.length;
        for (int i = 0; i < n; i++)
            for (int j = 0; j < n; j++)
            {
                L_re[i][j] = A_re[i][j];
                L_im[i][j] = A_im[i][j];
            }
        double[] det = decompose(L_re, L_im);
        Unmerge.triangular(L_re, L_im, R_re, R_im);
        return det;
    }
}
