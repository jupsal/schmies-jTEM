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
 * Computes for a square matrix A the matrices Q and R complying A=Q�R.
 * The algorithm used is equivalent to <em>Gram Schmidt orthonormalising</em>.
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
 * Method versions for A=Q�R:
 * <blockquote>
 * <table border="2" cellpadding="6">
 * <tr>
 * <th></th><th>optimized</th><th>easy to use</th>
 * </tr>
 * <tr>
 * <th>real valued</th>
 * <td>{@link #decompose(double[][],double[][]) decompose(Q,R)}</td>
 * <td>{@link #decompose(double[][],double[][],double[][]) decompose(A,Q,R)}</td>
 * </tr>
 * <tr>
 * <th>complex valued</th>
 * <td>{@link #decompose(double[][],double[][],double[][],double[][]) decompose(Q,R)}</td>
 * <td>{@link #decompose(double[][],double[][],double[][],double[][],double[][],double[][]) decompose(A,Q,R)}</td>
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
public class QR
{
    /**
     * Decomposes a real matrix.
     * <p />
     * @param   Q   The matrix to decompose and afterwards the orthogonal
     *              matrix.
     * @param   R   A matrix taking the right triangular matrix.
     * @return  ||det(A)||
     */
    public static double decompose(double[][] Q, double[][] R)
    {
        int n = Q.length;
        double det = 1.;
        for (int j = 0; j < n; j++)
        {
            for (int k = 0; k < j; k++)
            {
                R[j][k] = 0.;
                double t = 0.;
                for (int i = 0; i < n; i++)
                    t += Q[i][k] * Q[i][j];
                R[k][j] = t;
                for (int i = 0; i < n; i++)
                    Q[i][j] -= t * Q[i][k];
            }
            double t = 0.;
            for (int i = 0; i < n; i++)
                t += Q[i][j] * Q[i][j];
            t = Math.sqrt(t);
            det *= t;
            R[j][j] = t;
            if (t != 0.)
                for (int i = 0; i < n; i++)
                    Q[i][j] /= t;
            else
                for (int i = 0; i < n; i++)
                    Q[i][j] = 0.;
        }
        return det;
    }

    /**
     * Decomposes a real matrix.
     * <p />
     * @param   A   The matrix to decompose.
     * @param   Q   A matrix taking the orthogonal matrix.
     * @param   R   A matrix taking the right triangular matrix.
     * @return  ||det(A)||
     */
    public static double decompose(double[][] A, double[][] Q, double[][] R)
    {
        int n = A.length;
        for (int i = 0; i < n; i++)
            for (int j = 0; j < n; j++)
                Q[i][j] = A[i][j];
        return decompose(Q, R);
    }

    /**
     * Decomposes a complex matrix.
     * <p />
     * @param   Q_re    The real part of the matrix to decompose and afterwards
     *                  the real part of the unitary matrix.
     * @param   Q_im    The imaginary part of the matrix to decompose and
     *                  afterwards the imaginary part of the unitary matrix.
     * @param   R_re    A matrix taking the real part of the right triangular
     *                  matrix.
     * @param   R_im    A matrix taking the imaginary part of the right
     *                  triangular matrix.
     * @return  ||det(A)||
     */
    public static double[] decompose(double[][] Q_re, double[][] Q_im, double[][] R_re, double[][] R_im)
    {
        // A=(a(0)|..|a(n-1)), Q=(q(0)|..|q(n-1)), R=(r(ij))
        // A=Q�R <=> R=Q*�A
        // a(j)=sum[i=0..j]q(i)�r(ij)=sum[i=0..j-1]q(i)�r(ij)+q(j)�r(jj)
        // =>
        // q(j)=(a(j)-sum[i=0..j-1]q(i)�r(ij))/r(jj)
        // r(ij)=q(i)*�a(j)
        // r(jj)=||a(j)-sum[i=0..j-1]q(i)�r(ij)||
        int n = Q_re.length;
        double det = 1.;
        for (int j = 0; j < n; j++)
        {
            // q(j)'=a(j)-sum[i=0..j-1](q(i)�r(ij))
            for (int i = 0; i < j; i++)
            {
                // r(ji)=0
                R_re[j][i] = 0.;
                R_im[j][i] = 0.;
                // r(ij)=q(i)*�a(j)
                double t_re = 0.;
                double t_im = 0.;
                for (int k = 0; k < n; k++)
                {
                    double a_re = Q_re[k][i];
                    double a_im = - Q_im[k][i];
                    double b_re = Q_re[k][j];
                    double b_im = Q_im[k][j];
                    t_re += a_re * b_re - a_im * b_im;
                    t_im += a_im * b_re + a_re * b_im;
                }
                R_re[i][j] = t_re;
                R_im[i][j] = t_im;
                // q(j)'-=q(i)�r(ij)
                for (int k = 0; k < n; k++)
                {
                    double b_re = Q_re[k][i];
                    double b_im = Q_im[k][i];
                    Q_re[k][j] -= t_re * b_re - t_im * b_im;
                    Q_im[k][j] -= t_im * b_re + t_re * b_im;
                }
            }
            // r(jj)=||q(j)'||
            double t = 0.;
            for (int i = 0; i < n; i++)
                t += Q_re[i][j] * Q_re[i][j] + Q_im[i][j] * Q_im[i][j];
            t = Math.sqrt(t);
            det *= t;
            R_re[j][j] = t;
            R_im[j][j] = 0.;
            // q(j)=q(j)'/r(jj)
            if (t != 0.)
                for (int i = 0; i < n; i++)
                {
                    Q_re[i][j] /= t;
                    Q_im[i][j] /= t;
                }
            else
                for (int i = 0; i < n; i++)
                {
                    Q_re[i][j] = 0.;
                    Q_im[i][j] = 0.;
                }
        }
        return new double[] {det, 0.};
    }

    /**
     * Decomposes a complex matrix.
     * <p />
     * @param   A_re    The real part of the matrix to decompose.
     * @param   A_im    The imaginary part of the matrix to decompose.
     * @param   Q_re    A matrix taking the real part of the unitary matrix.
     * @param   Q_im    A matrix taking the imaginary part of the unitary
     *                  matrix.
     * @param   R_re    A matrix taking the real part of the right triangular
     *                  matrix.
     * @param   R_im    A matrix taking the imaginary part of the right
     *                  triangular matrix.
     * @return  ||det(A)||
     */
    public static double[] decompose(double[][] A_re, double[][] A_im, double[][] Q_re, double[][] Q_im, double[][] R_re, double[][] R_im)
    {
        int n = A_re.length;
        for (int i = 0; i < n; i++)
            for (int j = 0; j < n; j++)
            {
                Q_re[i][j] = A_re[i][j];
                Q_im[i][j] = A_im[i][j];
            }
        return decompose(Q_re, Q_im, R_re, R_im);
    }
}
