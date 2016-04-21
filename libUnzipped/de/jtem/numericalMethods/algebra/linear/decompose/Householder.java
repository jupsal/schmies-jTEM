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
 * Computes for a matrix A the matrices Q and R complying A=Q�R.
 * The algorithm used is the <em>Householder algorithm</em>.
 * <p />
 * <dl>
 * <dt>Advantages</dt>
 * <dd>� Always good-natured</dd>
 * <dd>� Few additional resources</dd>
 * <dd>� Works on non-square matrices (more rows than columns)</dd>
 * <dd>� Gives determinant</dd>
 * <dt>Disadvantages</dt>
 * <dd>� Not very fast</dd>
 * </dl>
 * <p />
 * This algorithm is implemented for both real and complex valued equation
 * systems.
 * And each method is again available in three versions.
 * <p />
 * One version is an optimized version that consumes as few resources as
 * possible.
 * Moreover it is possible to pass objects that this method can store temporary
 * results during computation in them.
 * This is for example useful when calling this method multiple times so there
 * is no need for this method to allocate memory again and again.
 * If null is passed for the temporary objects then suitable objects are
 * allocated by the method itself.
 * <p />
 * The third version is an easy to use version that takes care of for example
 * not changing any arguments where not expected.
 * But this version can sometimes be very resource consuming!
 * <p />
 * The second version is kind of a mix between the first and the third version
 * giving the orthogonal and triangular matrix while destroying the original
 * matrix.
 * Additionally it is possible to reduce the number of columns of the
 * orthogonal matrix so that only some leading orthonormal vectors are given
 * back.
 * <p />
 * Method versions for decompose:
 * <blockquote>
 * <table border="2" cellpadding="6">
 * <tr>
 * <th></th><th>optimized</th><th>semi optimized</th><th>easy to use</th>
 * </tr>
 * <tr>
 * <th>real valued</th>
 * <td>{@link #decompose(double[][],double[],double[]) decompose(A,d,...)}</td>
 * <td>{@link #decompose(double[][],double[][],double[],double[]) decompose(A,Q,...)}</td>
 * <td>{@link #decompose(double[][],double[][],double[][]) decompose(A,Q,R)}</td>
 * </tr>
 * <tr>
 * <th>complex valued</th>
 * <td>{@link #decompose(double[][],double[][],double[],double[],double[],double[]) decompose(A,d,...)}</td>
 * <td>{@link #decompose(double[][],double[][],double[][],double[][],double[],double[],double[],double[]) decompose(A,Q,...)}</td>
 * <td>{@link #decompose(double[][],double[][],double[][],double[][],double[][],double[][]) decompose(A,Q,R)}</td>
 * </tr>
 * </table>
 * </blockquote>
 * <em>
 * This class does not check any arguments for validity!
 * <br />
 * For identical arguments results are undefined!
 * </em>
 *
 * @author  Samy Khadem-Al-Charieh
 * @see de.jtem.numericalMethods.algebra.linear.decompose decompose
 */
public class Householder
{
    /**
     * Decomposes a real matrix storing the result in a merged way.
     * <p />
     * After decomposition the matrix contains the merged right triangular
     * matrix and the producing vectors of the orthogonal matrix.
     * The diagonal belongs to the orthogonal matrix and the diagonal of the
     * right triangular matrix is stored in an extra vector.
     * <p />
     * For multiplying the orthogonal matrix Q stored by its producing vectors
     * with another matrix A the method
     * {@link #qTimes(double[][],double[][],boolean,double[]) qTimes}
     * can be used:
     * <blockquote>
     * Q�A &equiv; <code>qTimes(Q,A,false,null)</code>
     * <br />
     * Q<sup>T</sup>�A &equiv; <code>qTimes(Q,A,true,null)</code>
     * </blockquote>
     * <p />
     * @param   A   The matrix to decompose and afterwards the merged
     *              producing vectors of the orthogonal matrix and of the
     *              right triangular matrix.
     * @param   d   The diagonal of the right triangular matrix.
     * @param   t   A double array of same length as the matrix has columns
     *              for storing temporary results during the computation.
     * @return  a double with det(A)
     * @see #qTimes(double[][],double[][],boolean,double[]) qTimes
     */
    public static double decompose(double[][] A, double[] d, double[] t)
    {
        int r = A.length;
        int c = A[0].length;
        if (t == null)
            t = new double[c];
        double det = 1.;
        for (int k = 0; k < c && k < r; k++)
        {
            double sgn = A[k][k];
            double s = sgn * sgn;
            double akk = Math.sqrt(s);
            if (akk != 0.)
                sgn /= akk;
            else
                sgn = 1.;
            for (int i = k + 1; i < r; i++)
            {
                double s_ = A[i][k];
                s += s_ * s_;
            }
            s = Math.sqrt(s);
            double b = Math.sqrt(1. / (2. * s * (s + akk)));
            double a = sgn * s;
            d[k] = -a;
            det *= a;
            A[k][k] += a;
            for (int i = k; i < r; i++)
                A[i][k] *= b;
            for (int j = k + 1; j < c; j++)
            {
                double s_ = 0.;
                for (int i = k; i < r; i++)
                    s_ += A[i][k] * A[i][j];
                t[j] = 2. * s_;
            }
            for (int j = k + 1; j < c; j++)
                for (int i = k; i < r; i++)
                    A[i][j] -= A[i][k] * t[j];
        }
        return det;
    }

    /**
     * Decomposes a real matrix while separating the orthogonal matrix.
     * <p />
     * @param   A   The matrix to decompose and afterwards the right triangular
     *              matrix.
     * @param   Q   A matrix taking (a part of) the orthogonal matrix.
     *              If one doesn't want all orthonormal vectors it is possible
     *              to pass a matrix beeing not square.
     *              That is still having as much rows as A but having possibly
     *              less columns.
     * @param   t0  A double array of same length as the matrix has columns
     *              for storing temporary results during the computation.
     * @param   t1  A double array of same length as the matrix has rows
     *              for storing temporary results during the computation.
     * @return  a double with det(A)
     */
    public static double decompose(double[][] A, double[][] Q, double[] t0, double[] t1)
    {
        int r = A.length;
        int c = A[0].length;
        int n = Q[0].length;
        // call decompose
        double det = decompose(A, t0, t1);
        // reset Q to (part of) E
        for (int j = 0; j < n; j++)
            Q[0][j] = 0.;
        for (int i = 1; i < r; i++)
            System.arraycopy(Q[0], 0, Q[i], 0, n);
        for (int i = 0; i < n; i++)
            Q[i][i] = 1.;
        // compute Q�E
        qTimes(A, Q, false, t1);
        // reset lower triangle of A
        for (int j = 0; j < c; j++)
            A[c - 1][j] = 0.;
        for (int i = 1; i < c - 1; i++)
            System.arraycopy(A[c - 1], 0, A[i], 0, i);
        for (int i = c; i < r; i++)
            System.arraycopy(A[c - 1], 0, A[i], 0, c);
        // copy diagonal into A
        for (int i = 0; i < c; i++)
            A[i][i] = t0[i];
        return det;
    }

    /**
     * Decomposes a real matrix.
     * <p />
     * This method temporary allocates:
     * <ul>
     * <li>2 double arrays</li>
     * </ul>
     * <p />
     * @param   A   The matrix to decompose.
     * @param   Q   A matrix taking the orthogonal matrix.
     *              It must be square and must have the same number of rows as
     *              A has.
     * @param   R   A matrix taking the right triangular matrix.
     *              It must have the same number of rows and columns as A has.
     * @return  a double with det(A)
     */
    public static double decompose(double[][] A, double[][] Q, double[][] R)
    {
        int r = A.length;
        int c = A[0].length;
        // copy A into R
        for (int i = 0; i < r; i++)
            System.arraycopy(A[i], 0, R[i], 0, c);
        // call decompose
        double[] t0 = new double[c];
        double[] t1 = new double[r];
        double det = decompose(R, Q, t0, t1);
        return det;
    }

    /**
     * Decomposes a complex matrix storing the result in a merged way.
     * <p />
     * After decomposition the matrix contains the merged right triangular
     * matrix and the producing vectors of the unitary matrix.
     * The diagonal belongs to the unitary matrix and the diagonal of the
     * right triangular matrix is stored in an extra vector.
     * <p />
     * For multiplying the unitary matrix Q stored by its producing vectors
     * with another matrix A the method
     * {@link #qTimes(double[][],double[][],double[][],double[][],boolean,double[],double[]) qTimes}
     * can be used:
     * <blockquote>
     * Q�A &equiv; <code>qTimes(Q<sub>re</sub>,Q<sub>im</sub>,A<sub>re</sub>,A<sub>im</sub>,false,null,null)</code>
     * <br />
     * Q<sup>&lowast;</sup>�A &equiv; <code>qTimes(Q<sub>re</sub>,Q<sub>im</sub>,A<sub>re</sub>,A<sub>im</sub>,true,null,null)</code>
     * <br />
     * <i>(Note: Q<sup>&lowast;</sup> means transposed <b>and</b> conjugated.)</i>
     * </blockquote>
     *
     * @param   A_re    The real part of the matrix to decompose and afterwards
     *                  the real parts of the merged producing vectors of the
     *                  unitary matrix and of the right triangular matrix.
     * @param   A_im    The imaginary part of the matrix to decompose and
     *                  afterwards the imaginary parts of the merged producing
     *                  vectors of the unitary matrix and of the right
     *                  triangular matrix.
     * @param   d_re    The real part of the diagonal of the right triangular
     *                  matrix.
     * @param   d_im    The imaginary part of the diagonal of the right
     *                  triangular matrix.
     * @param   t0      A double array of same length as the matrix has columns
     *                  for storing temporary results during the computation.
     * @param   t1      A double array of same length as the matrix has columns
     *                  for storing temporary results during the computation.
     * @return  a double[] with det(A)
     * @see #qTimes(double[][],double[][],double[][],double[][],boolean,double[],double[]) qTimes
     */
    public static double[] decompose(double[][] A_re, double[][] A_im, double[] d_re, double[] d_im, double[] t0, double[] t1)
    {
        int r = A_re.length;
        int c = A_re[0].length;
        if (t0 == null)
            t0 = new double[c];
        if (t1 == null)
            t1 = new double[c];
        double det[] = new double[] {1., 0.};
        for (int k = 0; k < c && k < r; k++)
        {
            // a=|Akk|, sgn(Akk)=Akk/a
            double sgn_re = A_re[k][k];
            double sgn_im = A_im[k][k];
            double s = sgn_re * sgn_re + sgn_im * sgn_im;
            double a = Math.sqrt(s);
            if (a != 0.)
            {
                sgn_re /= a;
                sgn_im /= a;
            }
            else
            {
                sgn_re = 1.;
                sgn_im = 0.;
            }
            // s=||(Akk,...,Ark)||
            for (int i = k + 1; i < r; i++)
            {
                double s_re = A_re[i][k];
                double s_im = A_im[i][k];
                s += s_re * s_re + s_im * s_im;
            }
            s = Math.sqrt(s);
            // b�=1/||(Akk+sgn�s,...,Ark)||�=1/(2�s�(s+a))
            double b = Math.sqrt(1. / (2. * s * (s + a)));
            // Akk=-sgn�s (Akk is stored in d)
            double akk_re = sgn_re * s;
            double akk_im = sgn_im * s;
            d_re[k] = -akk_re;
            d_im[k] = -akk_im;
            a = det[0] * akk_re - det[1] * akk_im;
            det[1] = det[0] * akk_im + det[1] * akk_re;
            det[0] = a;
            // w=b�(Akk+sgn�s,...,Ark)
            // (w is stored in places of new zeroes in A)
            A_re[k][k] += akk_re;
            A_im[k][k] += akk_im;
            for (int i = k; i < r; i++)
            {
                A_re[i][k] *= b;
                A_im[i][k] *= b;
            }
            // Update A:
            // H�A=(I-2�w�w*)�A
            //    =A-w�2�w*�A
            // v*=2�w*�A
            for (int j = k + 1; j < c; j++)
            {
                double s_re = 0.;
                double s_im = 0.;
                for (int i = k; i < r; i++)
                {
                    double a_re = A_re[i][k];
                    double a_im = - A_im[i][k];
                    double b_re = A_re[i][j];
                    double b_im = A_im[i][j];
                    s_re += a_re * b_re - a_im * b_im;
                    s_im += a_im * b_re + a_re * b_im;
                }
                t0[j] = 2. * s_re;
                t1[j] = 2. * s_im;
            }
            // A-=w�v*
            for (int j = k + 1; j < c; j++)
                for (int i = k; i < r; i++)
                {
                    double a_re = A_re[i][k];
                    double a_im = A_im[i][k];
                    double b_re = t0[j];
                    double b_im = t1[j];
                    A_re[i][j] -= a_re * b_re - a_im * b_im;
                    A_im[i][j] -= a_im * b_re + a_re * b_im;
                }
        }
        return det;
    }

    /**
     * Decomposes a complex matrix while separating the orthogonal matrix.
     * <p />
     * @param   A_re    The real part of the matrix to decompose and afterwards
     *                  the real part of the right triangular matrix.
     * @param   A_im    The imaginary part of the matrix to decompose and
     *                  afterwards the imaginary part of the right triangular
     *                  matrix.
     * @param   Q_re    A matrix taking (a part of) the real part of the
     *                  unitary matrix.
     *                  If one doesn't want all orthonormal vectors it is
     *                  possible to pass a matrix beeing not square.
     *                  That is still having as much rows as A but having
     *                  possibly less columns.
     * @param   Q_im    A matrix taking (a part of) the imaginary part of the
     *                  unitary matrix.
     *                  If one doesn't want all orthonormal vectors it is
     *                  possible to pass a matrix beeing not square.
     *                  That is still having as much rows as A but having
     *                  possibly less columns.
     * @param   t0      A double array of same length as the matrix has columns
     *                  for storing temporary results during the computation.
     * @param   t1      A double array of same length as the matrix has columns
     *                  for storing temporary results during the computation.
     * @param   t2      A double array of same length as the matrix has rows
     *                  for storing temporary results during the computation.
     * @param   t3      A double array of same length as the matrix has rows
     *                  for storing temporary results during the computation.
     * @return  a double[] with det(A)
     */
    public static double[] decompose(double[][] A_re, double[][] A_im, double[][] Q_re, double[][] Q_im, double[] t0, double[] t1, double[] t2, double[] t3)
    {
        int r = A_re.length;
        int c = A_re[0].length;
        int n = Q_re[0].length;
        double[] det = decompose(A_re, A_im, t0, t1, t2, t3);
        for (int j = 0; j < n; j++)
            Q_re[0][j] = 0.;
        System.arraycopy(Q_re[0], 0, Q_im[0], 0, n);
        for (int i = 1; i < r; i++)
        {
            System.arraycopy(Q_re[0], 0, Q_re[i], 0, n);
            System.arraycopy(Q_im[0], 0, Q_im[i], 0, n);
        }
        for (int i = 0; i < n; i++)
            Q_re[i][i] = 1.;
        qTimes(A_re, A_im, Q_re, Q_im, false, t2, t3);
        for (int j = 0; j < c; j++)
            A_re[c - 1][j] = 0.;
        System.arraycopy(A_re[c - 1], 0, A_im[c - 1], 0, c);
        for (int i = 1; i < c - 1; i++)
        {
            System.arraycopy(A_re[c - 1], 0, A_re[i], 0, i);
            System.arraycopy(A_im[c - 1], 0, A_im[i], 0, i);
        }
        for (int i = c; i < r; i++)
        {
            System.arraycopy(A_re[c - 1], 0, A_re[i], 0, c);
            System.arraycopy(A_im[c - 1], 0, A_im[i], 0, c);
        }
        for (int i = 0; i < c; i++)
        {
            A_re[i][i] = t0[i];
            A_im[i][i] = t1[i];
        }
        return det;
    }

    /**
     * Decomposes a complex matrix.
     * <p />
     * This method temporary allocates:
     * <ul>
     * <li>4 double arrays</li>
     * </ul>
     * <p />
     * @param   A_re    The real part of the matrix to decompose.
     * @param   A_im    The imaginary part of the matrix to decompose.
     * @param   Q_re    The real part of a matrix taking the unitary matrix.
     *                  It must be square and must have the same number of rows
     *                  as A has.
     * @param   Q_im    The imaginary part of a matrix taking the unitary
     *                  matrix.
     *                  It must be square and must have the same number of rows
     *                  as A has.
     * @param   R_re    The real part of a matrix taking the right triangular
     *                  matrix.
     *                  It must have the same number of rows and columns as A
     *                  has.
     * @param   R_im    The imaginary part of a matrix taking the right
     *                  triangular matrix.
     *                  It must have the same number of rows and columns as A
     *                  has.
     * @return  a double[] with det(A)
     */
    public static double[] decompose(double[][] A_re, double[][] A_im, double[][] Q_re, double[][] Q_im, double[][] R_re, double[][] R_im)
    {
        int r = A_re.length;
        int c = A_re[0].length;
        for (int i = 0; i < r; i++)
        {
            System.arraycopy(A_re[i], 0, R_re[i], 0, c);
            System.arraycopy(A_im[i], 0, R_im[i], 0, c);
        }
        double[] t0_re = new double[c];
        double[] t0_im = new double[c];
        double[] t1_re = new double[r];
        double[] t1_im = new double[r];
        double[] det = decompose(R_re, R_im, Q_re, Q_im, t0_re, t0_im, t1_re, t1_im);
        return det;
    }

    /**
     * Computes the product of a matrix with the producing vectors of an
     * orthogonal matrix and a vector.
     * <p />
     * @param   Q       The matrix with the producing vectors of the orthogonal
     *                  matrix.
     * @param   a       A vector to be multiplied with and afterwards the
     *                  result.
     * @param   conj    If false compute Q�a and if true Q<sup>T</sup>�a.
     * @see #decompose(double[][],double[],double[]) decompose
     */
    public static void qTimes(double[][] Q, double[] a, boolean conj)
    {
        int r = Q.length;
        int c = Q[0].length;
        // (I-2�w�w*)�a
        for (int l = 0; l < c; l++)
        {
            int k = conj ? l : (c - 1) - l;
            // t=2�w*�a
            double t = 0.;
            for (int i = k; i < r; i++)
                t += Q[i][k] * a[i];
            t *= 2.;
            // a-=w�t
            for (int i = k; i < r; i++)
                a[i] -= Q[i][k] * t;
        }
    }

    /**
     * Computes the product of a matrix with the producing vectors of an
     * unitary matrix with a vector.
     * <p />
     * @param   Q_re    The real part of the matrix with the producing vectors
     *                  of the unitary matrix.
     * @param   Q_im    The imaginary part of the matrix with the producing
     *                  vectors of the unitary matrix.
     * @param   a_re    The real part of a vector to be multiplied with and
     *                  afterwards the result.
     * @param   a_im    The imaginary part of a vector to be multiplied with
     *                  and afterwards the result.
     * @param   conj    If false compute Q�a and if true Q<sup>&lowast;</sup>�a.
     *                  <br />
     *                  <i>(Note: Q<sup>&lowast;</sup> means transposed <b>and</b> conjugated.)</i>
     * @see #decompose(double[][],double[][],double[],double[],double[],double[]) decompose
     */
    public static void qTimes(double[][] Q_re, double[][] Q_im, double[] a_re, double[] a_im, boolean conj)
    {
        int r = Q_re.length;
        int c = Q_re[0].length;
        // (I-2�w�w*)�a
        for (int l = 0; l < c; l++)
        {
            int k = conj ? l : (c - 1) - l;
            // t=2�w*�a
            double t_re = 0.;
            double t_im = 0.;
            for (int i = k; i < r; i++)
            {
                double a__re =   Q_re[i][k];
                double a__im = - Q_im[i][k];
                double b_re = a_re[i];
                double b_im = a_im[i];
                t_re += a__re * b_re - a__im * b_im;
                t_im += a__im * b_re + a__re * b_im;
            }
            t_re *= 2.;
            t_im *= 2.;
            // a-=w�t
            for (int i = k; i < r; i++)
            {
                double a__re = Q_re[i][k];
                double a__im = Q_im[i][k];
                a_re[i] -= a__re * t_re - a__im * t_im;
                a_im[i] -= a__im * t_re + a__re * t_im;
            }
        }
    }

    /**
     * Computes the product of a matrix with the producing vectors of an
     * orthogonal matrix and a matrix.
     * <p />
     * @param   Q       The matrix with the producing vectors of the orthogonal
     *                  matrix.
     * @param   A       A matrix to be multiplied with and afterwards the
     *                  result.
     * @param   conj    If false compute Q�A and if true Q<sup>T</sup>�A.
     * @param   t       A double array of same length as the matrix A has
     *                  columns for storing temporary results during the
     *                  computation.
     * @see #decompose(double[][],double[],double[]) decompose
     */
    public static void qTimes(double[][] Q, double[][] A, boolean conj, double[] t)
    {
        int r = Q.length;
        int c = Q[0].length;
        int n = A[0].length;
        if (t == null)
            t = new double[n];
        // (I-2�w�w*)�A
        for (int l = 0; l < c; l++)
        {
            int k = conj ? l : (c - 1) - l;
            // t*=2�w*�A
            for (int j = 0; j < n; j++)
            {
                double s = 0.;
                for (int i = k; i < r; i++)
                    s += Q[i][k] * A[i][j];
                t[j] = 2. * s;
            }
            // A-=w�t*
            for (int i = k; i < r; i++)
                for (int j = 0; j < n; j++)
                    A[i][j] -= Q[i][k] * t[j];
        }
    }

    /**
     * Computes the product of a matrix with the producing vectors of an
     * unitary matrix with a matrix.
     * <p />
     * @param   Q_re    The real part of the matrix with the producing vectors
     *                  of the unitary matrix.
     * @param   Q_im    The imaginary part of the matrix with the producing
     *                  vectors of the unitary matrix.
     * @param   A_re    The real part of a matrix to be multiplied with and
     *                  afterwards the result.
     * @param   A_im    The imaginary part of a matrix to be multiplied with
     *                  and afterwards the result.
     * @param   conj    If false compute Q�A and if true Q<sup>&lowast;</sup>�A.
     *                  <br />
     *                  <i>(Note: Q<sup>&lowast;</sup> means transposed <b>and</b> conjugated.)</i>
     * @param   t0      A double array of same length as the matrix A has
     *                  columns for storing temporary results during the
     *                  computation.
     * @param   t1      A double array of same length as the matrix A has
     *                  columns for storing temporary results during the
     *                  computation.
     * @see #decompose(double[][],double[][],double[],double[],double[],double[]) decompose
     */
    public static void qTimes(double[][] Q_re, double[][] Q_im, double[][] A_re, double[][] A_im, boolean conj, double[] t0, double[] t1)
    {
        int r = Q_re.length;
        int c = Q_re[0].length;
        int n = A_re[0].length;
        if (t0 == null)
            t0 = new double[n];
        if (t1 == null)
            t1 = new double[n];
        // (I-2�w�w*)�A
        for (int l = 0; l < c; l++)
        {
            int k = conj ? l : (c - 1) - l;
            // t*=2�w*�A
            for (int j = 0; j < n; j++)
            {
                double s_re = 0.;
                double s_im = 0.;
                for (int i = k; i < r; i++)
                {
                    double a_re = Q_re[i][k];
                    double a_im = - Q_im[i][k];
                    double b_re = A_re[i][j];
                    double b_im = A_im[i][j];
                    s_re += a_re * b_re - a_im * b_im;
                    s_im += a_im * b_re + a_re * b_im;
                }
                t0[j] = 2. * s_re;
                t1[j] = 2. * s_im;
            }
            // A-=w�t*
            for (int i = k; i < r; i++)
                for (int j = 0; j < n; j++)
                {
                    double a_re = Q_re[i][k];
                    double a_im = Q_im[i][k];
                    double b_re = t0[j];
                    double b_im = t1[j];
                    A_re[i][j] -= a_re * b_re - a_im * b_im;
                    A_im[i][j] -= a_im * b_re + a_re * b_im;
                }
        }
    }
}
