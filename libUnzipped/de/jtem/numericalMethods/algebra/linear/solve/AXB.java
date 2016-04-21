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

package de.jtem.numericalMethods.algebra.linear.solve;

import de.jtem.numericalMethods.algebra.linear.decompose.Householder;

/**
 * Computes for a system of equations A�X=B with given matrices A and B the
 * solving or minimizing matrix X.
 * The algorithms used are the
 * {@link de.jtem.numericalMethods.algebra.linear.decompose.Householder Householder decomposition}
 * and the solver of {@link RXB R�X=B}.
 * <p />
 * If the left side matrix is square and regular then the system is solved
 * exactly and if not then it is minimized.
 * <p />
 * This algorithm is implemented for both real and complex valued equation
 * systems.
 * And each method is again available in two version.
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
 * The other version is an easy to use version that takes care of for example
 * not changing any arguments where not expected.
 * But this version can sometimes be very resource consuming!
 * <p />
 * Method versions for A�x=b:
 * <blockquote>
 * <table border="2" cellpadding="6">
 * <tr>
 * <th></th><th>optimized</th><th>easy to use</th>
 * </tr>
 * <tr>
 * <th>real valued</th>
 * <td>{@link #solve(double[][],double[],double[],double[]) solve(A,b,...)}</td>
 * <td>{@link #solve(double[][],double[],double[]) solve(A,x,b)}</td>
 * </tr>
 * <tr>
 * <th>complex valued</th>
 * <td>{@link #solve(double[][],double[][],double[],double[],double[],double[],double[],double[]) solve(A,b,...)}</td>
 * <td>{@link #solve(double[][],double[][],double[],double[],double[],double[]) solve(A,x,b)}</td>
 * </tr>
 * </table>
 * </blockquote>
 * Method versions for A�X=B:
 * <blockquote>
 * <table border="2" cellpadding="6">
 * <tr>
 * <th></th><th>optimized</th><th>easy to use</th>
 * </tr>
 * <tr>
 * <th>real valued</th>
 * <td>{@link #solve(double[][],double[][],double[],double[]) solve(A,B,...)}</td>
 * <td>{@link #solve(double[][],double[][],double[][]) solve(A,X,B)}</td>
 * </tr>
 * <tr>
 * <th>complex valued</th>
 * <td>{@link #solve(double[][],double[][],double[][],double[][],double[],double[],double[],double[]) solve(A,B,...)}</td>
 * <td>{@link #solve(double[][],double[][],double[][],double[][],double[][],double[][]) solve(A,X,B)}</td>
 * </tr>
 * </table>
 * </blockquote>
 * <em>
 * This class does not check any arguments for validity!<br />
 * For identical arguments results are undefined!
 * </em>
 *
 * @author  Samy Khadem-Al-Charieh
 * @see de.jtem.numericalMethods.algebra.linear.decompose.Householder Householder
 * @see RXB
 */
public class AXB
{
    /**
     * Minimizes a system of real valued equations.
     * <p />
     * @param   A   The left side and destroyed afterwards.
     * @param   b   The right side and the solution afterwards (for non square
     *              left sides the lower part is unchanged).
     * @param   t0  A double array of same length as the matrix A has columns
     *              for storing temporary results during the computation.
     * @param   t1  A double array of same length as the matrix A has columns
     *              for storing temporary results during the computation.
     * @return  True if system was solvable or minimizable.
     */
    public static boolean solve(double[][] A, double[] b, double[] t0, double[] t1)
    {
        int c = A[0].length;
        if (t0 == null)
            t0 = new double[c];
        if (Householder.decompose(A, t0, t1) == 0.)
            return false;
        Householder.qTimes(A, b, true);
        for (int i = 0; i < c; i++)
            A[i][i] = t0[i];
        return RXB.solve(A, b);
    }

    /**
     * Minimizes a system of real valued equations.
     * <p />
     * This method temporary allocates:
     * <ul>
     * <li>3 double arrays</li>
     * <li>1 dual indexed double array</li>
     * </ul>
     * <p />
     * @param   A   The left side.
     * @param   x   A vector taking the solution.
     * @param   b   The right side.
     * @return  True if system was solvable or minimizable.
     */
    public static boolean solve(double[][] A, double[] x, double[] b)
    {
        int r = A.length;
        int c = A[0].length;
        double[] t = new double[r];
        double[][] T = new double[r][c];
        System.arraycopy(b, 0, t, 0, r);
        for (int i = 0; i < r; i++)
            System.arraycopy(A[i], 0, T[i], 0, c);
        if (!solve(T, t, new double[c], new double[c]))
            return false;
        System.arraycopy(t, 0, x, 0, c);
        return true;
    }

    /**
     * Minimizes a system of complex valued equations.
     * <p />
     * @param   A_re    The real part of the left side and destroyed afterwards.
     * @param   A_im    The imaginary part of the left side and destroyed
     *                  afterwards.
     * @param   b_re    The real part of the right side and the solution
     *                  afterwards (for non square left sides the lower part
     *                  is unchanged).
     * @param   b_im    The imaginary part of the right side and the solution
     *                  afterwards (for non square left sides the lower part
     *                  is unchanged).
     * @param   t0      A double array of same length as the matrix A has
     *                  columns for storing temporary results during the
     *                  computation.
     * @param   t1      A double array of same length as the matrix A has
     *                  columns for storing temporary results during the
     *                  computation.
     * @param   t2      A double array of same length as the matrix A has
     *                  columns for storing temporary results during the
     *                  computation.
     * @param   t3      A double array of same length as the matrix A has
     *                  columns for storing temporary results during the
     *                  computation.
     * @return  True if system was solvable or minimizable.
     */
    public static boolean solve(double[][] A_re, double[][] A_im, double[] b_re, double[] b_im, double[] t0, double[] t1, double[] t2, double[] t3)
    {

        int c = A_re[0].length;
        if (t0 == null)
            t0 = new double[c];
        if (t1 == null)
            t1 = new double[c];
        double[] t = Householder.decompose(A_re, A_im, t0, t1, t2, t3);
        if (t[0] == 0. && t[1] == 0.)
            return false;
        Householder.qTimes(A_re, A_im, b_re, b_im, true);
        for (int i = 0; i < c; i++)
        {
            A_re[i][i] = t0[i];
            A_im[i][i] = t1[i];
        }
        return RXB.solve(A_re, A_im, b_re, b_im);
    }

    /**
     * Minimizes a system of complex valued equations.
     * <p />
     * This method temporary allocates:
     * <ul>
     * <li>6 double arrays</li>
     * <li>2 dual indexed double arrays</li>
     * </ul>
     * <p />
     * @param   A_re    The real part of the left side.
     * @param   A_im    The imaginary part of the left side.
     * @param   x_re    A vector taking the real part of the solution.
     * @param   x_im    A vector taking the imaginary part of the solution.
     * @param   b_re    The real part of the right side.
     * @param   b_im    The imaginary part of the right side.
     * @return  True if system was solvable or minimizable.
     */
    public static boolean solve(double[][] A_re, double[][] A_im, double[] x_re, double[] x_im, double[] b_re, double[] b_im)
    {
        int r = A_re.length;
        int c = A_re[0].length;
        double[] t_re = new double[r];
        double[] t_im = new double[r];
        double[][] T_re = new double[r][c];
        double[][] T_im = new double[r][c];
        System.arraycopy(b_re, 0, t_re, 0, r);
        System.arraycopy(b_im, 0, t_im, 0, r);
        for (int i = 0; i < r; i++)
        {
            System.arraycopy(A_re[i], 0, T_re[i], 0, c);
            System.arraycopy(A_im[i], 0, T_im[i], 0, c);
        }
        if (!solve(T_re, T_im, t_re, t_im, new double[c], new double[c], new double[c], new double[c]))
            return false;
        System.arraycopy(t_re, 0, x_re, 0, c);
        System.arraycopy(t_im, 0, x_im, 0, c);
        return true;
    }

    /**
     * Minimizes a system of real valued equations.
     * <p />
     * @param   A   The left side and destroyed afterwards.
     * @param   B   The right side and the solution afterwards (for non square
     *              left sides the lower part is unchanged).
     * @param   t0  A double array of same length as the matrix A has columns
     *              for storing temporary results during the computation.
     * @param   t1  A double array of same length as the matrix A has columns
     *              for storing temporary results during the computation.
     * @return  True if system was solvable or minimizable.
     */
    public static boolean solve(double[][] A, double[][] B, double[] t0, double[] t1)
    {
        int c = A[0].length;
        if (t0 == null)
            t0 = new double[c];
        if (Householder.decompose(A, t0, t1) == 0.)
            return false;
        Householder.qTimes(A, B, true, t1);
        for (int i = 0; i < A[0].length; i++)
            A[i][i] = t0[i];
        return RXB.solve(A, B);
    }

    /**
     * Minimizes a system of real valued equations.
     * <p />
     * This method temporary allocates:
     * <ul>
     * <li>2 double arrays</li>
     * <li>2 dual indexed double arrays</li>
     * </ul>
     * <p />
     * @param   A   The left side.
     * @param   X   A matrix taking the solution.
     * @param   B   The right side.
     * @return  True if system was solvable or minimizable.
     */
    public static boolean solve(double[][] A, double[][] X, double[][] B)
    {
        int r = A.length;
        int c = A[0].length;
        int n = B[0].length;
        double[][] T0 = new double[r][c];
        double[][] T1 = new double[r][n];
        for (int i = 0; i < r; i++)
        {
            System.arraycopy(A[i], 0, T0[i], 0, c);
            System.arraycopy(B[i], 0, T1[i], 0, n);
        }
        if (!solve(T0, T1, new double[c], new double[c]))
            return false;
        for (int i = 0; i < c; i++)
            System.arraycopy(T1[i], 0, X[i], 0, n);
        return true;
    }

    /**
     * Minimizes a system of complex valued equations.
     * <p />
     * @param   A_re    The real part of the left side and destroyed afterwards.
     * @param   A_im    The imaginary part of the left side and destroyed
     *                  afterwards.
     * @param   B_re    The real part of the right side and the solution
     *                  afterwards (for non square left sides the lower part
     *                  is unchanged).
     * @param   B_im    The imaginary part of the right side and the solution
     *                  afterwards (for non square left sides the lower part
     *                  is unchanged).
     * @param   t0      A double array of same length as the matrix A has
     *                  columns for storing temporary results during the
     *                  computation.
     * @param   t1      A double array of same length as the matrix A has
     *                  columns for storing temporary results during the
     *                  computation.
     * @param   t2      A double array of same length as the matrix A has
     *                  columns for storing temporary results during the
     *                  computation.
     * @param   t3      A double array of same length as the matrix A has
     *                  columns for storing temporary results during the
     *                  computation.
     * @return  True if system was solvable or minimizable.
     */
    public static boolean solve(double[][] A_re, double[][] A_im, double[][] B_re, double[][] B_im, double[] t0, double[] t1, double[] t2, double[] t3)
    {
        int c = A_re[0].length;
        if (t0 == null)
            t0 = new double[c];
        if (t1 == null)
            t1 = new double[c];
        double[] t = Householder.decompose(A_re, A_im, t0, t1, t2, t3);
        if (t[0] == 0. && t[1] == 0.)
            return false;
        Householder.qTimes(A_re, A_im, B_re, B_im, true, t2, t3);
        for (int i = 0; i < A_re[0].length; i++)
        {
            A_re[i][i] = t0[i];
            A_im[i][i] = t1[i];
        }
        return RXB.solve(A_re, A_im, B_re, B_im);
    }

    /**
     * Minimizes a system of complex valued equations.
     * <p />
     * This method temporary allocates:
     * <ul>
     * <li>4 double arrays</li>
     * <li>4 dual indexed double arrays</li>
     * </ul>
     * <p />
     * @param   A_re    The real part of the left side.
     * @param   A_im    The imaginary part of the left side.
     * @param   X_re    A matrix taking the real part of the solution.
     * @param   X_im    A matrix taking the imaginary part of the solution.
     * @param   B_re    The real part of the right side.
     * @param   B_im    The imaginary part of the right side.
     * @return  True if system was solvable or minimizable.
     */
    public static boolean solve(double[][] A_re, double[][] A_im, double[][] X_re, double[][] X_im, double[][] B_re, double[][] B_im)
    {
        int r = A_re.length;
        int c = A_re[0].length;
        int n = B_re[0].length;
        double[][] A0_re = new double[r][c];
        double[][] A0_im = new double[r][c];
        double[][] B0_re = new double[r][n];
        double[][] B0_im = new double[r][n];
        for (int i = 0; i < r; i++)
        {
            System.arraycopy(A_re[i], 0, A0_re[i], 0, c);
            System.arraycopy(A_im[i], 0, A0_im[i], 0, c);
            System.arraycopy(B_re[i], 0, B0_re[i], 0, n);
            System.arraycopy(B_im[i], 0, B0_im[i], 0, n);
        }
        if (!solve(A0_re, A0_im, B0_re, B0_im, new double[c], new double[c], new double[c], new double[c]))
            return false;
        for (int i = 0; i < c; i++)
        {
            System.arraycopy(B0_re[i], 0, X_re[i], 0, n);
            System.arraycopy(B0_im[i], 0, X_im[i], 0, n);
        }
        return true;
    }
}
