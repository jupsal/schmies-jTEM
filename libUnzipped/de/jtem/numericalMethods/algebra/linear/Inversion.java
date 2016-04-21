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

package de.jtem.numericalMethods.algebra.linear;

import de.jtem.numericalMethods.algebra.linear.solve.AXB;

/**
 * Computes the inverse of a square matrix.
 * <p>
 * This algorithm is implemented for both real and complex valued equation
 * systems.
 * And each method is again available in two version.
 * <p>
 * One version is an optimized version that consumes as few resources as
 * possible.
 * Moreover it is possible to pass objects that this method can store temporary
 * results during computation in them.
 * This is for example useful when calling this method multiple times so there
 * is no need for this method to allocate memory again and again.
 * If null is passed for the temporary objects then suitable objects are
 * allocated by the method itself.
 * <p>
 * The other version is an easy to use version that takes care of for example
 * not changing any arguments where not expected.
 * But this version can sometimes be very resource consuming!
 * <p>
 * Method versions for Inversion:
 * <blockquote>
 * <table border="2" cellpadding="6">
 * <tr>
 * <th></th><th>optimized</th><th>easy to use</th>
 * </tr>
 * <tr>
 * <th>real valued</th>
 * <td>{@link #compute(double[][],double[][],double[],double[]) compute(A,...)}</td>
 * <td>{@link #compute(double[][],double[][]) compute(A,I)}</td>
 * </tr>
 * <tr>
 * <th>complex valued</th>
 * <td>{@link #compute(double[][],double[][],double[][],double[][],double[],double[],double[],double[]) compute(A,...)}</td>
 * <td>{@link #compute(double[][],double[][],double[][],double[][]) compute(A,I)}</td>
 * </tr>
 * </table>
 * </blockquote>
 * <p>
 * <em>
 * This class does not check any arguments for validity!<br />
 * For identical arguments results are undefined!
 * </em>
 *
 * @author  Samy Khadem-Al-Charieh
 */
public class Inversion
{
    /**
     * Inverts a real valued square matrix.
     * <p>
     * @param   A   The square matrix to invert and its inverse afterwards.
     * @param   T   A dual indexed double array of same size as matrix A for
     *              storing temporary results during the computation.
     * @param   t0  A double array of same size as matrix A has columns for
     *              storing temporary results during the computation.
     * @param   t1  A double array of same size as matrix A has columns for
     *              storing temporary results during the computation.
     * @return  True if matrix was regular and false if it was singular.
     */
    public static boolean compute(double[][] A, double[][] T, double[] t0, double[] t1)
    {
        int n = A[0].length;
        if (T == null)
            T = new double[n][n];
        for (int j = 0; j < n; j++)
            T[0][j] = 0.;
        for (int i = 1; i < n; i++)
            System.arraycopy(T[0], 0, T[i], 0, n);
        for (int i = 0; i < n; i++)
            T[i][i] = 1.;
        if (!AXB.solve(A, T, t0, t1))
            return false;
        for (int i = 0; i < n; i++)
            System.arraycopy(T[i], 0, A[i], 0, n);
        return true;
    }

    /**
     * Inverts a real valued square matrix.
     * <p>
     * This method temporary allocates:
     * <ul>
     * <li>2 double arrays</li>
     * <li>1 dual indexed double array</li>
     * </ul>
     * <p>
     * @param   A   The square matrix.
     * @param   I   A matrix taking the inverse.
     * @return  True if matrix was regular and false if it was singular.
     */
    public static boolean compute(double[][] A, double[][] I)
    {
        int n = A[0].length;
        if( A != I ) {
        	for (int i = 0; i < n; i++)
        		System.arraycopy(A[i], 0, I[i], 0, n);
        }
        return compute(I, (double[][]) null, (double[]) null, (double[]) null);
    }

    /**
     * Inverts a complex valued square matrix.
     * <p>
     * @param   A_re    The real part of the square matrix to invert and its
     *                  inverse afterwards.
     * @param   A_im    The imaginary part of the square matrix to invert and
     *                  its inverse afterwards.
     * @param   T0      A dual indexed double array of same size as matrix A
     *                  for storing temporary results during the computation.
     * @param   T1      A dual indexed double array of same size as matrix A
     *                  for storing temporary results during the computation.
     * @param   t0      A double array of same size as matrix A has columns for
     *                  storing temporary results during the computation.
     * @param   t1      A double array of same size as matrix A has columns for
     *                  storing temporary results during the computation.
     * @param   t2      A double array of same size as matrix A has columns for
     *                  storing temporary results during the computation.
     * @param   t3      A double array of same size as matrix A has columns for
     *                  storing temporary results during the computation.
     * @return  True if matrix was regular and false if it was singular.
     */
    public static boolean compute(double[][] A_re, double[][] A_im, double[][] T0, double[][] T1, double[] t0, double[] t1, double[] t2, double[] t3)
    {
        int n = A_re[0].length;
        if (T0 == null)
            T0 = new double[n][n];
        if (T1 == null)
            T1 = new double[n][n];
        for (int j = 0; j < n; j++)
            T1[0][j] = 0.;
        System.arraycopy(T1[0], 0, T0[0], 0, n);
        T0[0][0] = 1.;
        for (int i = 1; i < n; i++)
        {
            System.arraycopy(T1[0], 0, T0[i], 0, n);
            System.arraycopy(T1[0], 0, T1[i], 0, n);
            T0[i][i] = 1.;
        }
        if (!AXB.solve(A_re, A_im, T0, T1, t0, t1, t2, t3))
            return false;
        for (int i = 0; i < n; i++)
        {
            System.arraycopy(T0[i], 0, A_re[i], 0, n);
            System.arraycopy(T1[i], 0, A_im[i], 0, n);
        }
        return true;
    }

    /**
     * Inverts a complex valued square matrix.
     * <p>
     * This method temporary allocates:
     * <ul>
     * <li>4 double arrays</li>
     * <li>2 dual indexed double arrays</li>
     * </ul>
     * <p>
     * @param   A_re    The real part of the square matrix.
     * @param   A_im    The imaginary part of the square matrix.
     * @param   I_re    A matrix taking the real part of the inverse.
     * @param   I_im    A matrix taking the imaginary part of the inverse.
     * @return  True if matrix was regular and false if it was singular.
     */
    public static boolean compute(double[][] A_re, double[][] A_im, double[][] I_re, double[][] I_im)
    {
        int n = A_re[0].length;
        for (int i = 0; i < n; i++)
        {
            System.arraycopy(A_re[i], 0, I_re[i], 0, n);
            System.arraycopy(A_im[i], 0, I_im[i], 0, n);
        }
        return compute(I_re, I_im, null, null, null, null, null, null);
    }
}
