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

import de.jtem.numericalMethods.algebra.linear.decompose.PLR;

/**
 * Computes the determinant of a square matrix.
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
 * Method versions for Determinant:
 * <blockquote>
 * <table border="2" cellpadding="6">
 * <tr>
 * <th></th><th>optimized</th><th>easy to use</th>
 * </tr>
 * <tr>
 * <th>real valued</th>
 * <td>{@link #compute(double[][],int[]) compute(A,...)}</td>
 * <td>{@link #compute(double[][]) compute(A)}</td>
 * </tr>
 * <tr>
 * <th>complex valued</th>
 * <td>{@link #compute(double[][],double[][],int[]) compute(A,...)}</td>
 * <td>{@link #compute(double[][],double[][]) compute(A)}</td>
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
public class Determinant
{
    /**
     * Computes the determinant of a real valued square matrix.
     * <p>
     * @param   A   The square matrix and destroyed afterwards.
     * @param   t   An int array of same size as matrix A has rows for storing
     *              temporary results during the computation.
     * @return  The determinant.
     */
    public static double compute(double[][] A, int[] t)
    {
        if (t == null)
            t = new int[A.length];
        return PLR.decompose(t, A);
    }

    /**
     * Computes the determinant of a real valued square matrix.
     * <p>
     * This method temporary allocates:
     * <ul>
     * <li>1 int array</li>
     * <li>1 dual indexed double array</li>
     * </ul>
     * <p>
     * @param   A   The square matrix.
     * @return  The determinant.
     */
    public static double compute(double[][] A)
    {
        int n = A.length;
        double[][] T = new double[n][n];
        for (int i = 0; i < n; i++)
            System.arraycopy(A[i], 0, T[i], 0, n);
        return compute(T, (int[]) null);
    }

    /**
     * Computes the determinant of a complex valued square matrix.
     * <p>
     * @param   A_re    The real part of the square matrix and destroyed
     *                  afterwards.
     * @param   A_im    The imaginary part of the square matrix and destroyed
     *                  afterwards.
     * @param   t       An int array of same size as matrix A has rows for
     *                  storing temporary results during the computation.
     * @return  The determinant.
     */
    public static double[] compute(double[][] A_re, double[][] A_im, int[] t)
    {
        if (t == null)
            t = new int[A_re.length];
        return PLR.decompose(t, A_re, A_im);
    }

    /**
     * Computes the determinant of a complex valued square matrix.
     * <p>
     * This method temporary allocates:
     * <ul>
     * <li>1 int array</li>
     * <li>2 dual indexed double arrays</li>
     * </ul>
     * <p>
     * @param   A_re    The real part of the square matrix.
     * @param   A_im    The imaginary part of the square matrix.
     * @return  The determinant.
     */
    public static double[] compute(double[][] A_re, double[][] A_im)
    {
        int n = A_re.length;
        double[][] T_re = new double[n][n];
        double[][] T_im = new double[n][n];
        for (int i = 0; i < n; i++)
        {
            System.arraycopy(A_re[i], 0, T_re[i], 0, n);
            System.arraycopy(A_im[i], 0, T_im[i], 0, n);
        }
        return compute(T_re, T_im, null);
    }
}
