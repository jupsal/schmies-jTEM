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
 * Seperates merged matrices.
 * <p />
 * <em>
 * This class does not check any arguments for validity!<br />
 * For identical arguments results are undefined!
 * </em>
 *
 * @author  Samy Khadem-Al-Charieh
 * @see LR
 * @see PLR
 * @see Householder
 */
public class Unmerge
{
    /**
     * Separates the left and right triangular part out of a merged real matrix.
     * The diagonal goes to the right triangular matrix and the diagonal
     * elements of the left triangular matrix are set to one.
     * Additional rows are ignored.
     * <p />
     * @param   L   The matrix to unmerge and afterwards the left triangular
     *              matrix.
     * @param   R   A matrix to take the right triangular matrix.
     */
    public static void triangular(double[][] L, double[][] R)
    {
        int c = R[0].length;
        for (int i = 0; i < c; i++)
        {
            for (int j = 0; j < i; j++)
                R[i][j] = 0.;
            R[i][i] = L[i][i];
            L[i][i] = 1.;
            for (int j = i + 1; j < c; j++)
            {
                R[i][j] = L[i][j];
                L[i][j] = 0.;
            }
        }
        for (int i = c; i < R.length; i++)
            for (int j = 0; j < c; j++)
                R[i][j] = 0.;
    }

    /**
     * Separates the left and right triangular part out of a merged complex
     * matrix.
     * The diagonal goes to the right triangular matrix and the diagonal
     * elements of the left triangular matrix are set to one.
     * Additional rows are ignored.
     * <p />
     * @param   L_re    The real part of the matrix to unmerge and afterwards
     *                  the left triangular matrix.
     * @param   L_im    The imaginary part of the matrix to unmerge and
     *                  afterwards the left triangular matrix.
     * @param   R_re    A matrix to take the real part of the right triangular
     *                  matrix.
     * @param   R_im    A matrix to take the imaginary part of the right
     *                  triangular matrix.
     */
    public static void triangular(double[][] L_re, double[][] L_im, double[][] R_re, double[][] R_im)
    {
        int c = R_re[0].length;
        for (int i = 0; i < c; i++)
        {
            for (int j = 0; j < i; j++)
            {
                R_re[i][j] = 0.;
                R_im[i][j] = 0.;
            }
            R_re[i][i] = L_re[i][i];
            R_im[i][i] = L_im[i][i];
            L_re[i][i] = 1.;
            L_im[i][i] = 0.;
            for (int j = i + 1; j < c; j++)
            {
                R_re[i][j] = L_re[i][j];
                R_im[i][j] = L_im[i][j];
                L_re[i][j] = 0.;
                L_im[i][j] = 0.;
            }
        }
        for (int i = c; i < R_re.length; i++)
            for (int j = 0; j < c; j++)
            {
                R_re[i][j] = 0.;
                R_im[i][j] = 0.;
            }
    }
}
