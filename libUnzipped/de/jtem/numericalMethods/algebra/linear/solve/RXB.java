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

/**
 * Computes for a system of equations RX=B with given right triangular matrix
 * R and matrix B the solving matrix X.
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
 * Method versions for Rx=b:
 * <blockquote>
 * <table border="2" cellpadding="6">
 * <tr>
 * <th></th><th>optimized</th><th>easy to use</th>
 * </tr>
 * <tr>
 * <th>real valued</th>
 * <td>{@link #solve(double[][],double[]) solve(R,b)}</td>
 * <td>{@link #solve(double[][],double[],double[]) solve(R,x,b)}</td>
 * </tr>
 * <tr>
 * <th>complex valued</th>
 * <td>{@link #solve(double[][],double[][],double[],double[]) solve(R,b)}</td>
 * <td>{@link #solve(double[][],double[][],double[],double[],double[],double[]) solve(R,x,b)}</td>
 * </tr>
 * </table>
 * </blockquote>
 * Method versions for RX=B:
 * <blockquote>
 * <table border="2" cellpadding="6">
 * <tr>
 * <th></th><th>optimized</th><th>easy to use</th>
 * </tr>
 * <tr>
 * <th>real valued</th>
 * <td>{@link #solve(double[][],double[][]) solve(R,B)}</td>
 * <td>{@link #solve(double[][],double[][],double[][]) solve(R,X,B)}</td>
 * </tr>
 * <tr>
 * <th>complex valued</th>
 * <td>{@link #solve(double[][],double[][],double[][],double[][]) solve(R,B)}</td>
 * <td>{@link #solve(double[][],double[][],double[][],double[][],double[][],double[][]) solve(R,X,B)}</td>
 * </tr>
 * </table>
 * </blockquote>
 * <em>
 * This class does not check any arguments for validity!<br />
 * For identical arguments results are undefined!
 * </em>
 *
 * @author  Samy Khadem-Al-Charieh
 */
public class RXB
{
    /**
     * Solves a system of real equations.
     * <p />
     * @param   R   The right triangular matrix.
     * @param   b   The right side and the solution afterwards (for non square
     *              left sides the lower part is ignored).
     * @return  True if system was solvable.
     */
    public static boolean solve(double[][] R, double[] b)
    {
        int c = R[0].length;
        for (int i = c - 1; i >= 0; i--)
        {
            for (int k = i + 1; k < c; k++)
                b[i] -= b[k] * R[i][k];
            double t = R[i][i];
            if (t == 0.)
                return false;
            b[i] /= t;
        }
        return true;
    }

    /**
     * Solves a system of real equations.
     * <p />
     * @param   R   The right triangular matrix.
     * @param   x   A vector taking the solution.
     * @param   b   The right side (for non square left sides the lower part is
     *              ignored).
     * @return  True if system was solvable.
     */
    public static boolean solve(double[][] R, double[] x, double[] b)
    {
        System.arraycopy(b, 0, x, 0, x.length);
        return solve(R, x);
    }

    /**
     * Solves a system of complex equations.
     * <p />
     * @param   R_re    The real part of the right triangular matrix.
     * @param   R_im    The imaginary part of the right triangular matrix.
     * @param   b_re    The real part of the right side and the solution
     *                  afterwards (for non square left sides the lower part is
     *                  ignored).
     * @param   b_im    The imaginary part of the right side and the solution
     *                  afterwards (for non square left sides the lower part is
     *                  ignored).
     * @return  True if system was solvable.
     */
    public static boolean solve(double[][] R_re, double[][] R_im, double[] b_re, double[] b_im)
    {
        int c = R_re[0].length;
        for (int i = c - 1; i >= 0; i--)
        {
            for (int k = i + 1; k < c; k++)
            {
                double x_re = b_re[k];
                double x_im = b_im[k];
                double y_re = R_re[i][k];
                double y_im = R_im[i][k];
                b_re[i] -= x_re * y_re - x_im * y_im;
                b_im[i] -= x_im * y_re + x_re * y_im;
            }
            double y_re = R_re[i][i];
            double y_im = - R_im[i][i];
            double t = y_re * y_re + y_im * y_im;
            if (t == 0.)
                return false;
            double x_re = b_re[i];
            double x_im = b_im[i];
            b_re[i] = (x_re * y_re - x_im * y_im) / t;
            b_im[i] = (x_im * y_re + x_re * y_im) / t;
        }
        return true;
    }

    /**
     * Solves a system of complex equations.
     * <p />
     * @param   R_re    The real part of the right triangular matrix.
     * @param   R_im    The imaginary part of the right triangular matrix.
     * @param   x_re    A vector taking the real part of the solution.
     * @param   x_im    A vector taking the imaginary part of the solution.
     * @param   b_re    The real part of the right side (for non square left
     *                  sides the lower part is ignored).
     * @param   b_im    The imaginary part of the right side (for non square
     *                  left sides the lower part is ignored).
     * @return  True if system was solvable.
     */
    public static boolean solve(double[][] R_re, double[][] R_im, double[] x_re, double[] x_im, double[] b_re, double[] b_im)
    {
        int n = x_re.length;
        System.arraycopy(b_re, 0, x_re, 0, n);
        System.arraycopy(b_im, 0, x_im, 0, n);
        return solve(R_re, R_im, x_re, x_im);
    }

    /**
     * Solves a system of real equations.
     * <p />
     * @param   R   The right triangular matrix.
     * @param   B   The right side and the solution afterwards (for non square
     *              left sides the lower part is ignored).
     * @return  True if system was solvable.
     */
    public static boolean solve(double[][] R, double[][] B)
    {
        int c = R[0].length;
        int n = B[0].length;
        for (int i = c - 1; i >= 0; i--)
        {
            for (int k = i + 1; k < c; k++)
                for (int j = 0; j < n; j++)
                    B[i][j] -= B[k][j] * R[i][k];
            double t = R[i][i];
            if (t == 0.)
                return false;
            for (int j = 0; j < n; j++)
                B[i][j] /= t;
        }
        return true;
    }

    /**
     * Solves a system of real equations.
     * <p />
     * @param   R   The right triangular matrix.
     * @param   X   A matrix taking the solution.
     * @param   B   The right side (for non square left sides the lower part is
     *              ignored).
     * @return  True if system was solvable.
     */
    public static boolean solve(double[][] R, double[][] X, double[][] B)
    {
        int n = X[0].length;
        for (int i = 0; i < X.length; i++)
            System.arraycopy(B[i], 0, X[i], 0, n);
        return solve(R, X);
    }

    /**
     * Solves a system of complex equations.
     * <p />
     * @param   R_re    The real part of the right triangular matrix.
     * @param   R_im    The imaginary part of the right triangular matrix.
     * @param   B_re    The real part of the right side and the solution
     *                  afterwards (for non square left sides the lower part is
     *                  ignored).
     * @param   B_im    The imaginary part of the right side and the solution
     *                  afterwards (for non square left sides the lower part is
     *                  ignored).
     * @return  True if system was solvable.
     */
    public static boolean solve(double[][] R_re, double[][] R_im, double[][] B_re, double[][] B_im)
    {
        int c = R_re[0].length;
        int n = B_re[0].length;
        for (int i = c - 1; i >= 0; i--)
        {
            for (int k = i + 1; k < c; k++)
                for (int j = 0; j < n; j++)
                {
                    double a_re = B_re[k][j];
                    double a_im = B_im[k][j];
                    double b_re = R_re[i][k];
                    double b_im = R_im[i][k];
                    B_re[i][j] -= a_re * b_re - a_im * b_im;
                    B_im[i][j] -= a_im * b_re + a_re * b_im;
                }
            double b_re = R_re[i][i];
            double b_im = - R_im[i][i];
            double t = b_re * b_re + b_im * b_im;
            if (t == 0.)
                return false;
            for (int j = 0; j < n; j++)
            {
                double a_re = B_re[i][j];
                double a_im = B_im[i][j];
                B_re[i][j] = (a_re * b_re - a_im * b_im) / t;
                B_im[i][j] = (a_im * b_re + a_re * b_im) / t;
            }
        }
        return true;
    }

    /**
     * Solves a system of complex equations.
     * <p />
     * @param   R_re    The real part of the right triangular matrix.
     * @param   R_im    The imaginary part of the right triangular matrix.
     * @param   X_re    A matrix taking the real part of the solution.
     * @param   X_im    A matrix taking the imaginary part of the solution.
     * @param   B_re    The real part of the right side (for non square left
     *                  sides the lower part is ignored).
     * @param   B_im    The imaginary part of the right side (for non square
     *                  left sides the lower part is ignored).
     * @return  True if system was solvable.
     */
    public static boolean solve(double[][] R_re, double[][] R_im, double[][] X_re, double[][] X_im, double[][] B_re, double[][] B_im)
    {
        int n = X_re[0].length;
        for (int i = 0; i < X_re.length; i++)
        {
            System.arraycopy(B_re[i], 0, X_re[i], 0, n);
            System.arraycopy(B_im[i], 0, X_im[i], 0, n);
        }
        return solve(R_re, R_im, X_re, X_im);
    }
}
