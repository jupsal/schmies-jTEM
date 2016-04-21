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

package de.jtem.numericalMethods.calculus.functionApproximation;

import de.jtem.numericalMethods.function.DoubleParametrized;
import de.jtem.numericalMethods.function.DoubleValued;

/**
 * A class that provides static methods related to the approximation of
 * functions with Chebyshev polynomials.  
 *
 * <p> The Chebyshev polynomial of degree <code>n</code> is defined by </p>
 *
 * <p align=center> <code>T<sub>n</sub>(cos t) = cos(nt)</code>. </p>
 *
 * <p> 
 *   The Chebyshev polynomials satisfy the recursion formulas
 * </p>
 *
 * <a name=recursion>
 * <p>
 *   <table align=center border="0">
 *     <tr>
 *       <td>
 *         <code>T<sub>0</sub>(x)</code>
 *       </td>
 *       <td>
 *         <code>=&ensp;</code>
 *       </td>
 *       <td>
 *         <code>1</code>
 *       </td>
 *     </tr>
 *     <tr>
 *       <td>
 *         <code>T<sub>1</sub>(x)</code>
 *       </td>
 *       <td>
 *         <code>=&ensp;</code>
 *       </td>
 *       <td>
 *         <code>x</code>
 *       </td>
 *     </tr>
 *     <tr>
 *       <td>
 *         <code>T<sub>n+2</sub>(x)</code>
 *       </td>
 *       <td>
 *         <code>=&ensp;</code>
 *       </td>
 *       <td>
 *         <code>2x T<sub>n+1</sub>(x) - T<sub>n</sub>(x)</code>
 *       </td>
 *     </tr>
 *   </table>
 * </p>
 * </a>
 *
 * <p>
 *   They are orthogonal over the interval <code>[-1,1]</code> with weight 
 *   function <code>(1 - x<sup>2</sup>)<sup>-1/2</sup></code>. More precisely,
 * </p>
 * <p>
 *   <table border="0" align=center>
 *     <tr>
 *       <td>
 *         &nbsp;
 *       </td>
 *       <td>
 *         <code>0</code>,
 *       </td>
 *       <td>
 *         if <code>m&ne;n</code>
 *       </td>
 *     </tr>
 *     <tr>
 *       <td>
 *         <code>
 *           &int;<sub>-1</sub><sup>1</sup> 
 *           (1 - x<sup>2</sup>)<sup>-1/2</sup> 
 *           T<sub>n</sub>(x) T<sub>m</sub>(x) dx =
 *         </code>&emsp;
 *       </td>
 *       <td>
 *         <code>&pi;/2</code>,&emsp; 
 *       </td>
 *       <td>
 *         if <code>m=n&ne;0</code>&ensp;.
 *       </td>
 *     </tr>
 *     <tr>
 *       <td>
 *         &nbsp;
 *       </td>
 *       <td>
 *         <code>&pi</code>,
 *       </td>
 *       <td>
 *         if <code>m=n=0</code>
 *       </td>
 *     </tr>
 *   </table>
 * </p>
 *
 * <p> The Chebyshev polynomials also satisfy a discrete orthogonality
 * relation. If <code>x<sub>j</sub> = cos(&pi;(2j-1)/2n)</code>,
 * <code>j=1,&hellip; n</code>, are the <code>n</code> zeroes of
 * <code>T<sub>n</sub>(x)</code>, then </p>
 *
 * <p>
 *   <table align=center>
 *     <tr>
 *       <td>
 *         &nbsp;
 *       </td>
 *       <td>
 *         0,
 *       </td>
 *       <td>
 *         if <code>m&ne;k</code> 
 *     </tr>
 *     <tr>
 *       <td>
 *         <code>
 *           &sum<sub>j=1</sub><sup>n</sup> 
 *           T<sub>m</sub>(x<sub>j</sub>)
 *           T<sub>k</sub>(x<sub>j</sub>) = &emsp;
 *         </code>
 *       </td>
 *       <td>
 *         <code>&pi;/2</code>,&emsp;
 *       </td>
 *       <td>
 *         if <code>m=k&ne;0</code>&ensp;.
 *       </td>
 *     </tr>
 *     <tr>
 *       <td>
 *         &nbsp;
 *       </td>
 *       <td>
 *         <code>&pi;</code>,
 *       </td>
 *       <td>
 *         if <code>m=k=0</code> 
 *     </tr>
 *   </table>
 * </p>
 *
 * <p> Suppose <code>f(x)</code> is a continuous function on the interval
 * <code>[0,1]</code>. Let the <code><a name=n>n</a></code> coefficients <a
 * name=cj><code>c<sub>j</sub></code></a>, <code>j=0,&hellip; n-1</code> be
 * defined by </p>
 *
 * <p align=center>
 *   <code>
 *     c<sub>j</sub> = (2/n) &sum;<sub>k=1</sub><sup>n</sup>
 *     f(x<sub>k</sub>) T<sub>j</sub>(x<sub>k</sub>)
 *   </code>.
 * </p>
 * 
 * <p>Then</p>
 *
 * <a name=series>
 *   <p align=center>
 *     <code>
 *       f(x) &asymp; c<sub>0</sub>/2 
 *       + &sum;<sub>j=1</sub><sup>n-1</sup> c<sub>j</sub> T<sub>j</sub>(x).
 *     </code>
 *   </p>
 * </a>
 *
 * <p> This approximation is exact for the zeros <code>x<sub>j</sub></code>
 * of <code>T<sub>n</sub>.</p>
 * 
 * @see "Numerical Recipes in C, chapters 5.8&ndash;5.9"
 *
 * @author Boris Springborn */
public class ChebyshevApproximation {

    private ChebyshevApproximation() {}

    /**
     * Calculates the coefficients <a
     * href=#cj><code>c<sub>j</sub></code></a>, with <code><a href=#n>n</a> =
     * c.length</code>, for the function implemented by <code>p</code> and
     * <code>v</code>. The arguments <code>p</code> and <code>v</code> may
     * (and usually will) be the same object.
     *
     * @param c a <code>double</code> array.  When the method returns, it
     *          holds the coefficients <a
     *          href=#cj><code>c<sub>j</sub></code></a>, where <code><a
     *          href=#n>n</a> = c.length</code>.
     * @param p an object implementing the interface 
     *                 {@link DoubleParametrized}. 
     *          Only the method 
     *          {@link DoubleParametrized#setDoubleParameter(double) setDoubleParameter(double x)} is called with <code>-1 &le; x &le; 1</code>.
     * @param v an object implementing the interface {@link DoubleValued}. 
     *          Only the method {@link DoubleValued#getDoubleValue() getDoubleValue()} 
     *          is called.  */
    public static void fit(double[] c, RealFunction f ) {

        final int n = c.length;
        final double twoOverN = 2.0 / n;
        final double piOverN = Math.PI / n;

        double[] valueTable = new double[n];

        for (int i = 0; i < n; i++) {
            valueTable[i] = f.valueAt(Math.cos(piOverN * (i + 0.5)) );
        } 

        for (int j = 0; j < n; j++) {
            double sum = 0.0;
            for (int i = 0; i < n; i++) {
                sum += valueTable[i] * Math.cos(piOverN * j * (i + .5));
            } 
            c[j] = twoOverN * sum;
        } 

    } 


    /**
     * Evaluates the <a href=#series>Chebyshev series</a> with coefficients 
     * <code>c<sub>j</sub></code>. The <code>x</code>-value at which
     * the series is to be evaluated has to be in the range <nobr><code>-1
     * &le; x &le; 1</code></nobr>. Otherwise, an 
     * {@link java.lang.IllegalArgumentException IllegalArgumentException}
     * is thrown. Note that the coefficient of <code>T<sub>0</sub></code> 
     * is <code>c<sub>0</sub>/2</code>. This makes this method compatible 
     * with the method {@link #fit fit}.
     * 
     * <p>The <a href=#recursion>recursion formula</a> for the Chebyshev
     * polynomials is utilized to calculate the series without having to
     * calculate each Chebyshev polynomial individually. </p>
     *
     * @param c a <code>double</code> array holding the coefficents 
     *          <code>c<sub>j</sub></code>.
     * @param x a <code>double</code>. The <code>x</code>-value at which the
     *          series is to be evaluated. Must be in the range <code>-1 &le;
     *          x &le; 1</code>.
     * @return  the sum of the 
     *          <a href=#series>Chebyshev series</a>. 
     * @throws  IllegalArgumentException 
     *            if <code>x</code> is not in the 
     *            range <code>-1 &le; x &le; 1</code>.  */
    public static double evaluate(double[] c, double x) 
        throws java.lang.IllegalArgumentException {

        if ((x < -1.0) || (x > 1.0)) {
            throw new java.lang.IllegalArgumentException 
                ("Argument x = " + x + " not in range: " +
                 "-1.0 <= x <= 1.0 expected.");
        } 

        double a = 0.0;
        double b = 0.0;

        for (int i = c.length - 1; i >= 0; i--) {
            final double bNew = c[i] - a;
            a = 2.0 * x * a + b;
            b = bNew;
        }

        return a * x + b - c[0] / 2.0;

    } 

    /**
     * Simply calls {@link #integrate(double[], double[], double) integate(c, cInt, factor)} with <code>double factor = 1.0</code>. 
     */
    public static void integrate(double [] c, double [] cInt) {
        integrate(c, cInt, 1.0);
    }

    /**
     * Given the coefficients <code>c[j]</code> of a Chebyshev series
     * approximating some function <code>f(x)</code>, calculates the
     * coefficients <code>cInt[j]</code> of a Chebyshev series approximating
     * <nobr><code>factor &times; &int;&thinsp;f(x)dx</code></nobr>. The
     * constant of integration is chosen so that the value at <code>x =
     * 0</code> is <code>0</code>.
     *
     * <p> The arrays <code>c</code> and <code>cInt</code> must have the same
     * length, which must be at least <code>2</code>. </p> 
     *
     * @param c a <code>double</code> array of length at least 2. Holds the
     *        coefficients of a Chebyshev series approximating some function
     *        <code>f(x)</code>.
     * @param cInt a <code>double</code> array of the same length as
     *        <code>c</code>. After the method returns, it holds the
     *        coeffiecents of a Chebyshev series approximating
     *        <nobr><code>factor &times; * &int;&thinsp;f(x)dx</code></nobr>.
     * @param factor a <code>double</code> with which the integral is
     *        multiplied. Useful in the case of a change of variable.  
     * @throws IllegalArgumentException if 
     *        <code>c.length &ne; cInt.length</code>.  
     */
    public static void integrate(double [] c, double [] cInt, double factor)
        throws java.lang.IllegalArgumentException {

        final int n = c.length;

        if (cInt.length != n) {
            throw new IllegalArgumentException
                ("Incompatible Arguments: c.length == cInt.length expected.");
        }
        
        final double factorHalf = 0.5 * factor;

        for (int i = 1; i < n - 1; i++) {
            cInt [i] = factorHalf * (c [i - 1] - c [i + 1]) / i;
        } 

        cInt[n - 1] = factorHalf * c[n - 2] / (n - 1);

        cInt [0] -= 2 * evaluate(cInt, 0.0);
        
    }


    /**
     * Called with coefficients <code>c[j]</code> of a Chebyshev series which
     * sums to a function <code>f(x)</code> with <code>f(0) = 0</code>,
     * returns with <code>c</code> holding the coefficients of the Chebyshev
     * expansion of <code>f(x)/x</code>. Hence, when it returns, <code>c[n-1]
     * = 0.0</code>, where <code>n = c.length</code>.
     *
     * @param c a <code>double</code> array of length at least 1, holding the
     *          coefficients of a Chebyshev series which sums to a function
     *          <code>f(x)</code> with <code>f(0) = 0</code>. When the method
     *          returns, <code>c</code> holds the coefficients of the
     *          Chebyshev expansion of <code>f(x)/x.</code> */
    public static void divideByX(double [] c) {

        double cIPlusTwo = 0.0;
        double cIPlusOne = 0.0;

        for (int i = c.length - 1; i > 0; i--) {
            c [i] = 2 * c [i] - cIPlusTwo;
            cIPlusTwo = cIPlusOne;
            cIPlusOne = c [i];
        }

        for (int i = 0; i < c.length - 1; i++) {
            c [i] = c [i + 1];
        }

        c [c.length - 1] = 0;

    }

}
