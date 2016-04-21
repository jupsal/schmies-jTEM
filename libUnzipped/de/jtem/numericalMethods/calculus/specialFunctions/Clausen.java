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

package de.jtem.numericalMethods.calculus.specialFunctions;

import de.jtem.numericalMethods.calculus.functionApproximation.ChebyshevApproximation;
import de.jtem.numericalMethods.calculus.functionApproximation.RealFunction;

/**
 * <p>A class providing Clausen's integral and Catalans constant. <a
 * name=clausen>Clausen's integral</a> is the <code>2&pi;</code>-periodic
 * function</p>
 *
 * <p align=center>
 *   <code>
 *     Cl<sub>2</sub>(x) = &sum;<sub>1</sub><sup>&infin;</sup>
 *     sin(nx)/ n<sup>2</sup> 
 *   </code>.
 * </p>
 *
 * <p><a name=catalan>Catalan's constant</a> is</p>
 *
 * <p align=center>
 *   <code>
 *     <table>
 *       <tr>
 *         <td> G </td>
 *         <td>= Cl<sub>2</sub>(&pi;/2)</td>
 *       </tr>
 *       <tr>
 *         <td>&nbsp;</td>
 *         <td>= 1/1<sup>2</sup> - 1/3<sup>2</sup> 
 *               + 1/5<sup>2</sup> - 1/7<sup>2</sup> +- &hellip;</td>
 *       </tr>
 *       <tr>
 *         <td>&nbsp;</td>
 *         <td>= 0.915965594177&hellip;</td>
 *       </tr>
 *     </table>
 *   </code>
 * </p>
 *
 * <p>For <code>0 &le; x &le; 2&pi</code>, the following formula (which
 * explains the name "Clausen's <em>integral</em>") holds.</p>
 *
 * <p align=center>
 *   <code>
 *     Cl<sub>2</sub>(x) = - &int;<sub>0</sub><sup>x</sup>
 *     log 2sin (&xi;/2) d&xi;
 *   </code>
 * </p>
 *
 * <p>Clausen's integral is almost the same as Milnor's Lobachevski
 * function</p>
 *
 * <p align=center>
 *   <code>
 *     &#x041B;(x) = - &int;<sub>0</sub><sup>x</sup>
 *     log |2 sin (&xi;)| d&xi; = Cl<sub>2</sub>(2x)/ 2 
 *   </code>.
 * </p>
 * 
 *
 * @see "Leonard Lewin. <em>Polylogarithms and Associated Functions.</em>
 * North Holland, New York, 1981."
 *
 * @author Boris Springborn */
public class Clausen {

    private Clausen() {}

    /* 
     * These Chebyshev coefficients were calculated by main with n=100.
     */
    private static final double [] coeffs ={
        0.4482875649180977,
        0.22723902480478247,
        0.003193192528191099,
        1.0195254591199122E-4,
        4.188220138620993E-6,
        1.951947102494479E-7,
        9.833314612744994E-9,
        5.226269539038555E-10,
        2.889166021995073E-11,
        1.6465290606931503E-12,
        9.57374364530389E-14
    };

    /**
     * <a href=#catalan>Catalan's constant</a> <code>G</code>.
     */
    public static final double CATALAN = Clausen.cl2 (0.5 * Math.PI);

    /**
     * <a href=#clausen>Clausens integral</a> <code>Cl<sub>2</sub>.
     *
     * @param x a double
     * @return <code> Cl<sub>2</sub>(x)</code>
     */
    public static double cl2 (double x) {
        /* Since Clausen's integral is 2pi-periodic, we first shift x into
         * the range -pi<x<=pi, using Math.IEEEremainder. Around x=0,
         * Clausens integral is asymptotically x(log|x|-1). We subtract this
         * term and approximate the rest by a Chebyshev series. */
        x = Math.IEEEremainder (x, 2.0 * Math.PI);

        if (x == 0.0) {
             return 0.0;
        }
        else {
            final double xOverPi = x / Math.PI;
            final double y = xOverPi * ChebyshevApproximation.evaluate 
                (coeffs, 2*xOverPi*xOverPi - 1);
            return y - x * (Math.log (Math.abs (x)) - 1);
        }
    }


    /* Called by main */
    private static void printChebychevCoeffificients (int n) {

        class MinusLogSinOverX 
            implements RealFunction {
            	
            public double valueAt(double x){
                if (x == 0.0) {
                    return 0;
                }
                else {
                    final double xx = 0.5 * Math.PI * x;
                    return - Math.log (Math.sin (xx) / xx);
                }
            }
        }

        MinusLogSinOverX function = new MinusLogSinOverX ();

        double [] c1 = new double [2 * n];

        // Calculate Chebyshev coefficients for 'function'.
        ChebyshevApproximation.fit (c1, function);

       // Set odd coefficients (which should be zero) to exactly zero.
        for (int i = 1; i < c1.length; i+=2) {
            c1 [i] = 0.0;
        }

        // Integrate Chebyshev series.
        double [] c2 = new double [2 * n];;
        ChebyshevApproximation.integrate (c1, c2, Math.PI);

        // Calculate coefficients for 'integral/x' to get even series. 
        ChebyshevApproximation.divideByX (c2);

        // Output
        System.out.println ("Paste these coefficients into the code.");
        System.out.println ("{");
        for (int i = 0; i < c2.length - 2; i += 2) {
            System.out.println (c2 [i] + ",");
        }
        System.out.println (c2 [c2.length - 2]);
        System.out.println ("}");
        
    }


    /**
     * Prints Chebyshev coefficients to <code>stdout</code>, which are to be
     * pasted into the code. The number of coefficents which are
     * calculated may be passed as argument. The default is 25.
     * 
     * @param args a <code>String</code> array. If <code>args.length >
     * 0</code>, then <code>args[0]</code> is expected to have an integer
     * value, which determines the number of coefficients calculated. */
    public static void main(String[] args) {
        int n;
        if (args.length == 0) {
            n = 25;
        } else {
            n = (Integer.valueOf(args[0])).intValue();
        }
        printChebychevCoeffificients (n);
    }

}

