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

/** 
 * Implementation of the real (incomplete) Gamma function.
 * <p>
 * This class provides an implementaion of the real Gamma function
 * <p align=center>
 * <code>&Gamma;( x ) = &int;<sub>0</sub><sup>&infin;</sup> t<sup>x-1</sup>e<sup>-t</sup>dt</code>
 * </p>
 * and its sibling, the <i>"upper" incomplete</i> Gamma function:
 * <p align=center>
 * <code>&Gamma;( a, x ) = &int;<sub>x</sub><sup>&infin;</sup>t<sup>a-1</sup>e<sup>-t</sup>dt</code>
 * </p><p>
 * The Implementation follows the ideas of C.Lanczos and their realizations
 * presented in <em>Numerical Recipes in C</em>.
 * </p><p> 
 * Further reading:
 * <table>
 * <tr><td> - </td><td>C.Lanczos, 1964, SIAM Journal of Numerical Analysis, ser. B, vol. 1,
 * pp. 86-96.</td></tr>
 * <tr><td> - </td><td>W.H.Press, S.A.Teukolsky, W.T.Vetterling, B.P.Flannery.
 * <em>Numerical Recipes in C</em>. Cambridge University Press, 1992. </td></tr>
 * </table>
 *</p>
 */
public class Gamma {

    private static final double a = Math.sqrt( 2*Math.PI ) *    1.000000000190015;
    private static final double b = Math.sqrt( 2*Math.PI ) *   76.18009172947146;
    private static final double c = Math.sqrt( 2*Math.PI ) *  -86.50532032941677;
    private static final double d = Math.sqrt( 2*Math.PI ) *   24.01409824083091;
    private static final double e = Math.sqrt( 2*Math.PI ) *   -1.231739572450155;
    private static final double f = Math.sqrt( 2*Math.PI ) *    0.1208650973866179e-2;
    private static final double g = Math.sqrt( 2*Math.PI ) *   -0.5395239384953e-5;

    private Gamma() {}

    static double gamma_( double z ) {

	if( z <= 0 )
	    throw new IllegalArgumentException("argument must be positive");

	double s = Math.exp( Math.log( z + 5.5 ) * ( z + 0.5 ) - ( z + 5.5 ) );
	
	return s * ( a + 
		     b/(z+1) + c/(z+2) + d/(z+3) + 
		     e/(z+4) + f/(z+5) + g/(z+6)  ) / z;
    }

    /**
     * Returns the gamma function at <code>x</code>.
     * @param x real number
     * @return <code> &Gamma;( x ) = &int;<sub>0</sub><sup>&infin;</sup> 
     * t<sup>x-1</sup>e<sup>-t</sup>dt</code>
     */
    public static double gamma( double x ) {

	return x <= 0.5 ? Math.PI / gamma_ ( 1 - x ) / Math.sin( Math.PI * x ) : gamma_ ( x );
    }

    static double logOfGamma_( double z ) {

	if( z <= 0 )
	    throw new IllegalArgumentException("argument must be positive");

	double s = Math.log( z + 5.5 ) * ( z + 0.5 ) - ( z + 5.5 );
	
	return s + Math.log(  ( a + 
				b/(z+1) + c/(z+2) + d/(z+3) + 
				e/(z+4) + f/(z+5) + g/(z+6)  ) / z );
    }


    /**
     * Returns the logarithm of the gamma function at <code>x</code>.
     * @param x real number
     * @return <code>log( &Gamma;( x ) ) = log ( &int;<sub>0</sub><sup>&infin;</sup> 
     * t<sup>x-1</sup>e<sup>-t</sup>dt )</code>
     */
    public static double logOfGamma( double x ) {

	return x <= 0.5 
	    ? Math.log( Math.PI / Math.sin( Math.PI * x ) ) - logOfGamma_ ( 1 - x )  
	    : logOfGamma_ ( x );
    
    }


    private static double gamma( double a, double x, double gammaOfA, boolean computeGammaOfA ) {

	if( a<=0 )
	    throw new IllegalArgumentException( "a is not positive" );

	if( x<0 )
	    throw new IllegalArgumentException( "x is negative" );

	if( x < a + 1 ) {
	    if( computeGammaOfA )
		gammaOfA = gamma ( a );

	    return computeViaSeries( a, x, gammaOfA );
	
	} else
	    return computeViaContinuedFraction( a, x );
    }

    /**
     * Returns the <i>"upper" incomplete</i> gamma function at <code>a, x</code>.
     * The algorithms for the incomplete gamma function needs
     * to compute <code>&Gamma;(a)</code>; providing this value
     * will therefor lead to an optimization if you need
     * this value also for a diffrent porpuse.
     * @param a positive real number
     * @param x non-negative real number
     * @return <code>&Gamma;( a, x ) = &int;<sub>x</sub><sup>&infin;
     * </sup>t<sup>a-1</sup>e<sup>-t</sup>dt</code>
     */
    public static double gamma( double a, double x, double gammaOfA ) {
    	return gamma ( a, x, gammaOfA, false );
    }


    /**
     * Returns the <i>"upper" incomplete</i> gamma function at <code>a, x</code>.
     * @param a positive real number
     * @param x non-negative real number
     * @return <code>&Gamma;( a, x ) = &int;<sub>x</sub><sup>&infin;
     * </sup>t<sup>a-1</sup>e<sup>-t</sup>dt</code>
     */
    public static double gamma( double a, double x ) {
	return gamma ( a, x, 123, true );
    }
    
    static int MAX_NUMBER_OF_ITERATIONS = 200;

    static double EPS = 1e-16;

    static double computeViaSeries( double a, double x, double gammaOfA ) {

	double factor = a;
	double delta  = 1 / a;
	double sum    = delta;
	
	int counter   = 0;

	while( Math.abs(delta) > Math.abs( sum ) * EPS ) {
	    
	    delta *= x/++factor;

	    sum   += delta;

	    if( counter > MAX_NUMBER_OF_ITERATIONS )
		throw new RuntimeException( "exeeded max number of iterations="+
					    MAX_NUMBER_OF_ITERATIONS );
	}

	return gammaOfA  - sum * Math.exp( a*Math.log( x ) - x );
    }

    static double MIN_DP = 1e-30;

    static double computeViaContinuedFraction( double a, double x ) {

	    double b = x + 1 - a;
	    double c = 1/MIN_DP;
	    double d = 1/b;
	    double h = d;

	    for( int i=1; i<MAX_NUMBER_OF_ITERATIONS; i++ ) {

		double an = -i * (i-a);

		b += 2;

		d  = an*d + b;
		c  = an/c + b;

		if( Math.abs( c ) < MIN_DP )
		    c = MIN_DP;
		if( Math.abs( d ) < MIN_DP )
		    d = MIN_DP;

		d = 1 / d;

		double delta = c * d;

		h *= delta;

		if( Math.abs( delta - 1 ) < EPS )
		    return h * Math.exp( a * Math.log( x ) - x );
	    }

	    throw new RuntimeException( "exeeded max number of iterations="+
					MAX_NUMBER_OF_ITERATIONS  );
    }
    
}

