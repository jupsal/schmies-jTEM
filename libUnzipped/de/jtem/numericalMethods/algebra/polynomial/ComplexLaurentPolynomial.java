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

package de.jtem.numericalMethods.algebra.polynomial;

  
/** This is provides routines for complex laurent polynomials.
    
    @autor Markus Schmies
*/

public final class ComplexLaurentPolynomial {
  
    private static final double EPS = 1.e-14;

    /** evaluates the polynomial <code>a</code> with the coefficients 
	<code>a_k = aRe[k] + i aIm[k]</code> 
	and degree <code>degOfA</code> at
	<code>z = x + i y</code> and returns the result in <code>aAtZ</code>.
	The polynomial is evaluated using Horners method.	
    */
    public static void eval( final double [] aRe, 
			     final double [] aIm, 
			     final int lowerDegOfA, final int degOfA,
			     final double x, final double y, final double [] aAtZ ) {


	ComplexPolynomial.eval( aRe, aIm, -lowerDegOfA, degOfA, x, y, aAtZ );

	if( lowerDegOfA >= 0 )
	    return;

	final double  xInvOfZ =   x / ( x*x + y*y );
	final double  yInvOfZ = - y / ( x*x + y*y );
    
	double resultRe = aRe[ 0 ];
	double resultIm = aIm[ 0 ];
		
	int i = lowerDegOfA + 1;

	if( lowerDegOfA < -1 ) {

	    for( int index=1; i < 0 && i <= degOfA; i++, index++ ) {

		final double bRe = resultRe; 
		final double bIm = resultIm;
	
		resultRe = bRe * xInvOfZ - bIm * yInvOfZ + aRe[ index ];
		resultIm = bRe * yInvOfZ + bIm * xInvOfZ + aIm[ index ];
	    }
	}

	for( ; i <=0 ; i++ ) {

	    final double bRe = resultRe; 
	    final double bIm = resultIm;

	    resultRe = bRe * xInvOfZ - bIm * yInvOfZ;
	    resultIm = bRe * yInvOfZ + bIm * xInvOfZ;
	}
	
	aAtZ[0] += resultRe;
	aAtZ[1] += resultIm;
    }

    /** evaluates the <code>n</code>-the derivative of the polynomial 
	<code>a</code> with the coefficients 
	<code>a_k = aRe[k + offsetOfA ] + i aIm[k + offsetOfA]</code> 
	and degree <code>degOfA</code> at
	<code>z = x + i y</code> and returns the result in <code>dAAtZ</code>.
	The polynomial is evaluated using Horner�s method. */
    public static void evalDerivative( final double [] aRe, 
				       final double [] aIm, 
				       final int lowerDegOfA, final int degOfA, final int n,
				       final double x, final double y, final double [] dAAtZ ) {

	ComplexPolynomial.evalDerivative
	    ( aRe, aIm, -lowerDegOfA, degOfA, n, x, y, dAAtZ );

	if( lowerDegOfA >= 0 )
	    return;

	final double  xInvOfZ =   x / ( x*x + y*y );
	final double  yInvOfZ = - y / ( x*x + y*y );

	double factor = 1;
	for( int i = 0, j = lowerDegOfA; i<n; i++, j-- )
	    factor *= j;
    
	double resultRe = factor * aRe[ 0 ];
	double resultIm = factor * aIm[ 0 ];
		
	int i = lowerDegOfA + 1;

	if( lowerDegOfA < -1 ) {

	    for( int index=1; i < 0 && i <= degOfA; i++, index++ ) {
		
		factor /= i-n;		
		factor *= i;

		final double bRe = resultRe; 
		final double bIm = resultIm;
	
		resultRe = bRe * xInvOfZ - bIm * yInvOfZ + factor * aRe[ index ];
		resultIm = bRe * yInvOfZ + bIm * xInvOfZ + factor * aIm[ index ];
	    }
	}

	for( ; i <=n ; i++ ) {

	    final double bRe = resultRe; 
	    final double bIm = resultIm;

	    resultRe = bRe * xInvOfZ - bIm * yInvOfZ;
	    resultIm = bRe * yInvOfZ + bIm * xInvOfZ;
	}
	
	dAAtZ[0] += resultRe;
	dAAtZ[1] += resultIm;
    }



    /** computes the nth derivative of given polynomial.
     ** <code>aRe(aIm)</code> and <code>aIm(bIm</code>) may coinside. 
     ** @param aRe real part of coefficients of polynomial that is to differentiate
     ** @param aIm imaginary part of coefficients of polynomial that is to differentiate
     ** @param degOfA degree of polynomial which is to differentiate
     ** @param n number of differentiations
     ** @param bRe real part of the derivative
     ** @param bIm imaginary part of the derivative
     */
    public static void derivative( final double [] aRe, final double [] aIm, 
				   final int lowerDegOfA, final int degOfA,
				   final int n,
				   final double [] bRe, final double [] bIm ) {

	if( lowerDegOfA < 0 ) {

	    double factor = 1;
	    for( int i = 0, j = lowerDegOfA; i<n; i++, j-- )
		factor *= j;
    
	    bRe[ 0 ] = aRe[ 0 ] * factor;
	    bIm[ 0 ] = aIm[ 0 ] * factor;

	    int index = 1;
		
	    for( int i = lowerDegOfA + 1; i < 0 && i <= degOfA; i++, index++ ) {
		
		factor *= i;		
		factor /= i-n;

		 bRe[ index ] = aRe[ index ] * factor;
		 bIm[ index ] = aIm[ index ] * factor;
	    }

	    for( ; index < n - lowerDegOfA && index < degOfA - lowerDegOfA ; index++ ) {
		bRe[ index ] = 0;
		bIm[ index ] = 0;
	    }
	}

	if( degOfA >= n ) {

	    double factor = 1;
	    for( int i = degOfA; i > degOfA-n; i-- )
		factor *= i;

	    int index = degOfA - lowerDegOfA;

	    bRe[ index ] = aRe[ index ] * factor;
	    bIm[ index ] = aIm[ index ] * factor;

	    index--;

	    for( int i=degOfA - 1; i >= n && index >= 0; i--, index-- ) {

		factor *= i+1-n;		
		factor /= i+1;

		bRe[index] = aRe[index] * factor;
		bIm[index] = aIm[index] * factor;
	    }
	}
    }
    
    public static String toString( final double [] re, 
				   final double [] im, final int lowerDegree, final int degree ) {
	return toString( re, im, 0, degree, EPS );
    }
 
    public static String toString( final double [] re, 
				   final double [] im, 
				   final int lowerDegree, final int degree, final double eps ) {

	return ComplexPolynomial.toString( re, im, lowerDegree, degree, eps );
    }
}








