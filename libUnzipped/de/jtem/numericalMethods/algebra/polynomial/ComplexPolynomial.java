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

/** This is provides routines for complex polynomials.

    @autor Piotr Greniv, Markus Schmies
*/

public final class ComplexPolynomial {

    private static final double EPS = 1.e-14;

    /** evaluates the polynomial <code>a</code> with the coefficients
	<code>a_k = aRe[k] + i aIm[k]</code>
	and degree <code>degOfA</code> at
	<code>z = x + i y</code> and returns the result in <code>aAtZ</code>.
	The polynomial is evaluated using Horner\uFFFDs method. */
    public static void eval( final double [] aRe,
			     final double [] aIm, final int degOfA,
			     final double x, final double y, final double [] aAtZ ) {
	eval( aRe, aIm, 0, degOfA, x, y, aAtZ );
    }

    /** evaluates the polynomial <code>a</code> with the coefficients
	<code>a_k = aRe[k + offsetOfA ] + i aIm[k + offsetOfA]</code>
	and degree <code>degOfA</code> at
	<code>z = x + i y</code> and returns the result in <code>aAtZ</code>.
	The polynomial is evaluated using Horner\uFFFDs method.
	If the <code>offsetOfA</code> is negative the occuring negative array entries
	are assumed as zeros. */
    public static void eval( final double [] aRe,
			     final double [] aIm,
			     final int offsetOfA, final int degOfA,
			     final double x, final double y, final double [] aAtZ ) {

	if( degOfA < -offsetOfA ) {
	    aAtZ[0] = 0;
	    aAtZ[1] = 0;
	    return;
	}

	double resultRe = aRe[ offsetOfA + degOfA ];
	double resultIm = aIm[ offsetOfA + degOfA ];

	if( degOfA > 0 ) {

	    int i = degOfA-1;

	    for( int index = i+offsetOfA; index >= 0 && i >= 0 ; i--, index--) {

		final double bRe = resultRe;
		final double bIm = resultIm;

		resultRe = bRe * x - bIm * y + aRe[ index ];
		resultIm = bRe * y + bIm * x + aIm[ index ];
	    }

	    for( ; i >= 0 ; i-- ) {

		final double bRe = resultRe;
		final double bIm = resultIm;

		resultRe = bRe * x - bIm * y;
		resultIm = bRe * y + bIm * x;
	    }
	}

	aAtZ[0] = resultRe;
	aAtZ[1] = resultIm;
    }

    /** evaluates the <code>n</code>-the derivative of the polynomial
	<code>a</code> with the coefficients
	<code>a_k = aRe[k + offsetOfA ] + i aIm[k + offsetOfA]</code>
	and degree <code>degOfA</code> at
	<code>z = x + i y</code> and returns the result in <code>dAAtZ</code>.
	The polynomial is evaluated using Horner\uFFFDs method. */
    public static void evalDerivative( final double [] aRe,
				       final double [] aIm,
				       final int degOfA, final int n,
				       final double x, final double y, final double [] dAAtZ ) {
	evalDerivative( aRe, aIm, 0, degOfA, n, x, y, dAAtZ );
    }

    /** evaluates the <code>n</code>-the derivative of the polynomial
	<code>a</code> with the coefficients
	<code>a_k = aRe[k + offsetOfA ] + i aIm[k + offsetOfA]</code>
	and degree <code>degOfA</code> at
	<code>z = x + i y</code> and returns the result in <code>dAAtZ</code>.
	The polynomial is evaluated using Horner\uFFFDs method.
	If the <code>offsetOfA</code> is negative the occuring negative array entries
	are assumed as zeros. */
    public static void evalDerivative( final double [] aRe,
				       final double [] aIm,
				       final int offsetOfA, final int degOfA, final int n,
				       final double x, final double y, final double [] dAAtZ ) {

	if( degOfA < -offsetOfA || degOfA < n ) {
	    dAAtZ[0] = 0;
	    dAAtZ[1] = 0;
	    return;
	}

	double factor = 1;
	for( int i = degOfA; i > degOfA-n; i-- )
	    factor *= i;

	double resultRe = factor * aRe[offsetOfA + degOfA];
	double resultIm = factor * aIm[offsetOfA + degOfA];

	if( degOfA > 0 ) {

	    int e = degOfA-1;

	    for( int index = e + offsetOfA; index >= 0 && e >= n ; e--, index--) {

		factor *= e+1-n;
		factor /= e+1;

		final double bRe = resultRe;
		final double bIm = resultIm;

		resultRe = bRe * x - bIm * y + factor * aRe[ index ];
		resultIm = bRe * y + bIm * x + factor * aIm[ index ];
	    }

	    for( ; e >= n ; e-- ) {

		final double bRe = resultRe;
		final double bIm = resultIm;

		resultRe = bRe * x - bIm * y;
		resultIm = bRe * y + bIm * x;
	    }
	}

	dAAtZ[0] = resultRe;
	dAAtZ[1] = resultIm;
    }





    public static int calcDeg( final double[] aRe, final double[] aIm, int degOfA )
    {
	for( int i = degOfA; i >= 0; i--)
	    {
		if( aRe[ i ] * aRe[ i ] + aIm[ i ] * aIm[ i ] >= EPS * EPS ) return i;
	    }
	return -1;
    }

    /** This methods multiplies a complex polynomial with constant complex number.
    */
    public static void multiplyToConst( final double lr, final double li,
					final double[] aRe, final double[] aIm, int degOfA,
					final double[] cRe, final double[] cIm )
					// C = ( lr + i * li) *A
    {
	for( int i = 0; i <= degOfA; i++)
	    {
		double ar = aRe[i]; double ai = aIm[i];
		cRe[ i ] = lr * ar - li * ai;
		cIm[ i ] = lr * ai + li * ar;
	    }
    }

    public static void plus( final double[] aRe, final double[] aIm, int degOfA,
			     final double[] bRe, final double[] bIm, int degOfB,
			     final double[] cRe, final double[] cIm )
    {
	if( degOfA == degOfB )
	    {
		for(int i=0; i<= degOfA; i++)
		    {
			cRe[i] = aRe[i]+bRe[i];
			cIm[i] = aIm[i]+bIm[i];
		    }
	    }
	else
	    {
		if(degOfA < degOfB)
		    {
			for(int i=0; i<= degOfA; i++)
			    {
				cRe[i] = aRe[i]+bRe[i];
				cIm[i] = aIm[i]+bIm[i];
			    }
			System.arraycopy( bRe, degOfA+1, cRe, degOfA+1, degOfB-degOfA );
			System.arraycopy( bIm, degOfA+1, cIm, degOfA+1, degOfB-degOfA );
		    }
		else
		    {
			for(int i=0; i<= degOfB; i++)
			    {
				cRe[i] = aRe[i]+bRe[i];
				cIm[i] = aIm[i]+bIm[i];
			    }
			System.arraycopy( aRe, degOfB+1, cRe, degOfB+1, degOfA-degOfB );
			System.arraycopy( aIm, degOfB+1, cIm, degOfB+1, degOfA-degOfB );
		    }
	    }
    }

    public static void linComb( double lr, double li,
				final double[] aRe, final double[] aIm, int degOfA,
				double mr, double mi,
				final double[] bRe, final double[] bIm, int degOfB,
				final double[] cRe, final double[] cIm )
				// C= ( lr + i*li ) * A + ( mr + i*mi ) * B
    {
	double ar,ai,br,bi;
	if( degOfA == degOfB )
	    {
		for(int i=0; i<= degOfA; i++)
		    {
			ar = aRe[ i ]; ai = aIm[ i ]; br = bRe[ i ]; bi = bIm[ i ];
			cRe[i] = lr * ar - li * ai + mr * br - mi * bi ;
			cIm[i] = lr * ai + li * ar + mr * bi + mi * br ;
		    }
	    }
	else
	    {
		if(degOfA < degOfB)
		    {
			for(int i = 0; i <= degOfA; i++)
			    {
				ar = aRe[ i ]; ai = aIm[ i ]; br = bRe[ i ]; bi = bIm[ i ];
				cRe[i] = lr * ar - li * ai + mr * br - mi * bi ;
				cIm[i] = lr * ai + li * ar + mr * bi + mi * br ;
			    }
			for(int i = degOfA+1; i <= degOfB; i++)
			    {
				br = bRe[ i ]; bi = bIm[ i ];
				cRe[i] = mr * br - mi * bi ;
				cIm[i] = mr * bi + mi * br ;
			    }

		    }
		else
		    {
			for(int i=0; i<= degOfB; i++)
			    {
				ar = aRe[ i ]; ai = aIm[ i ]; br = bRe[ i ]; bi = bIm[ i ];
				cRe[i] = lr * ar - li * ai + mr * br - mi * bi ;
				cIm[i] = lr * ai + li * ar + mr * bi + mi * br ;
			    }
			for(int i=degOfB+1; i<= degOfA; i++)
			    {
				ar = aRe[ i ]; ai = aIm[ i ];
				cRe[i] = lr * ar - li * ai;
				cIm[i] = lr * ai + li * ar;
			    }
		    }
	    }
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
				   final int degOfA,
				   final int n,
				   final double [] bRe, final double [] bIm ) {
	int factor = 1;
	for( int i=2; i<=n; i++ )
	    factor *= i;

	for( int i=0, j=n; j<=degOfA; i++, j++ ) {
	    bRe[i] = aRe[j] * factor;
	    bIm[i] = aIm[j] * factor;

	    factor *= (j+1);
	    factor /= (i+1);
	}
    }

    public static void times( final double[] aRe, final double[] aIm, int degOfA,
			      final double[] bRe, final double[] bIm, int degOfB,
			      final double[] cRe, final double[] cIm ) {

	times( aRe, aIm, degOfA,
	       bRe, bIm, degOfB,
	       cRe, cIm, degOfA + degOfB );
    }

    public static void times( final double[] aRe, final double[] aIm, int degOfA,
			      final double[] bRe, final double[] bIm, int degOfB,
			      final double[] cRe, final double[] cIm, int degOfC ) {


	int k = degOfC;

	for( ; k > degOfA + degOfB; k-- ) {
	    cRe[ k ] = 0;
	    cIm[ k ] = 0;
	}

	for( ; k >= 0; k-- ) {
	    double cr = 0.0;
	    double ci = 0.0;

	    final int maxI = Math.min( k, degOfA );
	    for( int i = Math.max(0,k-degOfB); i <= maxI; i++ ) {

		cr = cr + aRe[ i ] * bRe[ k-i ] - aIm[ i ] * bIm[ k-i ];
		ci = ci + aIm[ i ] * bRe[ k-i ] + aRe[ i ] * bIm[ k-i ];
	    }

	    cRe[ k ] = cr;
	    cIm[ k ] = ci;
	}
    }

    public static void divide( final double[] aRe, final double[] aIm, int degOfA,
			       final double[] bRe, final double[] bIm, int degOfB,
			       final double[] cRe, final double[] cIm,
			       final double[] dRe, final double[] dIm  )
    {
	System.arraycopy( aRe, 0, dRe, 0, degOfA+1 );
	System.arraycopy( aIm, 0, dIm, 0, degOfA+1 );
	if( degOfA >= degOfB )
	    {
		double br=bRe[degOfB]; double bi=bIm[degOfB];
		double b1r = br/( br * br + bi * bi);
		double b1i = - bi/( br * br + bi * bi);
		for( int k = degOfA - degOfB ; k >= 0; k--)
		    {
			double b2r = b1r * dRe [ k + degOfB] -  b1i * dIm [ k + degOfB];
			double b2i = b1i * dRe [ k + degOfB] +  b1r * dIm [ k + degOfB];
			cRe[k] = b2r; cIm[k] = b2i;
			for( int i = degOfB; i >= 0; i--)
			    {
				dRe [ k + i] -=  (b2r * bRe [ i ] -  b2i * bIm [ i ]);
				dIm [ k + i] -=  (b2i * bRe [ i ] +  b2r * bIm [ i ]);
			    }
		    }
	    }
	else
	    {
		cRe[ 0 ] =0.0; cIm[ 0 ] =0.0;
	    }
    }

    public static void laguerOneRoot( final double[] Ar, final double[] Ai, int degOfA,
				      double calcError , final double[] root )
    {
	final int MR = 8;
	final int MT = 1000;
	final int MAXIT = MR * MT;
	double abx, abp, abm, err;
	double xr, xi, dxr, dxi, x1r, x1i, br, bi, dr, di, fr, fi, bufr, bufi , bs  ;
	double gr , gi , hr , hi , sqr , sqi, gpr , gpi , gmr , gmi , g2r , g2i ;

	double[] frac = {0.0,0.5,0.25,0.75,0.13,0.38,0.62,0.88,1.0};

	xr = root[ 0 ]; xi = root[ 1 ];
	for( int iter = 1; iter <= MAXIT; iter++ )
	    {
		br = Ar[degOfA]; bi = Ai[degOfA];
		err = Math.sqrt(br * br + bi * bi );
		dr = 0.0; di = 0.0; fr = 0.0; fi = 0.0;
		abx = Math.sqrt(xr * xr + xi * xi );
		for( int j = degOfA-1; j >= 0; j --)
		    {
			bufr = fr; bufi = fi;
			fr = xr * bufr - xi * bufi + dr;
			fi = xr * bufi + xi * bufr + di;
			bufr = dr; bufi = di;
			dr = xr * bufr - xi * bufi + br;
			di = xr * bufi + xi * bufr + bi;
			bufr = br; bufi = bi;
			br = xr * bufr - xi * bufi + Ar[ j ];
			bi = xr * bufi + xi * bufr + Ai[ j ];
			err = Math.sqrt( br * br + bi * bi ) + abx * err;
		    }
		bs = br * br + bi * bi;
		err *= calcError;

		if( bs <= (err * err) )
		    {
			root[ 0 ] = xr; root[ 1 ] = xi;
			return;
		    }
		gr = (dr * br + di * bi ) / bs; gi = (di * br - dr * bi ) / bs;
		g2r = (gr * gr - gi * gi ); g2i = 2. * gi * gr ;
		bufr = (fr * br + fi * bi ) / bs; bufi = (fi * br - fr * bi ) / bs;
		hr = g2r - 2. * bufr; hi = g2i - 2. * bufi;
		bufr = ( degOfA -1. ) * ( degOfA * hr - g2r );
		bufi = ( degOfA -1. ) * ( degOfA * hi - g2i );
		if( bufr == 0.0  && bufi == 0)
		    {
			sqr = 0.0;
			sqi = 0.0;
		    }
		else
		    {
			if ( bufr >= 0 )
			    {
				sqr = Math.sqrt( ( Math.sqrt(bufr * bufr + bufi * bufi )  + bufr )/2. );
				sqi = bufi / ( 2. * sqr );
			    }
			else
			    {
				sqi = Math.sqrt( ( Math.sqrt(bufr * bufr + bufi * bufi )  - bufr )/2. );
				sqr = bufi / ( 2. * sqi );
			    }
		    }

		gpr = gr + sqr; gpi = gi + sqi;
		gmr = gr - sqr; gmi = gi - sqi;
		abp = ( gpr * gpr + gpi * gpi );
		abm = ( gmr * gmr + gmi * gmi );
		if( abp < abm )
		    {
			gpr = gmr;
			gpi = gmi;
			abp = abm;
		    }
		if( abp > 0.0 )
		    {
			dxr =   gpr * degOfA / abp;
			dxi = - gpi * degOfA / abp;
		    }
		else
		    {
			dxr = (1. + abx) * Math.cos((double)iter );
			dxi = (1. + abx) * Math.sin((double)iter ) ;
		    }
		x1r = xr - dxr ; x1i = xi - dxi ;
		if( x1r == xr && x1i == xi )
		    {
			root[ 0 ] = xr; root[ 1 ] = xi;
			return;
		    }
		if( (iter % MT) != 0  )
		    {
			xr = x1r; xi = x1i;
		    }
		else
		    {
			xr = xr - dxr * frac[ iter / MT ] ;
			xi = xi - dxi * frac[ iter / MT ] ;
		    }
	    }
	throw new IllegalArgumentException(" The root was not found -- too many iterations " );
    }

    public static void laguerAllRoots( final double[] Ar, final double[] Ai, int degOfA,
				      double calcError , final double[] rootsRe , final double[] rootsIm )

    {
	double bufr , bufi  ;
	double[] root = new double[2];
	double[] adr  = new double[ degOfA + 1 ];
	double[] adi  = new double[ degOfA + 1 ];

	System.arraycopy( Ar, 0 , adr, 0 , degOfA+1 );
	System.arraycopy( Ai, 0 , adi, 0 , degOfA+1 );

	for( int j =  degOfA; j >= 1; j-- )
	    {
		root[ 0 ] = root [ 1 ] =0.0;
		laguerOneRoot( adr, adi, j , calcError , root  );

		double xr = rootsRe[ j - 1 ] = root[ 0 ];
		double xi = rootsIm[ j - 1 ] = root[ 1 ];
		double br = adr[ j ];
		double bi = adi[ j ];

		for ( int jj = j - 1; jj >= 0 ; jj-- )
		    {
			double cr = adr [ jj ];
			double ci = adi [ jj ];
			adr [ jj ] = br ;
			adi [ jj ] = bi;
			bufr = cr + xr * br  -  xi * bi;
			bufi = ci + xi * br  +  xr * bi;
			br = bufr ; bi = bufi;
		    }
	    }

	// The polishing procedure -- not activated
	// {
	//   for( int j = 0; j < degOfA; j++)
	//     {
	//	  root [ 0 ] = rootsRe[j] ; root [ 1 ] = rootsIm[j] ;
	//	  laguerOneRoot( Ar, Ai, degOfA , calcError , root );
	//	  rootsRe[j] = root [ 0 ] ; rootsIm[j] = root [ 1 ];
	//     }
	// }
    }


    public static void setZero(int deg, final double[] a)
    {
	for( int i=0; i<=deg; i++ )
	    a[i]=0.0;
    }

    public static String toString( final double [] re,
				   final double [] im, final int degree ) {
	return toString( re, im, 0, degree, EPS );
    }

    public static String toString( final double [] re,
				   final double [] im, final int degree, final double eps ) {
	return toString( re, im, 0, degree, eps );
    }

    public static String toString( final double [] re,
				   final double [] im,
				   final int lowerDegree, final int degree, final double eps ) {

	StringBuffer sb=new StringBuffer(300);

	boolean first = true;

	for( int i=lowerDegree; i<=degree; i++ ) {

	    final double x = re[i - lowerDegree ];
	    final double y = im[i - lowerDegree ];

	    if( x*x + y*y > EPS * EPS ) {

		if( !first )
		    sb.append( " + " );
		else
		    first = false;

		sb.append( "( " );
		sb.append( x );

		if( y < 0 ) {
		    sb.append( " - i" );
		    sb.append( (-y) );
		} else {
		    sb.append( " + i" );
		    sb.append(   y );
		}

		sb.append( " ) z^" );
		sb.append( i );
	    }
	}

	if( first )
	    sb.append( "0" );

	return sb.toString();
    }

}








