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

package de.jtem.mfc.polynomial;

import java.io.Serializable;

import de.jtem.mfc.field.Complex;
import de.jtem.mfc.field.Field;

/** This class implements Laurent polynomials in one variable with Field.Complex
    coefficients. */

public final class ComplexLaurentPolynomial implements Complex.FunctionOnComplex, Serializable {
    private static final long serialVersionUID = 1L;

    /** 1e-14; constant for zero tests and precision demands. */
    static final double EPS = 1.e-14;

    double[] re;
    double[] im;

    int degree;
    int lowerDegree;

    /** Constructs a monomial <pre>f(z) = a * z^n<\pre> with <pre>z = x + iy<\pre>*/
    public ComplexLaurentPolynomial( int n, double x, double y ) {
	re = new double[] { x };
	im = new double[] { y };

	lowerDegree = degree = n;
    }

    /** Constructs a monomial <pre>f(z) = a * z^n<\pre>*/
    public ComplexLaurentPolynomial( int n, Field.Complex a ) {
	this( n, a.getRe(), a.getIm() );
    }

    /** Constructs a monomial <pre>f(z) = z^n<\pre>*/
    public ComplexLaurentPolynomial( int n ) {
	this( n, 1, 0 );
    }

    /** Constructs the zero Laurent polynomial. */
    public ComplexLaurentPolynomial() {
	this( 0, 0, 0 );
    }

    /** Constructs Laurent polynomial with prescribed degrees, but
	with all coefficients equal to zero. */
    public ComplexLaurentPolynomial( int lowerDegree, int degree ) {
	this.degree      = degree;
	this.lowerDegree = lowerDegree;

	re = new double[ degree - lowerDegree + 1 ];
	im = new double[ degree - lowerDegree + 1 ];
    }

    /** Creates copy of <pre>p<\pre>. */
    public ComplexLaurentPolynomial( ComplexPolynomial p ) {
	re = (double [])p.re.clone();
	im = (double [])p.im.clone();
	degree = p.degree;
	lowerDegree = 0;
    }

    /** Creates copy of <pre>p<\pre>. */
    public ComplexLaurentPolynomial( ComplexLaurentPolynomial p ) {
	re = (double [])p.re.clone();
	im = (double [])p.im.clone();
	degree = p.degree;
	lowerDegree = p.lowerDegree;
    }

    /** Constructs a constant Laurent polynomial with value <pre>a<\pre> */

    public ComplexLaurentPolynomial(Field.Complex a) {
	this( 0, a.getRe(), a.getIm() );
    }

    /** Constructs a Laurent polynomial with lower degree <pre>lowerDeg<\pre>
	whose coefficients are copied from the array <pre>a<\pre> of complex
	numbers. */
    public ComplexLaurentPolynomial(int lowerDeg, Complex[] a) {
	assign(lowerDeg, a);
    }

    /** Constructs a polynomial with (<pre>lowerDeg = 0<\pre>)
	whose coefficients are copied from the array <pre>a<\pre> of complex
	numbers. */

    public ComplexLaurentPolynomial(Complex[] a) {
	this(0, a);
    }

    /** Constructs a Laurent polynomial with lower degree <pre>lowerDeg<\pre>
	described by two double arrays <pre>aRe<\pre> and <pre>aIm<\pre> which
	are stored (reference) as the real and imaginary parts of the coefficients.
	The length of <pre>aRe<\pre> and <pre>aIm<\pre> has to coincide. */
    public ComplexLaurentPolynomial(int lowDeg, double[] aRe, double[] aIm ) {

	if( aRe.length != aIm.length )
	    throw new IllegalArgumentException( "array lengths do not coincide" );

	lowerDegree = lowDeg;
	degree      = aRe.length + lowDeg - 1;

	re = aRe;
	im = aIm;
    }

    /** Returns reference to the array storing the real part of coefficents. */
    public double [] re() {
	return re;
    }

    /** Returns reference to the array storing the imaginary part of coefficents. */
    public double [] im() {
	return im;
    }

    /** Returns the lower degree of the Laurent polynomial in the sense that
	the coefficients of order less than the lower degree are exactly zero
	by definition. Note that all other coefficients could be zero
	<em>numerically<\em>. */

    public int getLowerDegree() {
	return lowerDegree;
    }

    /** Sets the lower degree to <pre>lowDeg<\pre>. Adds coefficients which are
	numerically zero if lowDeg is smaller than the current degree,
	truncates the Laurent polynomials from below to degree otherwise. */

    public void setLowerDegree(int lowDeg) {
	setDegrees( lowDeg, degree );
    }

    /** Truncates the Laurent from Below to lower degree <pre>lowDeg<\pre>.
	Does nothing if <pre>lowDeg<\pre> is smaller than the current
	degree. */

    public void truncateBelow(int lowDeg) {
	if ( lowDeg > lowerDegree ) {
	    setLowerDegree(lowDeg);
	}
    }

    /** Returns the degree of the Laurent polynomial in the sense that the
	coefficients of order bigger than the degree are exactly zero by
	definition. Note that for a Laurent polynomial say of degree seven all
	coefficients could be zero <em>numerically<\em>. */

    public int getDegree() {
	return degree;
    }

    /** Sets the degree to <pre>newDeg<\pre>. Adds coefficients which are
	numerically zero if newDeg is bigger than the current degree, truncates
	the Laurent polynomials to degree <pre>newDeg<\pre> otherwise. */

    public void setDegree(int newDeg) {
	setDegrees( lowerDegree, newDeg );
    }

    /** Sets the lower degree to <pre>lowDeg<\pre> and the degree to <pre>newDeg<\pre>.
	Adds coefficients which are
	numerically zero if new degrees exeed the current degrees, truncates
	the Laurent polynomials otherwise. */
    
    public void setDegrees(int lowDeg, int newDeg ) {
    	if (newDeg != degree || lowDeg != lowerDegree ) {
    		
    		final int newLength = newDeg - lowDeg + 1;
    		
    		final double[] newRe = new double[newLength];
    		final double[] newIm = new double[newLength];
    		
    		if( re != null ) {
    			
    			final int offset = lowDeg - lowerDegree;
    			
    			final int polyStart = Math.max(0,  offset);
    			final int newStart  = Math.max(0, -offset);
    			
    			final int len = Math.min( degree, newDeg ) - lowerDegree + 1 - polyStart;
    			
    			if( len > 0 ) {
    				System.arraycopy(re, polyStart, newRe, newStart, len);
    				System.arraycopy(im, polyStart, newIm, newStart, len);
    			}
    		}
    		
    		re = newRe;
    		im = newIm;
    		
    		lowerDegree = lowDeg;
    		degree      = newDeg;
    	}
    }

    /** Truncates the Laurent polynomial to degree <pre>newDeg<\pre>. Does
	nothing if <pre>newDeg<\pre> is bigger than the current degree. */

    public void truncate(int deg) {
	if ( deg < degree ) {
	    setDegree(deg);
	}
    }


    /** Truncates off all low and high coefficients with an absolute value
     ** smaller than eps.
     **/
    public void truncate( final double eps ) {
	final double sqrEps = eps * eps;
	int i, j, k;
	for( i=degree, k=re.length-1; i>=lowerDegree; i--, k-- )
	    if( re[k]*re[k] + im[k]*im[k] > sqrEps )
		break;

	if( i < lowerDegree ) {
	    assignMonomial( 0, 0, 0 );
	    return;
	}

	for( j=lowerDegree, k=0; j<=i; j++, k++ )
	    if( re[k]*re[k] + im[k]*im[k] > sqrEps ) {
		setDegrees( j, i );
		return;
	    }

	setDegrees( i, i );
    }

    /** Truncates off all low and high coefficients with an absolute value
     ** smaller than {#EPS}. */
    public void truncate() {
	truncate( EPS );
    }

    /** Assigns this to zero polynomial */
    public void assignZero() {
	assignMonomial( 0, 0, 0 );
    }

    /** Assigns this to monimial <pre>f(z) = a * z^n<\pre> with <pre>z = x + iy<\pre>. */
    public void assignMonomial( int n, double x, double y ) {
	if( degree != lowerDegree ) {
	    re = new double[] { x };
	    im = new double[] { y };
	} else {
	    re[0] = x;
	    im[0] = y;
	}

	lowerDegree = degree = n;
    }

    /** Assigns this to monimial <pre>f(z) = a * z^n<\pre>. */
    public void assignMonomial( int n, Field.Complex z ) {
	assignMonomial( n, z.getRe(), z.getIm() );
    }

    /** Assigns this to monimial <pre>f(z) = z^n<\pre>. */
    public void assignMonomial( int n ) {
	assignMonomial( n, 1, 0 );
    }

    /** Assign this with laurent polynomial <pre>p<\pre>. */
    public void assign( ComplexLaurentPolynomial p ) {
	assign( p.lowerDegree, p.re, p.im );
    }

    /** Reads the coefficients of this Laurent polynomial from the array of
	complex numbers <pre>a<\pre> an sets the lower degree to
	<pre>lowDeg>\pre> */

    public void assign( int lowDeg, Complex[] a ) {
	if ( a.length != degree - lowerDegree + 1 ) {
	    re = new double[a.length];
	    im = new double[a.length];
	}
	for ( int i=0; i< a.length; i++ ) {
	    re[i] = a[i].re;
	    im[i] = a[i].im;
	}
	degree      = lowDeg + a.length - 1;
	lowerDegree = lowDeg;
    }

    /** Sets the lower degree to lowDeg and reads the real and imaginary
	parts of the coefficients from the double arrays <pre>aRe<\pre> and
	<pre>aIm<\pre>. aRe and aIm are assumed to have the same length. */

    public void assign( int lowDeg, double[] aRe, double[] aIm ) {
	if ( aRe.length != degree - lowerDegree + 1 ) {
	    re = (double[])aRe.clone();
	    im = (double[])aIm.clone();
	} else {
	    System.arraycopy(aRe, 0, re, 0, degree - lowerDegree + 1);
	    System.arraycopy(aIm, 0, im, 0, degree - lowerDegree + 1);
	}
	degree = lowDeg + aRe.length - 1;
	lowerDegree = lowDeg;
    }

    /** Writes the coefficients starting with the coefficient of order
	startOrder to the array <pre>a<\pre> of complex numbers (as far as
	they fit in). Creates new instances of Field.Complex only as needed. */

    public void coefficients( int startOrder, Complex[] a ) {
	int minLength = Math.min(degree - startOrder + 1, a.length);
	int j;
	double aRe, aIm;

	for (int i=0; i< minLength; i++) {
	    j = i + startOrder - lowerDegree;
	    if ( a[i] == null ) {
		a[i] = new Complex(re[j], im[j]);
	    } else {
		a[i].assign(re[j], im[j]);
	    }
	}
    }

    /** Returns the coefficients starting  with the coefficient of order
	startOrder as a freshly created array of complex
	numbers. */
    public Complex[] coefficients(int startOrder) {
	Complex[] c = new Complex[degree + 1 - startOrder];
        coefficients(startOrder, c);
	return c;
    }

    /** Returns all coefficients as a freshly created array of complex
	numbers. */
    public Complex[] getCoefficients() {
	return coefficients(lowerDegree);
    }

    /** Returns the <pre>i<\pre>'th coefficient as a freshly created Field.Complex
	number. Returns zero if <pre>i<\pre> is smaller than the lower Degree
	or larger than the degree. */
    public Complex getCoefficient(int i) {
	Complex z = new Complex();
	getCoefficient(i,z);
	return z;
    }

    /** Writes the value of the <pre>i<\pre>'th coefficient (zero if
	<pre>i<\pre> is smaller than the lower Degree or larger than the
	degree) to the complex number <pre>z<\pre>. */
    public void getCoefficient(int i, Complex z) {
	if ( i >= lowerDegree && i <= degree ) {
	    z.assign( re[i - lowerDegree], im[i - lowerDegree] );
	} else {
	    z.assign(0,0);
	}
    }

    /** Reads the value of the <pre>i<\pre>'th coefficient from the complex
	number <pre>z<\pre>. <pre>i<\pre> is assumed to be non-negative and
	not bigger than the degree. */
    public void setCoefficient(int i, Field.Complex z) {
	setCoefficient( i, z.getRe(), z.getIm() );
    }

    /** Reads the value of the <pre>i<\pre>'th coefficient from the complex
	number <pre>x+iy<\pre>. <pre>i<\pre> is assumed to be non-negative and
	not bigger than the degree. */

    public void setCoefficient(int i, double x, double y) {
	if( i<lowerDegree )
	    setLowerDegree( i );
	if( i>     degree )
	    setDegree( i );

	re[i - lowerDegree] = x;
	im[i - lowerDegree] = y;
    }

    /** sets all coefficients to these in <pre>p<\pre> if they occur and to zero if not. */
    public void assignCoefficients( ComplexLaurentPolynomial p ) {

	final int minDegree      = Math.min(      degree, p.     degree );
	final int maxLowerDegree = Math.max( lowerDegree, p.lowerDegree );

	final int length = minDegree - maxLowerDegree + 1;

	System.arraycopy( p.re, maxLowerDegree - p.lowerDegree,
			  re, maxLowerDegree - lowerDegree, length );
	System.arraycopy( p.im, maxLowerDegree - p.lowerDegree,
			  im, maxLowerDegree - lowerDegree, length );

	for( int i = lowerDegree; i<maxLowerDegree; i++ ) {
	    re[ i - lowerDegree ] = 0;
	    im[ i - lowerDegree ] = 0;
	}

	for( int i = minDegree + 1; i <= degree; i++ ) {
	    re[ i - lowerDegree ] = 0;
	    im[ i - lowerDegree ] = 0;
	}
    }

    /** Returns the value of the polynomial at <pre>z<\pre> as a freshly
	created Field.Complex number. */
    public Complex eval(Complex z) {
	Complex w = new Complex();
	eval(z, w);
	return w;
    }

    private static double[] p = new double[2];

    /** Writes the value of the polynomial at <pre>z<\pre> to the Field.Complex
	number <pre>w<\pre>. */
    public void eval(Complex z, Complex w) {
	de.jtem.numericalMethods.algebra.polynomial.ComplexLaurentPolynomial.eval
	    ( re, im, lowerDegree, degree, z.re, z.im, p );
	w.assign(p[0],p[1]);
    }

    /** Returns the value of the n-the derivative of the polynomial
	at <pre>z<\pre> as a freshly created Field.Complex number. */
    public Complex evalDerivative(Complex z, int n) {
	Complex w = new Complex();
	evalDerivative(z, n, w );
	return w;
    }

    /** Writes the value of the polynomial at <pre>z<\pre> to the Field.Complex
	number <pre>w<\pre>. */
    public void evalDerivative(Complex z, int n, Complex w) {
	double[] p = new double[2];
	de.jtem.numericalMethods.algebra.polynomial.ComplexLaurentPolynomial.evalDerivative
	    ( re, im, lowerDegree, degree, n, z.re, z.im, p );
	w.assign( p[0],p[1] );
    }

    /** Returns the sum of this polynomial and the polynomial <pre>q<\pre> as
	 a freshly created ComplexLaurentPolynomial. */

    public ComplexLaurentPolynomial plus(ComplexLaurentPolynomial q) {
	ComplexLaurentPolynomial result = new ComplexLaurentPolynomial();
	result.assignPlus(this,q);
	return result;
    }

    /** Makes this polynomial equal to the sum of itself
	and the polynomial <pre>q<\pre>. */
    public void assignPlus(ComplexLaurentPolynomial q) {
	assignPlus(this, q);
    }

    /** Makes this Laurent polynomial equal to the sum of the Laurent
	polynomials <pre>q<\pre> and <pre>q<\pre>. */
    public void assignPlus( ComplexLaurentPolynomial p,
			    ComplexLaurentPolynomial q ) {
	int      newDegree = Math.max( p.     degree, q.     degree );
	int newLowerDegree = Math.min( p.lowerDegree, q.lowerDegree );

	double[] newRe = new double[newDegree - newLowerDegree + 1];
	double[] newIm = new double[newDegree - newLowerDegree + 1];

	for ( int i = p.lowerDegree; i <= p.degree; i++ ) {
	    newRe[i - newLowerDegree] = p.re[i - p.lowerDegree];
	    newIm[i - newLowerDegree] = p.im[i - p.lowerDegree];
	}
	for ( int i = q.lowerDegree; i <= q.degree; i++ ) {
	    newRe[i - newLowerDegree] += q.re[i - q.lowerDegree];
	    newIm[i - newLowerDegree] += q.im[i - q.lowerDegree];
	}

	re = newRe;
	im = newIm;

	degree      = newDegree;
	lowerDegree = newLowerDegree;

    }

    /** Returns the difference of this polynomial and the polynomial
	<pre>q<\pre> as a freshly created ComplexPolynomial. */
    public ComplexLaurentPolynomial minus(ComplexLaurentPolynomial q) {
	ComplexLaurentPolynomial result = new ComplexLaurentPolynomial();
	result.assignMinus(this,q);
	return result;
    }

    /** Makes this polynomial equal to the difference of itself and the
	Laurent Polynomial <pre>q<\pre>. */
    public void assignMinus(ComplexLaurentPolynomial q) {
	assignMinus(this, q);
    }

    /** Makes this polynomial equal to the difference of the polynomials
	<pre>q<\pre> and <pre>q<\pre>. */
    public void assignMinus( ComplexLaurentPolynomial p,
			     ComplexLaurentPolynomial q ) {
	int      newDegree = Math.max( p.     degree, q.     degree );
	int newLowerDegree = Math.min( p.lowerDegree, q.lowerDegree );

	double[] newRe = new double[newDegree - newLowerDegree + 1];
	double[] newIm = new double[newDegree - newLowerDegree + 1];

	for ( int i = p.lowerDegree; i <= p.degree; i++ ) {
	    newRe[i - newLowerDegree] = p.re[i - p.lowerDegree];
	    newIm[i - newLowerDegree] = p.im[i - p.lowerDegree];
	}
	for ( int i = q.lowerDegree; i <= q.degree; i++ ) {
	    newRe[i - newLowerDegree] -= q.re[i - q.lowerDegree];
	    newIm[i - newLowerDegree] -= q.im[i - q.lowerDegree];
	}

	re = newRe;
	im = newIm;

	degree      = newDegree;
	lowerDegree = newLowerDegree;

    }

    /** Returns the product of this Laurent polynomial and the Laurent
	polynomial <pre>q<\pre> as a freshly created Laurent polynomial. */
    public ComplexLaurentPolynomial times(ComplexLaurentPolynomial q) {
	ComplexLaurentPolynomial result = new ComplexLaurentPolynomial();
	result.assignTimes(this,q);
	return result;
    }

    /** Makes this Laurent polynomial equal to the product of itself
	and the Laurent polynomial <pre>q<\pre>. */
    public void assignTimes(ComplexLaurentPolynomial q) {
	assignTimes(this, q);
    }

    /** Makes this Laurent polynomial equal to the product of the Laurent
	polynomials <pre>p<\pre> and <pre>q<\pre>. */
    public void assignTimes( ComplexLaurentPolynomial p,
			     ComplexLaurentPolynomial q
			     ) {
	int newDegree      = p.     degree + q.     degree;
	int newLowerDegree = p.lowerDegree + q.lowerDegree;

	double[] newRe = new double[newDegree - newLowerDegree + 1 ];
	double[] newIm = new double[newDegree - newLowerDegree + 1 ];

	de.jtem.numericalMethods.algebra.polynomial.ComplexPolynomial.times
	    ( p.re, p.im, p.degree - p.lowerDegree,
	      q.re, q.im, q.degree - q.lowerDegree, newRe, newIm );

	re = newRe;
	im = newIm;

	degree      =      newDegree;
	lowerDegree = newLowerDegree;
    }

    /** Makes this Laurent polynomial equal to the product of the Laurent
	polynomial <pre>p<\pre> and the complex number <pre>z<\pre>. */
    public void assignTimes( ComplexLaurentPolynomial p, Complex z) {
	setNewDegrees(p);
	for ( int i=0; i<re.length; i++ ) {
	    final double pRe = p.re[i];
	    final double pIm = p.im[i];
	    re[i] = pRe * z.re - pIm * z.im;
	    im[i] = pRe * z.im + pIm * z.re;
	}
    }

    /** Makes this Laurent polynomial equal to the product of the Laurent
	polynomial <pre>p<\pre> and the complex number <pre>z<\pre>. */
    public void assignTimes(Complex z, ComplexLaurentPolynomial p) {
	assignTimes(p, z);
    }

    /** Makes this Laurent polynomial equal to the product of itself
	and the complex number <pre>z<\pre>. */
    public void assignTimes(Complex z) {
	assignTimes(this, z);
    }

    /** Returns the product of this Laurent polynomial and the complex
	number <pre>z<\pre> as a freshly created ComplexLaurentPolynomial. */
    public ComplexLaurentPolynomial times(Complex z) {
	ComplexLaurentPolynomial result = new ComplexLaurentPolynomial();
	result.assignTimes(this, z);
	return result;
    }

    /** Makes this polynomial equal to the product of the polynomial
	<pre>p<\pre> and the real number <pre>x<\pre>. */
    public void assignTimes(ComplexLaurentPolynomial p, double x) {
    	setNewDegrees(p);
    
	for ( int i=0; i<re.length; i++ ) {
	    re[i] = p.re[i] * x;
	    im[i] = p.im[i] * x;
	}
    }

    /** Makes this polynomial equal to the product of the polynomial
	<pre>p<\pre> and the real number <pre>x<\pre>. */
    public void assignTimes(double x, ComplexLaurentPolynomial p) {
	assignTimes(p, x);
    }

    /** Makes this polynomial equal to the product of itself
	and the real number <pre>z<\pre>. */
    public void assignTimes(double x) {
	assignTimes(this, x);
    }

    /** Returns the product of this Laurent polynomial and the real
	number <pre>x<\pre> as a freshly created ComplexLaurentPolynomial. */
    public ComplexLaurentPolynomial times(double x) {
	ComplexLaurentPolynomial result = new ComplexLaurentPolynomial();
	result.assignTimes(this, x);
	return result;
    }
    
    /**
     * Sets new degrees and allocates new arrays (if needed) without
     * copying any data, e.g. old coefficients are lost.
     */
    void setNewDegrees( int lowerDegree, int degree ) {
    	this.degree = degree;
    	this.lowerDegree = lowerDegree;
    	
    	if( re.length != degree - lowerDegree +1 ) {
    		re = new double[degree - lowerDegree +1];
    		im = new double[degree - lowerDegree +1];
    	};
    }
    
    /**
     * Sets new degrees and allocates new arrays (if needed) without
     * copying any data, e.g. old coefficients are lost.
     */
    void setNewDegrees(ComplexLaurentPolynomial p) {
    	setNewDegrees( p.lowerDegree, p.degree);
    }
    
    /** Makes this polynomial equal to the negative of the polynomial
     <pre>p<\pre> . */
    public void assignNeg(ComplexLaurentPolynomial p) {
    	setNewDegrees(p);

    	for ( int i=0; i< re.length; i++ ) {
    		re[i] = -p.re[i];
    		im[i] = -p.im[i];
    	}
    }

    /** Makes this polynomial equal to the product of itself
     and the real number <pre>z<\pre>. */
    public void assignNeg() {
    	assignNeg(this);
    }

    /** Returns the negative of this Laurent polynomial 
     * as a freshly created ComplexLaurentPolynomial. */
    public ComplexLaurentPolynomial neg() {
    	ComplexLaurentPolynomial result = 
    		new ComplexLaurentPolynomial(lowerDegree, degree );
    	result.assignNeg(this);
    	return result;
    }
    
    /** Assigns this Laurent polynomial with the star of p.
     * Starring a polynomial will conjugate the coefficients
     * and replaces the complex variable lambda by its
     * inverse.
     */
    public void assignStar(ComplexLaurentPolynomial p) {
  
    	if( p == this ) {
    		assignStar();
    		return;
    	}
    		
    	setNewDegrees( -p.degree,  -p.lowerDegree );
    
    	for( int i=0, j=re.length-1; i<re.length; i++, j-- ) {
    		re[i] = p.re[j];
    		im[i] = -p.im[j];
    	}
    }
    
    /**
     * Returns the star of this Laurent polynomrial.
     * Starring a polynomial will conjugate the coefficients
     * and replaces the complex variable lambda by its
     * inverse.
     */
    public ComplexLaurentPolynomial star() {
    	ComplexLaurentPolynomial result = new ComplexLaurentPolynomial(this);
    	result.assignStar();
    	return result;
    }
    
    /**
     * Starring a polynomial will conjugate the coefficients
     * and replaces the complex variable lambda by its
     * inverse.
     */
    public void assignStar() {
    	final int tmpDegree =degree;
    	degree = -lowerDegree;
    	lowerDegree = -tmpDegree;
    	
    	for( int i=0, j=re.length-1; i<j; i++, j--) {
    		final double tmpRe = re[i]; 
    		re[i] = re[j]; 
    		re[j] = tmpRe;
    		final double tmpIm = im[i]; 
    		im[i] = -im[j];
    		im[j] = -tmpIm;
    	}
    	
    	if( re.length % 2 == 1)
    		im[re.length/2] = -im[re.length/2];
    }
    
    /** assigns this with <code>n</code>th derivative of <code>p</code>.
	<code>p</code> and <code>this</code> my coincide. */
    public void assignDerivative( final ComplexLaurentPolynomial p, final int n ) {

	if( n<0 )
	    throw new IllegalArgumentException( "order of derivative is negative" );

	final double [] pRe = p.re;
	final double [] pIm = p.im;

	if( pRe.length != re.length ) {

	    re = new double[ degree - lowerDegree + 1];
	    im = new double[ degree - lowerDegree + 1];
	}

	de.jtem.numericalMethods.algebra.polynomial.ComplexLaurentPolynomial.derivative
	    ( pRe, pIm, p.lowerDegree, p.degree, n, re, im );

	degree      = p.degree      - n;
	lowerDegree = p.lowerDegree - n;
    }

    /** assigns this with its <code>n</code>th derivative. */
    public void assignDerivative( final int n ) {
	assignDerivative( this, n );
    }

    /** returns <code>n</code>th derivative of this. */
    public ComplexLaurentPolynomial derivative( final int n ) {
	ComplexLaurentPolynomial result
	    = new ComplexLaurentPolynomial
		( lowerDegree - n,
		  new double[ degree - lowerDegree + 1 ],
		  new double[ degree - lowerDegree + 1 ] );

	result.assignDerivative( this, n );

	return result;
    }

    /** assigns this with first derivative of <code>p</code>.
	<code>p</code> and <code>this</code> my coincide. */
    public void assignDerivative( final ComplexLaurentPolynomial p ) {
	this.assignDerivative( p, 1 );
    }

    /** assigns this with its first derivative. */
    public void assignDerivative() {
	assignDerivative( this, 1 );
    }

    /** returns first derivative of this. */
    public ComplexLaurentPolynomial derivative() {
	return derivative( 1 );
    }

    /** Returns the roots of the polynomial obtained by multiplying this
     Laurent Polynomial with the monomial <pre>f(x) = x^-lowerDeg<\pre>. *
    public Field.Complex[] getRoots() {
	double[] rootsRe = new double[degree];
	double[] rootsIm = new double[degree];
	Field.Complex[] roots = new Field.Complex[degree];
	de.jtem.numericalMethods.algebra.polynomial.ComplexPolynomial.laguerAllRoots
	    (re, im, degree, ROOT_PRECISION, rootsRe, rootsIm );
	for (int i=0; i<degree; i++) {
	    roots[i] = new Field.Complex(rootsRe[i], rootsIm[i]);
	}
	return roots;
    }

    /** Makes this a Laurent polynomial which, when multiplied with the
	monomial <pre>f(x) = x^-lowerDeg<\pre> becomes a polynomial having
	the prescribed roots. If the degree of this polynomial is
	equal to the number of roots and the highest coefficient is non-zero,
	the highest coefficient bis left unchanged. Otherwise the highest
	coefficient is set to one. */
    public void setRoots(Complex[] roots) {
	Field.Complex aN;
	if (roots.length == degree) {
	    aN = new Complex(re[degree],im[degree]);
	} else {
	    aN = new Complex( 1, 0 );
	}
	assignMonomial(0, aN);
	ComplexLaurentPolynomial q =
	    new ComplexLaurentPolynomial(1, 1, 0 );
	for (int i=0; i<roots.length; i++) {
	    q.setCoefficient(0, roots[i].neg());
	    assignTimes(q);
	}
    }

    /** Assigns all coefficients of this with random numbers. */
    public void assignRandom() {
	for( int i=lowerDegree; i<=degree; i++ ) {
	    setCoefficient( i, Math.random(), Math.random() );
	}
    }

    final private boolean isZero( int index, final double epsSqr ) {
	index -= lowerDegree;

	return re[index]*re[index] + im[index]*im[index] < epsSqr;
    }

    final private boolean equals( int index, ComplexLaurentPolynomial p, double epsSqr ) {
	final int pIndex = index - p.lowerDegree;

	index -= lowerDegree;

	final double dx = re[index]-p.re[pIndex];
	final double dy = im[index]-p.im[pIndex];

	return dx*dx + dy*dy < epsSqr;
    }

    /** test if coefficients of this and <pre>p<\pre> equals up to a precision of {#EPS}. */
    public boolean equals( ComplexLaurentPolynomial p ) {
	return equals( p, EPS );
    }

    /** test if coefficients of this and <pre>p<\pre> equals up to a precision of <pre>eps<\pre>. */
    public boolean equals( ComplexLaurentPolynomial p, double eps ) {
	final double epsSqr = eps*eps;
	if( p == this )
	    return true;

	if( lowerDegree < p.lowerDegree ) {
	    for( int i=lowerDegree; i<p.lowerDegree; i++ )
		if( !isZero( i, epsSqr ) )
		    return false;
	} else {
	    for( int i=p.lowerDegree; i<lowerDegree; i++ )
		if( !p.isZero( i, epsSqr ) )
		    return false;
	}

	if( degree > p.degree ) {
	    for( int i=p.degree+1; i<=degree; i++ )
		if( !isZero( i, epsSqr ) )
		    return false;
	} else {
	    for( int i=degree+1; i<=p.degree; i++ )
		if( !p.isZero( i, epsSqr ) )
		    return false;
	}

	final int min = Math.max( p.lowerDegree, lowerDegree );
	final int max = Math.min( p.     degree,      degree );

	for( int i=min; i<=max; i++ )
	    if( !equals( i, p, epsSqr ) )
		return false;

	return true;
    }

    /** test if coefficients of this and <pre>o<\pre> equals up to a precision of {#EPS}. */
    public boolean equals( Object o ) {
	if( o == this )
	    return true;

	try {
	    return equals( (ComplexLaurentPolynomial)o );
	} catch ( ClassCastException e ) {
	    return false;
	}
    }

    public String toString() {
	return 	de.jtem.numericalMethods.algebra.polynomial.ComplexLaurentPolynomial.toString
	    ( re, im, lowerDegree, degree, EPS );
    }


}


