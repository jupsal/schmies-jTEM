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
import de.jtem.mfc.field.ComplexConstant;
import de.jtem.mfc.field.Field;

/** This class implements polynomials in one variable with Field.Complex
    coefficients. Uses some algorithms like root finding and polynomial
    division from the class
    <pre>de.jtem.numericalMethods.algebra.polynomials.ComplexPolynomial<\pre> */

public class ComplexPolynomial implements Complex.FunctionOnComplex, Serializable {
    private static final long serialVersionUID = 1L;

    /** 1e-14; constant for zero tests and precision demands. */
    static final double EPS = 1.e-14;

    double [] re;
    double [] im;

    public int degree;

    /** Returns reference to the array storing the real part of coefficents. */
    public double [] re() {
	return re;
    }

    /** Returns reference to the array storing the imaginary part of coefficents. */
    public double [] im() {
	return im;
    }


    /** The default constructor returns the zero polynomial. */
    public ComplexPolynomial() {
	re = new double[1];
	im = new double[1];
    }

    /** Creates copy of <pre>p<\pre>. */
    public ComplexPolynomial( ComplexPolynomial p ) {
	re =(double []) p.re.clone();
	im =(double []) p.im.clone();
	degree = p.degree;
    }

    /** Constructs a constant polynomial with value <pre>a<\pre> */
    public ComplexPolynomial(Field.Complex a) {
	this(0, a);
    }

    /** Constructs a monomial <pre>f(z) = a * z^n<\pre>. */
    public ComplexPolynomial( final int deg, final Field.Complex a ) {
	this( deg, a.getRe(), a.getIm() );
    }

    /** Constructs a monomial <pre>f(z) = a * z^n<\pre> with <pre>z = x + iy<\pre>. */
    public ComplexPolynomial( final int deg, final double x, final double y ) {
	degree = deg;
	re = new double[degree +1];
	im = new double[degree +1];
	re[degree] = x;
	im[degree] = y;
    }

    /** Constructs a monomial <pre>f(z) = z^n<\pre>. */
    public ComplexPolynomial( final int deg ) {
	this( deg, 1, 0 );
    }

    /** Constructs a polynomial whose coefficients are copied from the array
	<pre>a<\pre> of complex numbers. */
    public ComplexPolynomial(Complex[] a) {
	assign(a);
    }


    /** Constructs a polynomial described by two double arrays <pre>aRe<\pre>
	and <pre>aIm<\pre> which are copied to the real and
	imaginary parts of the coefficients. */
    public ComplexPolynomial(double[] aRe, double[] aIm) {
	assign(aRe, aIm);
    }

    /** Returns the degree of the polynomial in the sense that the coefficients
	of order bigger than the degree are exactly zero by definition. Note
	that for a polynomial say of degree seven all coefficients could be
	zero <em>numerically<\em>. */

    public int getDegree() {
	return degree;
    }

    /** Sets the degree to <pre>newDeg<\pre>. Adds coefficients which are
	numerically zero if newDeg is bigger than the current degree, truncates
	the polynomials to degree <pre>newDeg<\pre> otherwise. */

    public void setDegree(int newDeg) {
	if (newDeg != degree) {

	    double[] newRe = new double[newDeg + 1];
	    double[] newIm = new double[newDeg + 1];
	    int minDeg = Math.min(degree,newDeg);
	    System.arraycopy( re, 0, newRe, 0, minDeg + 1 );
	    System.arraycopy( im, 0, newIm, 0, minDeg + 1 );
	    degree = newDeg;
	    re = newRe;
	    im = newIm;
	}
    }

    /** Truncates the polynomials to degree <pre>newDeg<\pre>. Does nothing if
	<pre>newDeg<\pre> is bigger than the current degree. */
    public void truncate(int deg) {
	if (deg < degree) {
	    setDegree(deg);
	}
    }

    /** Truncates off all high coefficients with an absolute value
     ** smaller than eps. */
    public void truncate( final double eps ) {
	final double sqrEps = eps * eps;
	int i;
	for( i=degree; i>0; i-- )
	    if( re[i]*re[i] + im[i]*im[i] > sqrEps )
		break;
	setDegree( i );
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

    /** Makes this polynomial equal to the monomial <pre>f(z) = a * z^n<\pre>. */
    public void assignMonomial( final int deg, final Field.Complex a ) {
	assignMonomial( deg, a.getRe(), a.getIm() );
    }

    /** Assigns this to monimial <pre>f(z) = a * z^n<\pre> with <pre>z = x + iy<\pre>. */
    public void assignMonomial( final int deg, final double x, final double y ) {
	if ( degree != deg) {
	    degree = deg;
	    re = new double[degree + 1];
	    im = new double[degree + 1];
	} else {
	    for( int i=0; i<degree; i++ )
		re[i] = im[i] = 0;
	}
	re[degree] = x;
	im[degree] = y;
    }

    /** Makes this polynomial equal to the monomial <pre>f(z) = z^n<\pre>. */

    public void assignMonomial(int deg) {
	assignMonomial(deg, 1, 0 );
    }

    /** Reads the coefficients of this polynomial from the array of
	complex numbers <pre>a<\pre> */
    public void assign(Complex[] a) {
	if ( degree != a.length - 1 || re == null ) {
	    degree = a.length - 1;
	    re = new double[degree +1];
	    im = new double[degree +1];
	}
	for (int i=0; i<= degree; i++) {
	    re[i] = a[i].re;
	    im[i] = a[i].im;
	}
    }


    /** Reads the coefficients of this polynomial from the array of
	complex numbers <pre>a<\pre> */
    public void assign(ComplexConstant[] a) {
	if ( degree != a.length - 1 || re == null ) {
	    degree = a.length - 1;
	    re = new double[degree +1];
	    im = new double[degree +1];
	}
	for (int i=0; i<= degree; i++) {
	    re[i] = a[i].getRe();
	    im[i] = a[i].getIm();
	}
    }

     /** Reads the coefficients of this polynomial from the array of
	complex numbers <pre>a<\pre> */
    public void assign(Field.Complex[] a) {
	if ( degree != a.length - 1 || re == null ) {
	    degree = a.length - 1;
	    re = new double[degree +1];
	    im = new double[degree +1];
	}
	for (int i=0; i<= degree; i++) {
	    re[i] = a[i].getRe();
	    im[i] = a[i].getIm();
	}
    }

   /** Assigns this with polynomial <pre>p<\pre>. */
    public void assign( ComplexPolynomial p ) {
	assign( p.re, p.im );
    }

    /** Reads the real and imaginary parts of the coefficients from the double
	arrays <pre>aRe<\pre> and <pre>aIm<\pre>. The degree is set to one plus
	the length of aRe. aRe and aIm are assumed to have the same length. */
    public void assign(double[] aRe, double[] aIm) {
	if (degree != aRe.length - 1 || re == null ) {
	    degree= aRe.length - 1;
	    re = (double[])aRe.clone();
	    im = (double[])aIm.clone();
	} else {
	    System.arraycopy(aRe,0, re, 0, degree + 1);
	    System.arraycopy(aIm,0, im, 0, degree + 1);
	}
    }

    /** Writes the coefficients (as far as they fit in) to the array
	<pre>aIm<\pre> of complex numbers. Creates new instances of Field.Complex
	only as needed. */
    public void getCoefficients(Complex[] a) {
	int minLength = Math.min(degree + 1, a.length);
	for (int i=0; i< minLength; i++) {
	    if ( a[i] == null ) {
		a[i] = new Complex(re[i], im[i]);
	    } else {
		a[i].assign(re[i], im[i]);
	    }
	}
    }

    /** Returns the coefficients as a freshly created array of complex
	numbers. */
    public Complex[] getCoefficients() {
	Complex[] c = new Complex[degree +1];
        getCoefficients(c);
	return c;
    }

    /** Returns the <pre>i<\pre>'th coefficient as a freshly created Field.Complex
	number. Returns zero if <pre>i<\pre> is negative or larger than the
	degree.*/
    public Complex getCoefficient(int i) {
	Complex z = new Complex();
	getCoefficient( i, z );
	return z;
    }

    /** Writes the value of the <pre>i<\pre>'th coefficient (zero if
	<pre>i<\pre> is negative or larger than the degree) to the complex
	number <pre>z<\pre>. */
    public void getCoefficient(int i, Complex z) {
	if ( i >= 0 && i <= degree ) {
	    z.assign( re[i], im[i] );
	} else {
	    z.assign(0,0);
	}
    }

    /** Reads the value of the <pre>i<\pre>'th coefficient from the complex
	number <pre>z<\pre>. <pre>i<\pre> is assumed to be non-negative and
	not bigger than the degree.*/
    public void setCoefficient(int i, Field.Complex z) {
	re[i] = z.getRe();
	im[i] = z.getIm();
    }


    /** Reads the value of the <pre>i<\pre>'th coefficient from the complex
	number <pre>x * iy<\pre>. <pre>i<\pre> is assumed to be non-negative and
	not bigger than the degree.*/
    public void setCoefficient(int i, double x, double y ) {
	re[i] = x;
	im[i] = y;
    }


    /** sets all coefficients to these in <pre>p<\pre> if they occur and to zero if not. */
    public void assignCoefficients( ComplexPolynomial p ) {

	final int minDegree = Math.min( degree, p.degree );

	System.arraycopy( p.re, 0, re, 0, minDegree + 1 );
	System.arraycopy( p.im, 0, im, 0, minDegree + 1 );

	for( int i = minDegree + 1; i <= degree; i++ ) {
	    re[ i ] = 0;
	    im[ i ] = 0;
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
	de.jtem.numericalMethods.algebra.polynomial.ComplexPolynomial.eval
	    ( re, im, degree, z.re, z.im, p );
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
	de.jtem.numericalMethods.algebra.polynomial.ComplexPolynomial.evalDerivative
	    ( re, im, degree, n, z.re, z.im, p );
	w.assign(p[0],p[1]);
    }

    /** Returns the sum of this polynomial and the polynomial <pre>q<\pre> as
	 a freshly created ComplexPolynomial. */
    public ComplexPolynomial plus(ComplexPolynomial q) {
	ComplexPolynomial result = new ComplexPolynomial( Math.max( degree, q.degree ) );
	result.assignPlus(this,q);
	return result;
    }

    /** Makes this polynomial equal to the sum of itself
	and the polynomial <pre>q<\pre>. */
    public void assignPlus(ComplexPolynomial q) {
	assignPlus(this, q);
    }

    /** Makes this polynomial equal to the sum of the polynomials <pre>q<\pre>
	and <pre>q<\pre>. */
    public void assignPlus( ComplexPolynomial p,
			    ComplexPolynomial q ) {
	int newDegree = Math.max(p.degree, q.degree);
	double[] newRe = new double[newDegree + 1];
	double[] newIm = new double[newDegree + 1];
	int i;
	for ( i=0; i <= p.degree; i++ ) {
	    newRe[i] = p.re[i];
	    newIm[i] = p.im[i];
	}
	for ( i=0; i<=q.degree; i++ ) {
	    newRe[i] += q.re[i];
	    newIm[i] += q.im[i];
	}

	re = newRe;
	im = newIm;

	degree = newDegree;
    }

    /** Returns the negative of this polynomial as a freshly
	created ComplexPolynomial. */
    public ComplexPolynomial neg() {
	double[] newRe = new double[degree +1];
	double[] newIm = new double[degree +1];
	for ( int i=0; i<= degree; i++) {
	    newRe[i] = -re[i];
	    newIm[i] = -im[i];
	}
	return new ComplexPolynomial( newRe, newIm );
    }

    /** Returns the difference of this polynomial and the polynomial
	<pre>q<\pre> as a freshly created ComplexPolynomial. */
    public ComplexPolynomial minus(ComplexPolynomial q) {
	ComplexPolynomial result = new ComplexPolynomial( Math.max( degree, q.degree ) );
	result.assignMinus(this,q);
	return result;
    }

    public void assignMinus(ComplexPolynomial q) {
	assignMinus(this, q);
    }

    /** Makes this polynomial equal to the difference of the polynomials
	<pre>q<\pre> and <pre>q<\pre>. */
    public void assignMinus( ComplexPolynomial p,
			     ComplexPolynomial q ) {
	int newDegree = Math.max(p.degree, q.degree);
	double[] newRe = new double[newDegree + 1];
	double[] newIm = new double[newDegree + 1];
	int i;
	for ( i=0; i <= p.degree; i++ ) {
	    newRe[i] = p.re[i];
	    newIm[i] = p.im[i];
	}
	for ( i=0; i<=q.degree; i++ ) {
	    newRe[i] -= q.re[i];
	    newIm[i] -= q.im[i];
	}
	re = newRe;
	im = newIm;
	degree = newDegree;
    }

    /** Returns the product of this polynomial and the polynomial <pre>q<\pre>
	as a freshly created ComplexPolynomial. */
    public ComplexPolynomial times(ComplexPolynomial q) {
	ComplexPolynomial result = new ComplexPolynomial( degree + q.degree );
	result.assignTimes(this,q);
	return result;
    }

    /** Makes this polynomial equal to the product of itself
	and the polynomial <pre>q<\pre>. */
    public void assignTimes(ComplexPolynomial q) {
	assignTimes(this, q);
    }

    /** Makes this polynomial equal to the product of the polynomials
	<pre>q<\pre> and <pre>q<\pre>. */
    public void assignTimes( ComplexPolynomial p,
			     ComplexPolynomial q ) {
	int newDegree = p.degree + q.degree;
	double[] newRe = new double[newDegree +1];
	double[] newIm = new double[newDegree +1];
	de.jtem.numericalMethods.algebra.polynomial.ComplexPolynomial.times
	    ( p.re, p.im, p.degree, q.re, q.im, q.degree, newRe, newIm );
	re = newRe;
	im = newIm;
	degree = newDegree;
    }

    /** Makes this polynomial equal to the product of the polynomial
	<pre>p<\pre> and the complex number <pre>z<\pre>. */
    public void assignTimes(ComplexPolynomial p, Complex z) {
	double pRe, pIm;
	setDegree(p.degree);
	for (int i=0; i<= degree; i++) {
	    pRe = p.re[i];
	    pIm = p.im[i];
	    re[i] = pRe * z.re - pIm * z.im;
	    im[i] = pRe * z.im + pIm * z.re;
	}
    }

    /** Makes this polynomial equal to the product of the polynomial
	<pre>p<\pre> and the complex number <pre>z<\pre>. */
    public void assignTimes(Complex z, ComplexPolynomial p) {
	assignTimes(p, z);
    }

    /** Makes this polynomial equal to the product of itself
	and the complex number <pre>z<\pre>. */
    public void assignTimes(Complex z) {
	assignTimes(this, z);
    }

    /** Returns the product of this polynomial and the complex
	number <pre>z<\pre> as a freshly created ComplexPolynomial. */
    public ComplexPolynomial times(Complex z) {
	ComplexPolynomial result = new ComplexPolynomial( degree );
	result.assignTimes(this, z);
	return result;
    }

    /** Makes this polynomial equal to the product of the polynomial
	<pre>p<\pre> and the real number <pre>x<\pre>. */
    public void assignTimes(ComplexPolynomial p, double x) {
	setDegree(p.degree);
	for (int i=0; i<= degree; i++) {
	    re[i] = p.re[i] * x;
	    im[i] = p.im[i] * x;
	}
    }

    /** Returns the product of this Laurent polynomial and the real
	number <pre>z<\pre> as a freshly created ComplexPolynomial. */
    public ComplexPolynomial times( double x) {
	ComplexPolynomial result = new ComplexPolynomial( degree );
	result.assignTimes(this, x);
	return result;
    }

    /** Makes this polynomial equal to the product of the polynomial
	<pre>p<\pre> and the real number <pre>x<\pre>. */
    public void assignTimes(double x, ComplexPolynomial p) {
	assignTimes(p, x);
    }

    /** Makes this polynomial equal to the product of itself
	and the real number <pre>z<\pre>. */
    public void assignTimes(double x) {
	assignTimes(this, x);
    }

    /** Assigns this with quotient of <pre>p<\pre> and <pre>q<\pre>. */
    public void assignDivide( ComplexPolynomial p, ComplexPolynomial q ) {
	double[] auxRe = new double[degree +1];
	double[] auxIm = new double[degree +1];
	double[] newRe = new double[p.degree -q.degree +1];
	double[] newIm = new double[p.degree -q.degree +1];
	de.jtem.numericalMethods.algebra.polynomial.ComplexPolynomial.divide(
	         p.re, p.im, p.degree, q.re, q.im, q.degree,
		 auxRe, auxIm, newRe, newIm);
	re = newRe;
	im = newIm;
	degree = p.degree -q.degree;
    }

    /** assigns this with <code>n</code>th derivative of <code>p</code>.
	<code>p</code> and <code>this</code> my coincide. */
    public void assignDerivative( final ComplexPolynomial p, final int n ) {

	if( n<0 )
	    throw new IllegalArgumentException( "order of derivative is negative" );

	if( degree < n )
	    assignZero();

	final double [] pRe = p.re; // this is for the case that this equals p
	final double [] pIm = p.im;

	final int degreeOfP = p.degree;

	if( degree != degreeOfP - n ) {

	    degree = degreeOfP - n;

	    re = new double[ degree + 1];
	    im = new double[ degree + 1];
	}

	de.jtem.numericalMethods.algebra.polynomial.ComplexPolynomial.derivative
	    ( pRe, pIm, degreeOfP, n, re, im );
    }

    /** assigns this with its <code>n</code>th derivative. */
    public void assignDerivative( final int n ) {
	assignDerivative( this, n );
    }

    /** returns <code>n</code>th derivative of this. */
    public ComplexPolynomial derivative( final int n ) {
	ComplexPolynomial result = new ComplexPolynomial( Math.max( degree - n, 0) );
	result.assignDerivative( this, n );
	return result;
    }

    /** assigns this with first derivative of <code>p</code>.
	<code>p</code> and <code>this</code> my coincide. */
    public void assignDerivative( final ComplexPolynomial p ) {
	this.assignDerivative( p, 1 );
    }

    /** assigns this with its first derivative. */
    public void assignDerivative() {
	assignDerivative( this, 1 );
    }

    /** returns first derivative of this. */
    public ComplexPolynomial derivative() {
	return derivative( 1 );
    }

    /** Returns the roots of this polynomial */
    public Complex[] getRoots() {
      Complex[] roots = new Complex[degree];
      getRoots( roots );
      return roots;
    }

  /** Returns the roots of this polynomial.
      @param roots array of length must coinside with degree of this. */
    public void getRoots( Complex [] roots ) {
	double[] rootsRe = new double[degree];
	double[] rootsIm = new double[degree];
	if( roots == null || roots.length != degree )
	  throw new IllegalArgumentException( "array has wrong size" );
	de.jtem.numericalMethods.algebra.polynomial.ComplexPolynomial.laguerAllRoots
	  (re, im, degree, EPS, rootsRe, rootsIm );
	for (int i=0; i<degree; i++) {
	  if( roots[i] == null ) {
	    roots[i] = new Complex( rootsRe[i], rootsIm[i] );
	  } else {
	    roots[i].assign( rootsRe[i], rootsIm[i] );
	  }
	}
    }

    /** Makes this a polynomial of degree given by the number of elements of
	the array roots, having the prescribed roots. If the current degree is
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
	ComplexPolynomial q = new ComplexPolynomial(1);
	for (int i=0; i<roots.length; i++) {
	    q.setCoefficient(0, roots[i].neg());
	    assignTimes(q);
	}
    }

    public void assignRandom() {
	for( int i=0; i<=degree; i++ ) {
	    re[i] = Math.random();
	    im[i] = Math.random();
	}
    }


    final private boolean isZero( final int index, final double epsSqr ) {
	return re[index]*re[index] + im[index]*im[index] < epsSqr;
    }

    final private boolean equals( final int index, ComplexPolynomial p, double epsSqr ) {
	final double dx = re[index]-p.re[index];
	final double dy = im[index]-p.im[index];

	return dx*dx + dy*dy < epsSqr;
    }


    /** test if coefficients of this and <pre>p<\pre> equals up to a precision of {#EPS}. */
    public boolean equals( ComplexPolynomial p ) {
	return equals( p, EPS );
    }

    /** test if coefficients of this and <pre>p<\pre> equals up to a precision of <pre>eps<\pre>. */
    public boolean equals( ComplexPolynomial p, double eps ) {
	final double epsSqr = eps*eps;
	if( p == this )
	    return true;

	if( degree > p.degree ) {
	    for( int i=p.degree+1; i<=degree; i++ )
		if( !isZero( i, epsSqr ) )
		    return false;
	} else {
	    for( int i=degree+1; i<=p.degree; i++ )
		if( !p.isZero( i, epsSqr ) )
		    return false;
	}

	final int max = Math.min( p.degree, degree );

	for( int i=0; i<=max; i++ )
	    if( !equals( i, p, epsSqr ) )
		return false;

	return true;
    }

    /** test if coefficients of this and <pre>o<\pre> equals up to a precision of {#EPS}. */
    public boolean equals( Object o ) {
	if( o == this )
	    return true;

	try {
	    return equals( (ComplexPolynomial)o );
	} catch ( ClassCastException e ) {
	    return false;
	}
    }

    public String toString() {
	return 	de.jtem.numericalMethods.algebra.polynomial.ComplexPolynomial.toString
	    ( re, im, degree, EPS );
    }


}

