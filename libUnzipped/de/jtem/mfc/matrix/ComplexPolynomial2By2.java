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

package de.jtem.mfc.matrix;

import java.io.Serializable;

import de.jtem.mfc.field.Complex;
import de.jtem.mfc.field.Field;
import de.jtem.mfc.polynomial.ComplexPolynomial;

/**
 * This class represents a 2-by-2 matrix of complex polynomials as entries.
 * The functionality closely resembles that of {@link de.jtem.mfc.matrix.Complex2By2}
 * matrices.
 */
public class ComplexPolynomial2By2 implements Serializable {
    private static final long serialVersionUID = 1L;

    /** 1e-14; constant for zero tests and precision demands. */
    static final double EPS = 1.e-14;

    public ComplexPolynomial a, b, c, d;

    static ComplexPolynomial p11 = new ComplexPolynomial();
    static ComplexPolynomial p12 = new ComplexPolynomial();
    static ComplexPolynomial p21 = new ComplexPolynomial();
    static ComplexPolynomial p22 = new ComplexPolynomial();
    static ComplexPolynomial q11 = new ComplexPolynomial();
    static ComplexPolynomial q12 = new ComplexPolynomial();
    static ComplexPolynomial q21 = new ComplexPolynomial();
    static ComplexPolynomial q22 = new ComplexPolynomial();


    public ComplexPolynomial2By2() {
	a = new ComplexPolynomial();
	b = new ComplexPolynomial();
	c = new ComplexPolynomial();
	d = new ComplexPolynomial();
    }

    public ComplexPolynomial2By2( ComplexPolynomial2By2 p ) {
	this( p.a, p.b, p.c, p.d ) ;
    }

    public ComplexPolynomial2By2( ComplexPolynomial a,
				  ComplexPolynomial b,
				  ComplexPolynomial c,
				  ComplexPolynomial d ) {

	this.a = new ComplexPolynomial( a );
	this.b = new ComplexPolynomial( b );
	this.c = new ComplexPolynomial( c );
	this.d = new ComplexPolynomial( d );
    }

    public ComplexPolynomial getA() {
	return a;
    }

    public void setA(ComplexPolynomial p) {
	a = p;
    }

    public ComplexPolynomial getB() {
	return b;
    }

    public void setB(ComplexPolynomial p) {
	b = p;
    }

    public ComplexPolynomial getC() {
	return c;
    }

    public void setC(ComplexPolynomial p) {
	c = p;
    }
    public ComplexPolynomial getD() {
	return d;
    }

    public void setD(ComplexPolynomial p) {
	d = p;
    }

    private static final Complex evalResult = new Complex();

    public void eval(Complex z, AbstractComplex2By2 m) {
	a.eval( z, evalResult ); m.setA( evalResult );
	b.eval( z, evalResult ); m.setB( evalResult );
	c.eval( z, evalResult ); m.setC( evalResult );
	d.eval( z, evalResult ); m.setD( evalResult );
    }

    public Complex2By2 eval(Complex z) {
	Complex2By2 m = new Complex2By2();
	eval( z, m );
	return m;
    }

    private static final Complex evalDerivativeResult = new Complex();

    public void evalDerivative( Complex z, int n, AbstractComplex2By2 m ) {
	a.evalDerivative( z, n, evalDerivativeResult ); m.setA( evalDerivativeResult );
	b.evalDerivative( z, n, evalDerivativeResult ); m.setB( evalDerivativeResult );
	c.evalDerivative( z, n, evalDerivativeResult ); m.setC( evalDerivativeResult );
	d.evalDerivative( z, n, evalDerivativeResult ); m.setD( evalDerivativeResult );
    }

    public Complex2By2 evalDerivative( Complex z, int n ) {
	Complex2By2 m = new Complex2By2();
	evalDerivative( z, n, m );
	return m;
    }

    public void assign( ComplexPolynomial2By2 p ) {
	a.assign( p.a );
	b.assign( p.b );
	c.assign( p.c );
	d.assign( p.d );
    }

    public static ComplexPolynomial2By2 constant(
        double aRe, double aIm, double bRe, double bIm,
        double cRe, double cIm, double dRe, double dIm) {
      ComplexPolynomial2By2 result = new ComplexPolynomial2By2();
      result.getA().setCoefficient(0, aRe, aIm);
      result.getB().setCoefficient(0, bRe, bIm);
      result.getC().setCoefficient(0, cRe, cIm);
      result.getD().setCoefficient(0, dRe, dIm);
      return result;
    }

    public static ComplexPolynomial2By2 constant(Complex a,
        Complex b, Complex c, Complex d) {
      return constant(a.re, a.im, b.re, b.im, c.re, c.im, d.re, d.im);
    }

    public static ComplexPolynomial2By2 constant(Field.Complex a,
        Field.Complex b, Field.Complex c, Field.Complex d) {
      return constant(a.getRe(), a.getIm(), b.getRe(), b.getIm(), 
		      c.getRe(), c.getIm(), d.getRe(), d.getIm());
    }

    public static ComplexPolynomial2By2 identity() {
      return constant(1, 0, 0, 0, 0, 0, 1, 0);
    }

    public void assignIdentity() {
      a.assignMonomial(0, 1, 0);
      b.assignMonomial(0, 0, 0);
      c.assignMonomial(0, 0, 0);
      d.assignMonomial(0, 1, 0);
    }

    public void assignConstant( Field.Complex a, Field.Complex b, Field.Complex c, Field.Complex d ) {
	this.a.assignMonomial(0, a );
	this.b.assignMonomial(0, b );
	this.c.assignMonomial(0, c );
	this.d.assignMonomial(0, d );
    }

    public void assignTimes( ComplexPolynomial2By2 l,
			     ComplexPolynomial2By2 m) {

	p11.assignTimes(l.a, m.a);
	q11.assignTimes(l.b, m.c);
	p12.assignTimes(l.a, m.b);
	q12.assignTimes(l.b, m.d);
	p21.assignTimes(l.c, m.a);
	q21.assignTimes(l.d, m.c);
	p22.assignTimes(l.c, m.b);
	q22.assignTimes(l.d, m.d);

	a.assignPlus(p11, q11);
	b.assignPlus(p12, q12);
	c.assignPlus(p21, q21);
	d.assignPlus(p22, q22);
    }

    public void assignTimes(ComplexPolynomial2By2 m) {
	assignTimes(this, m);
    }

    public ComplexPolynomial2By2 times(ComplexPolynomial2By2 m) {
	ComplexPolynomial2By2 result = new ComplexPolynomial2By2();
	result.assignTimes(this, m);
	return result;
    }

    public void assignPlus( ComplexPolynomial2By2 l,
			    ComplexPolynomial2By2 m) {
	a.assignPlus(l.a, m.a);
	b.assignPlus(l.b, m.b);
	c.assignPlus(l.c, m.c);
	d.assignPlus(l.d, m.d);
    }

    public void assignPlus( ComplexPolynomial2By2 m ) {
	assignPlus(this, m);
    }

    public ComplexPolynomial2By2 plus(ComplexPolynomial2By2 m) {
	ComplexPolynomial2By2 result = new ComplexPolynomial2By2();
	result.assignPlus(this, m);
	return result;
    }

    public ComplexPolynomial2By2 neg() {
	ComplexPolynomial2By2 result = new ComplexPolynomial2By2();
	result.setA(a.neg());
	result.setB(b.neg());
	result.setC(c.neg());
	result.setD(d.neg());
	return result;
    }

    public ComplexPolynomial2By2 minus(ComplexPolynomial2By2 q) {
	ComplexPolynomial2By2 result = new ComplexPolynomial2By2();
	result.assignMinus(this, q);
	return result;
    }

    public void assignMinus( ComplexPolynomial2By2 q ) {
	assignMinus(this, q );
    }

    public void assignMinus( ComplexPolynomial2By2 p,
			     ComplexPolynomial2By2 q ) {
	a.assignMinus(p.a, q.a);
	b.assignMinus(p.b, q.b);
	c.assignMinus(p.c, q.c);
	d.assignMinus(p.d, q.d);
    }

    public void getDeterminant(ComplexPolynomial p) {
	p.assignMinus(a.times(d) , b.times(c));
    }

    public ComplexPolynomial getDeterminant() {
	ComplexPolynomial result = new ComplexPolynomial();
	getDeterminant(result);
	return result;
    }

    public void getTrace(ComplexPolynomial p) {
	p.assignPlus(a,d);
    }

    public ComplexPolynomial getTrace() {
	ComplexPolynomial result = new ComplexPolynomial();
	getTrace(result);
	return result;
    }

    public int getDegree() {
	int deg = a.getDegree();
	int newDeg = b.getDegree();
	if (newDeg > deg) deg = newDeg;
	newDeg = c.getDegree();
	if (newDeg > deg) deg = newDeg;
	newDeg = d.getDegree();
	if (newDeg > deg) deg = newDeg;
	return deg;
    }

    public void setDegree(int deg) {
	a.setDegree(deg);
	b.setDegree(deg);
	c.setDegree(deg);
	d.setDegree(deg);
    }

    public void truncate(int deg) {
	a.truncate(deg);
	b.truncate(deg);
	c.truncate(deg);
	d.truncate(deg);
    }

    public void truncate( final double eps ) {
	a.truncate(eps);
	b.truncate(eps);
	c.truncate(eps);
	d.truncate(eps);
    }

    public void truncate() {
	truncate( EPS );
    }

    public Complex[] getPointsOfPeriod(int n) {
	ComplexPolynomial t = getTrace();
	t.assignTimes(t);

	ComplexPolynomial det = getDeterminant();
	int deg = Math.max(t.getDegree(), det.getDegree());
	Complex[] pts = new Complex[deg * (n/2 + 1)];
	int j,k;
	for ( k=0; k<= n/2; k++) {
	    double c = Math.cos( Math.PI * k / n );
	    double c2 = c*c;
	    Complex[] roots = t.minus(det.times( 4*c2 )).getRoots();
	    for ( j=0; j < deg; j++) {
		pts[k*deg + j] = roots[j];
	    }
	}
	return pts;
    }


    /** assigns this with <code>n</code>th derivative of <code>p</code>.
	<code>p</code> and <code>this</code> my coincide. */
    public void assignDerivative( final ComplexPolynomial2By2 p, final int n ) {

	a.assignDerivative( p.a, n );
	b.assignDerivative( p.b, n );
	c.assignDerivative( p.c, n );
	d.assignDerivative( p.d, n );
    }

    /** assigns this with its <code>n</code>th derivative. */
    public void assignDerivative( final int n ) {
	assignDerivative( this, n );
    }

    /** returns <code>n</code>th derivative of this. */
    public ComplexPolynomial2By2 derivative( final int n ) {
	ComplexPolynomial2By2 result = new ComplexPolynomial2By2( a.derivative( n ), b.derivative( n ),
								  c.derivative( n ), d.derivative( n ) );
	return result;
    }

    /** assigns this with first derivative of <code>p</code>.
	<code>p</code> and <code>this</code> my coincide. */
    public void assignDerivative( final ComplexPolynomial2By2 p ) {
	this.assignDerivative( p, 1 );
    }

    /** assigns this with its first derivative. */
    public void assignDerivative() {
	assignDerivative( this, 1 );
    }

    /** returns first derivative of this. */
    public ComplexPolynomial2By2 derivative() {
	return derivative( 1 );
    }

    public void assignRandom() {
	a.assignRandom();
	b.assignRandom();
	c.assignRandom();
	d.assignRandom();
    }

    /** test if coefficients of this and <pre>p<\pre> equals up to a precision of {#EPS}. */
    public boolean equals( ComplexPolynomial2By2 p ) {
	return equals( p, EPS );
    }

    /** test if coefficients of this and <pre>p<\pre> equals up to a precision of <pre>eps<\pre>. */
    public boolean equals( ComplexPolynomial2By2 p, double eps ) {
	return a.equals( p.a ) && b.equals( p.b ) && c.equals( p.c ) && d.equals( p.d );
    }

    /** test if coefficients of this and <pre>o<\pre> equals up to a precision of {#EPS}. */
    public boolean equals( Object o ) {
	if( o == this )
	    return true;

	try {
	    return equals( (ComplexPolynomial2By2)o );
	} catch ( ClassCastException e ) {
	    return false;
	}
    }

    public String toString() {
	StringBuffer sb=new StringBuffer(300);

	sb.append( "a = " + a + "\n" );
	sb.append( "b = " + b + "\n" );
	sb.append( "c = " + c + "\n" );
	sb.append( "d = " + d + "\n" );

	return sb.toString();
    }

}

