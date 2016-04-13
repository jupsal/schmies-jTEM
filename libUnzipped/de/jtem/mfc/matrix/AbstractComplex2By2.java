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
import de.jtem.mfc.vector.Complex2;

/** This class is the base for all classes that are or can be represented as
 ** a complex 2-by-2 matrix.
 **
 ** <p> Most of its functionallity is implemented
 ** in protected methods. This enables derived classes to specialize and omit
 ** certain functions. For example, it does not make any sense to add {@link
 ** de.jtem.mfc.group.Moebius} transformations.
 **
 ** <p> To minimize the object count its
 ** 4 complex entries <code>a,b,c,d</code> are represented by 8 double values
 ** <code>{@link #aRe}, {@link #aIm}, {@link #bRe}, {@link #bIm}, {@link
 ** #cRe}, {@link #cIm}, {@link #dRe}, {@link #dIm}</code>, e.g.:
 **
 ** <pre>
 ** [ a   b ]   [ aRe + i aIm    bRe + i bIm ]
 ** [       ] = [                            ]
 ** [ c   d ]   [ dRe + i dIm    dRe + i dIm ]
 ** </pre>
 **
 ** <p> All methods avoid the creation of temporarily used objects
 ** in their internal computations unless the opposite is stated in their
 ** descriptions.
 **
 ** <p>This class has <code>final</code> member objects, which are
 ** used as dummy variables. Thus, it is definitely not thread-safe.
 **
 ** @author schmies */

public abstract class AbstractComplex2By2 implements Serializable, Cloneable {

    private static final long serialVersionUID = 1L;


    /** @deprecated superceded by EPS */
    public static final double EPSILON = 1.e-14;

    /** threashold for zero tests; default is 1e-14 */
    public static final double EPS = 1.e-14;

    /** the square of {@link #EPS} */
    public static final double EPSSQR = EPS * EPS;

    /**
     * Dummy variable used in calculations.
     */
    static final Complex     dummyComplex     = new Complex();

    /**
     * Dummy variable used in calculations.
     */
    static final Complex2By2 dummyComplex2By2 = new Complex2By2();

    /** entry of the matrix */
    protected double aRe, aIm, bRe, bIm, cRe, cIm, dRe, dIm;

    /** creates an identity matrix */
    protected AbstractComplex2By2() {
	assignIdentity();
    }

    /** creates a matrix with the prescribed entries */
    protected AbstractComplex2By2( final double aRe, final double aIm,
				   final double bRe, final double bIm,
				   final double cRe, final double cIm,
				   final double dRe, final double dIm) {
	assign( aRe, aIm, bRe, bIm, cRe, cIm, dRe, dIm );
    }

    /** creates a matrix with the prescribed entries */
    protected AbstractComplex2By2( final Complex a, final Complex b,
                                   final Complex c, final Complex d) {
        assign( a.re, a.im, b.re, b.im, c.re, c.im, d.re, d.im);
    }

    /** creates a matrix with the prescribed entries */
    protected AbstractComplex2By2(final Field.Complex a,
                                  final Field.Complex b,
                                  final Field.Complex c,
                                  final Field.Complex d) {
      assign(a.getRe(), a.getIm(), b.getRe(), b.getIm(),
             c.getRe(), c.getIm(), d.getRe(), d.getIm());
    }

    /** creates a matrix equaling prescribed matrix */
    public AbstractComplex2By2( final AbstractComplex2By2 m ) {
	assign( m );
    }

    /** returns entry <code>a</code> of this */
    public Complex getA() {
	return new Complex( aRe, aIm );
    }

    /** assigns <code>a</code> with entry <code>a</code> of this
     ** @param a complex value which is assigned with entry <code>a</code>. */
    public void getA( final Complex a ) {
	a.assign( aRe, aIm );
    }

    /** sets entry <code>a</code> of this with <code>a</code> */
    protected void setA( final Complex a ) {
	aRe = a.re;
	aIm = a.im;
    }

    /** returns entry <code>b</code> of this */
    public Complex getB() {
	return new Complex( bRe, bIm );
    }

    /** assigns <code>b</code> with entry <code>b</code> of this
     ** @param b complex value which is assigned with entry <code>b</code>.*/
    public void getB( final Complex b ) {
	b.assign( bRe, bIm );
    }

    /** sets entry <code>b</code> of this with <code>b</code> */
    protected void setB( final Complex b ) {
	bRe = b.re;
	bIm = b.im;
    }

    /** returns entry <code>c</code> of this */
    public Complex getC() {
	return new Complex( cRe, cIm );
    }

    /** assigns <code>c</code> with entry <code>c</code> of this
     ** @param c complex value which is assigned with entry <code>c</code>. */
    public void getC( final Complex c ) {
	c.assign( cRe, cIm );
    }

    /** sets entry <code>c</code> of this with <code>c</code> */
    protected void setC( final Complex c ) {
	cRe = c.re;
	cIm = c.im;
    }

    /** returns entry <code>d</code> of this */
    public Complex getD() {
	return new Complex( dRe, dIm );
    }

    /** assigns <code>d</code> with entry <code>d</code> of this
     ** @param d complex value which is assigned with entry <code>d</code>. */
    public void getD( final Complex d ){
	d.assign( dRe, dIm );
    }

    /** sets entry <code>d</code> of this with <code>d</code> */
    protected void setD( final Complex d ) {
	dRe = d.re;
	dIm = d.im;
    }

    /** assigns this with the zero matrix */
    protected void assignZero() {
	assign( 0, 0, 0, 0, 0, 0, 0, 0 );
    }

    /** assigns this with the identity matrix */
    protected void assignIdentity() {
	assign( 1, 0, 0, 0, 0, 0, 1, 0 );
    }

    /** assigns this with the prescribed entries */
    protected void assign( double aRe, double aIm, double bRe, double bIm,
			   double cRe, double cIm, double dRe, double dIm) {
	this.aRe = aRe; this.aIm = aIm; this.bRe = bRe; this.bIm = bIm;
	this.cRe = cRe; this.cIm = cIm; this.dRe = dRe; this.dIm = dIm;
    }

    /** assigns this with the prescribed entries */
    protected void assign( final Complex a, final Complex b,
                              final Complex c, final Complex d ) {
	assign( a.re, a.im, b.re, b.im, c.re, c.im, d.re, d.im );
    }

    /** assigns this with the prescribed matrix */
    protected void assign( final AbstractComplex2By2 s ) {
	assign( s.aRe, s.aIm, s.bRe, s.bIm, s.cRe, s.cIm, s.dRe, s.dIm );
    }

    /** assigns this with the product of itself and a */
    protected void assignTimes( final AbstractComplex2By2 a ) {
	assignTimes(this, a );
    }

    /**
     * multiplies <code>a</code> and <code>b</code> and stores
     * the result in <code>c</code>.
     */
    public static void times( final AbstractComplex2By2 a,
			      final AbstractComplex2By2 b,
			      final AbstractComplex2By2 c ) {

	if( a == c || b == c) {
	    double d1r = a.aRe*b.aRe-a.aIm*b.aIm + a.bRe*b.cRe-a.bIm*b.cIm;
	    double d1i = a.aRe*b.aIm+a.aIm*b.aRe + a.bRe*b.cIm+a.bIm*b.cRe;
	    double d2r = a.aRe*b.bRe-a.aIm*b.bIm + a.bRe*b.dRe-a.bIm*b.dIm;
	    double d2i = a.aRe*b.bIm+a.aIm*b.bRe + a.bRe*b.dIm+a.bIm*b.dRe;
	    double d3r = a.cRe*b.aRe-a.cIm*b.aIm + a.dRe*b.cRe-a.dIm*b.cIm;
	    double d3i = a.cRe*b.aIm+a.cIm*b.aRe + a.dRe*b.cIm+a.dIm*b.cRe;
	    double d4r = a.cRe*b.bRe-a.cIm*b.bIm + a.dRe*b.dRe-a.dIm*b.dIm;
	    double d4i = a.cRe*b.bIm+a.cIm*b.bRe + a.dRe*b.dIm+a.dIm*b.dRe;
	    c.aRe=d1r; c.aIm=d1i; c.bRe=d2r; c.bIm=d2i; c.cRe=d3r; c.cIm=d3i; c.dRe=d4r; c.dIm=d4i;
	} else {
	    c.aRe = a.aRe*b.aRe-a.aIm*b.aIm + a.bRe*b.cRe-a.bIm*b.cIm;
	    c.aIm = a.aRe*b.aIm+a.aIm*b.aRe + a.bRe*b.cIm+a.bIm*b.cRe;
	    c.bRe = a.aRe*b.bRe-a.aIm*b.bIm + a.bRe*b.dRe-a.bIm*b.dIm;
	    c.bIm = a.aRe*b.bIm+a.aIm*b.bRe + a.bRe*b.dIm+a.bIm*b.dRe;
	    c.cRe = a.cRe*b.aRe-a.cIm*b.aIm + a.dRe*b.cRe-a.dIm*b.cIm;
	    c.cIm = a.cRe*b.aIm+a.cIm*b.aRe + a.dRe*b.cIm+a.dIm*b.cRe;
	    c.dRe = a.cRe*b.bRe-a.cIm*b.bIm + a.dRe*b.dRe-a.dIm*b.dIm;
	    c.dIm = a.cRe*b.bIm+a.cIm*b.bRe + a.dRe*b.dIm+a.dIm*b.dRe;
	}
    }

    /** assigns this with the product of <code>a</code> and <code>b</code> */
    protected void assignTimes( final AbstractComplex2By2 a, final AbstractComplex2By2 b ) {
	times( a, b, this );
    }

    /** assigns this with the product of itself and <code>r</code> */
    protected void assignTimes( final double r ) {
	assignTimes( this, r );
    }

    /** assigns this with the product of <code>m</code> and <code>r</code> */
    protected void assignTimes( final AbstractComplex2By2 m, final double r ) {
	if( m == this ) {

	    aRe *= r; aIm *= r;
	    bRe *= r; bIm *= r;
	    cRe *= r; cIm *= r;
	    dRe *= r; dIm *= r;
	} else {

	    aRe = m.aRe * r; aIm = m.aIm * r;
	    bRe = m.bRe * r; bIm = m.bIm * r;
	    cRe = m.cRe * r; cIm = m.cIm * r;
	    dRe = m.dRe * r; dIm = m.dIm * r;
	}
    }

    /** assigns this with the product of itself and <code>x+iy</code> */
    protected void assignTimes( final double x, final double y ) {
	assignTimes( this, x, y );
    }

    /** assigns this with the product of <code>m</code> and <code>x+iy</code> */
    protected void assignTimes( final AbstractComplex2By2 m, final double x, final double y ) {

	if( m == this ) {
	    double tmp;

	    tmp = aRe * x - aIm * y;
	    aIm = aIm * x + aRe * y;
	    aRe = tmp;

	    tmp = bRe * x - bIm * y;
	    bIm = bIm * x + bRe * y;
	    bRe = tmp;

	    tmp = cRe * x - cIm * y;
	    cIm = cIm * x + cRe * y;
	    cRe = tmp;

	    tmp = dRe * x - dIm * y;
	    dIm = dIm * x + dRe * y;
	    dRe = tmp;

	} else {
	    aRe = m.aRe * x - m.aIm * y;
	    aIm = m.aIm * x + m.aRe * y;

	    bRe = m.bRe * x - m.bIm * y;
	    bIm = m.bIm * x + m.bRe * y;

	    cRe = m.cRe * x - m.cIm * y;
	    cIm = m.cIm * x + m.cRe * y;

	    dRe = m.dRe * x - m.dIm * y;
	    dIm = m.dIm * x + m.dRe * y;
	}
    }

    /** returns product of this and <code>z</code> */
    protected void assignTimes( final Complex z ) {
	assignTimes( this, z.re, z.im );
    }

    /** assigns this with the product of <code>m</code> and <code>z</code> */
    protected void assignTimes( final AbstractComplex2By2 m, final Complex z ) {
	assignTimes( m, z.re, z.im );
    }

    /** assigns this with the product of itself and inverse of <code>m</code>*/
    protected void assignDivide( final AbstractComplex2By2 m ) {
	assignDivide( this, m );
    }

    /**
     * multiplies <code>a</code> and inverse of <code>b</code> and stores
     * the result in <code>c</code>.
     */
    public static void divide( final AbstractComplex2By2 a,
			       final AbstractComplex2By2 b,
			       final AbstractComplex2By2 c ) {

	double detRe = b.aRe * b.dRe - b.aIm * b.dIm  -  ( b.bRe * b.cRe - b.bIm * b.cIm );
	double detIm = b.aRe * b.dIm + b.aIm * b.dRe  -  ( b.bRe * b.cIm + b.bIm * b.cRe );

	double nn = detRe * detRe + detIm * detIm;

	if( nn == 0 )
	    throw new RuntimeException( "Matrix is not invertible. " );

	if( a==c || b==c ) {
	    double d1r = a.aRe*b.dRe-a.aIm*b.dIm - a.bRe*b.cRe+a.bIm*b.cIm;
	    double d1i = a.aRe*b.dIm+a.aIm*b.dRe - a.bRe*b.cIm-a.bIm*b.cRe;
	    double d2r = a.bRe*b.aRe-a.bIm*b.aIm - a.aRe*b.bRe+a.aIm*b.bIm;
	    double d2i = a.bRe*b.aIm+a.bIm*b.aRe - a.aRe*b.bIm-a.aIm*b.bRe;
	    double d3r = a.cRe*b.dRe-a.cIm*b.dIm - a.dRe*b.cRe+a.dIm*b.cIm;
	    double d3i = a.cRe*b.dIm+a.cIm*b.dRe - a.dRe*b.cIm-a.dIm*b.cRe;
	    double d4r = a.dRe*b.aRe-a.dIm*b.aIm - a.cRe*b.bRe+a.cIm*b.bIm;
	    double d4i = a.dRe*b.aIm+a.dIm*b.aRe - a.cRe*b.bIm-a.cIm*b.bRe;
	    c.aRe=d1r; c.aIm=d1i; c.bRe=d2r; c.bIm=d2i; c.cRe=d3r; c.cIm=d3i; c.dRe=d4r; c.dIm=d4i;
	} else {
	    c.aRe = a.aRe*b.dRe-a.aIm*b.dIm - a.bRe*b.cRe+a.bIm*b.cIm;
	    c.aIm = a.aRe*b.dIm+a.aIm*b.dRe - a.bRe*b.cIm-a.bIm*b.cRe;
	    c.bRe = a.bRe*b.aRe-a.bIm*b.aIm - a.aRe*b.bRe+a.aIm*b.bIm;
	    c.bIm = a.bRe*b.aIm+a.bIm*b.aRe - a.aRe*b.bIm-a.aIm*b.bRe;
	    c.cRe = a.cRe*b.dRe-a.cIm*b.dIm - a.dRe*b.cRe+a.dIm*b.cIm;
	    c.cIm = a.cRe*b.dIm+a.cIm*b.dRe - a.dRe*b.cIm-a.dIm*b.cRe;
	    c.dRe = a.dRe*b.aRe-a.dIm*b.aIm - a.cRe*b.bRe+a.cIm*b.bIm;
	    c.dIm = a.dRe*b.aIm+a.dIm*b.aRe - a.cRe*b.bIm-a.cIm*b.bRe;
	}

	c.assignDivide( nn );
    }


    /** assigns this with the product of <code>a</code> and inverse of <code>b</code> */
    protected void assignDivide( final AbstractComplex2By2 a, final AbstractComplex2By2 b ) {
	divide( a, b, this );
    }

    /** assigns this with the product of itself and <code>1/r</code> */
    protected void assignDivide( final double r ) {
	assignDivide( this, r );
    }

    /** assigns this with the product of <code>m</code> and <code>1/r</code> */
    protected void assignDivide( final AbstractComplex2By2 m, final double r ) {
	if( m == this ) {

	    aRe /= r; aIm /= r;
	    bRe /= r; bIm /= r;
	    cRe /= r; cIm /= r;
	    dRe /= r; dIm /= r;
	} else {

	    aRe = m.aRe / r; aIm = m.aIm / r;
	    bRe = m.bRe / r; bIm = m.bIm / r;
	    cRe = m.cRe / r; cIm = m.cIm / r;
	    dRe = m.dRe / r; dIm = m.dIm / r;
	}
    }

    /** assigns this with the product of itself and <code>1/(x+iy)</code> */
    protected void assignDivide( final double x, final double y ) {
	assignDivide( this, x, y );
    }

    /** assigns this with the product of <code>m</code> and <code>1/(x+iy)</code> */
    protected void assignDivide( final AbstractComplex2By2 m, final double x, final double y ) {

	double nn = x * x + y * y;

	if( nn != 0 ) {
	    assignTimes( m, x, -y );
	}

	assignDivide( nn );
    }

    /** assigns this with the product of itself and <code>1/z</code> */
    protected void assignDivide( final Complex z ) {
	assignDivide( this, z.re, z.im );
    }

    /** assigns this with the product of <code>m</code> and <code>1/z</code> */
    protected void assignDivide( final AbstractComplex2By2 m, final Complex z ) {

	double nn = z.re * z.re + z.im * z.im;

	if( nn != 0 ) {
	    assignTimes( m, z.re, -z.im );
	}

	assignDivide( nn );
    }



    /** assigns this with the sum of itself and a */
    protected void assignPlus( final AbstractComplex2By2 a ) {
	assignPlus(this, a );
    }

    /**
     * adds <code>a</code> and <code>b</code> and stores
     * the result in <code>c</code>.
     */
    public static void plus( final AbstractComplex2By2 a,
			     final AbstractComplex2By2 b,
			     final AbstractComplex2By2 c ) {

	c.aRe = a.aRe + b.aRe;
	c.aIm = a.aIm + b.aIm;

	c.bRe = a.bRe + b.bRe;
	c.bIm = a.bIm + b.bIm;

	c.cRe = a.cRe + b.cRe;
	c.cIm = a.cIm + b.cIm;

	c.dRe = a.dRe + b.dRe;
	c.dIm = a.dIm + b.dIm;

    }

    /** assigns this with the sum of <code>a</code> and <code>b</code> */
    protected void assignPlus( final AbstractComplex2By2 a, final AbstractComplex2By2 b ) {
	plus( a, b, this );
    }

    /** adds <code>x</code> to all entries of this */
    protected void assignPlus( final double x ) {
	assignPlus( this, x );
    }

    /** adds <code>x</code> to all entries of <code>m</code> and
	and assigns this with the result. */
    protected void assignPlus( final AbstractComplex2By2 m, final double x ) {

	aRe = m.aRe + x;
	bRe = m.bRe + x;
	cRe = m.cRe + x;
	dRe = m.dRe + x;
    }

    /** adds <code>x+iy</code> to all entries of this */
    protected void assignPlus( final double x, final double y ) {
	assignPlus( this, x, y );
    }

    /** adds <code>x+iy</code> to all entries of <code>m</code> and
	and assigns this with the result. */
    protected void assignPlus( final AbstractComplex2By2 m, final double x, final double y ) {

	aRe = m.aRe + x;
	aIm = m.aIm + y;

	bRe = m.bRe + x;
	bIm = m.bIm + y;

	cRe = m.cRe + x;
	cIm = m.cIm + y;

	dRe = m.dRe + x;
	dIm = m.dIm + y;
    }

    /** adds <code>z</code> to all entries of this */
    protected void assignPlus( final Complex z ) {
	assignPlus( this, z.re, z.im );
    }

    /** adds <code>z</code> to all entries of <code>m</code> and
	and assigns this with the result. */
    protected void assignPlus( final AbstractComplex2By2 m, final Complex z ) {
	assignPlus( m, z.re, z.im );
    }

    /** subtracts <code>a</code> from this */
    protected void assignMinus( final AbstractComplex2By2 a ) {
	assignMinus(this, a );
    }

    /**
     * subtracts <code>b</code> from <code>a</code> and stores
     * the result in <code>c</code>.
     */
    public static void minus( final AbstractComplex2By2 a,
			      final AbstractComplex2By2 b,
			      final AbstractComplex2By2 c ) {

	c.aRe = a.aRe - b.aRe;
	c.aIm = a.aIm - b.aIm;

	c.bRe = a.bRe - b.bRe;
	c.bIm = a.bIm - b.bIm;

	c.cRe = a.cRe - b.cRe;
	c.cIm = a.cIm - b.cIm;

	c.dRe = a.dRe - b.dRe;
	c.dIm = a.dIm - b.dIm;

    }

    /** subtracts <code>b</code> from <code>a</code> and stores the result it this */
    protected void assignMinus( final AbstractComplex2By2 a, final AbstractComplex2By2 b ) {
	minus( a, b, this );
    }

    /** subtracts <code>x</code> from all entries of this */
    protected void assignMinus( final double x ) {
	assignMinus( this, x );
    }

    /** subtracts <code>x</code> from all entries of <code>m</code> and
	and assigns this with the result. */
    protected void assignMinus( final AbstractComplex2By2 m, final double x ) {

	aRe = m.aRe - x;
	bRe = m.bRe - x;
	cRe = m.cRe - x;
	dRe = m.dRe - x;
    }

    /** subtracts <code>x+iy</code> to all entries of this */
    protected void assignMinus( final double x, final double y ) {
	assignMinus( this, x, y );
    }

    /** subtracts <code>x+iy</code> to all entries of <code>m</code> and
	and assigns this with the result. */
    protected void assignMinus( final AbstractComplex2By2 m, final double x, final double y ) {

	aRe = m.aRe - x;
	aIm = m.aIm - y;

	bRe = m.bRe - x;
	bIm = m.bIm - y;

	cRe = m.cRe - x;
	cIm = m.cIm - y;

	dRe = m.dRe - x;
	dIm = m.dIm - y;
    }

    /** subtracts <code>z</code> from all entries of this */
    protected void assignMinus( final Complex z ) {
	assignMinus( this, z.re, z.im );
    }

    /** subtracts <code>z</code> to all entries of <code>m</code> and
	and assigns this with the result. */
    protected void assignMinus( final AbstractComplex2By2 m, final Complex z ) {
	assignMinus( m, z.re, z.im );
    }


    /** assigns this with the negative of <code>a</code> */
    protected void assignNeg( final AbstractComplex2By2 a ) {

	aRe = -a.aRe;
	aIm = -a.aIm;

	bRe = -a.bRe;
	bIm = -a.bIm;

	cRe = -a.cRe;
	cIm = -a.cIm;

	dRe = -a.dRe;
	dIm = -a.dIm;
    }

    /** assigns this with its own negative */
    protected void assignNeg() {
	assignNeg( this );
    }

    /** assigns this with the conjugateative of <code>a</code> */
    protected void assignConjugate( final AbstractComplex2By2 a ) {

	if( a == this ) {

	    aIm = -a.aIm;
	    bIm = -a.bIm;
	    cIm = -a.cIm;
	    dIm = -a.dIm;

	} else {

	    aRe =  a.aRe;
	    aIm = -a.aIm;

	    bRe =  a.cRe;
	    bIm = -a.cIm;

	    cRe =  a.bRe;
	    cIm = -a.bIm;

	    dRe =  a.dRe;
	    dIm = -a.dIm;
	}
    }

    /** assigns this with its own conjugate */
    protected void assignConjugate() {
	assignConjugate( this );
    }


    /** assigns this with the transposed of <code>a</code> */
    protected void assignTranspose( final AbstractComplex2By2 a ) {

	if( a == this ) {
	    double tmp;

	    tmp = bRe; bRe = cRe; cRe = tmp;
	    tmp = bIm; bIm = cIm; cIm = tmp;

	} else {

	    aRe = a.aRe;
	    aIm = a.aIm;

	    bRe = a.cRe;
	    bIm = a.cIm;

	    cRe = a.bRe;
	    cIm = a.bIm;

	    dRe = a.dRe;
	    dIm = a.dIm;
	}
    }

    /** assigns this with its own transposed */
    protected void assignTranspose() {
	assignTranspose( this );
    }


    /** assigns this with the transposed and conjugated of <code>a</code> */
    protected void assignStar( final AbstractComplex2By2 a ) {

	if( a == this ) {
	    double tmp;

	    tmp = bRe; bRe =  cRe; cRe =  tmp;
	    tmp = bIm; bIm = -cIm; cIm = -tmp;

	    aIm = -aIm;
	    dIm = -dIm;

	} else {

	    aRe =  a.aRe;
	    aIm = -a.aIm;

	    bRe =  a.cRe;
	    bIm = -a.cIm;

	    cRe =  a.bRe;
	    cIm = -a.bIm;

	    dRe =  a.dRe;
	    dIm = -a.dIm;
	}
    }

    /** assigns this with its own transposed and conjugated*/
    protected void assignStar() {
	assignStar( this );
    }

    /** assigns this with the adjugate of <code>a</code> */
    protected void assignAdjugate( final AbstractComplex2By2 a ) {

	double tmp;

	tmp = a.aRe; aRe = a.dRe; dRe = tmp;
	tmp = a.aIm; aIm = a.dIm; dIm = tmp;

	bRe = -a.bRe;
	bIm = -a.bIm;

	cRe = -a.cRe;
	cIm = -a.cIm;
    }

    /** assigns this with its own adjugate */
    protected void assignAdjugate() {
	assignAdjugate( this );
    }

    /** assigns this with the inverse of <code>m</code> */
    protected void assignInvert( final AbstractComplex2By2 m ) {

	double detRe = m.aRe*m.dRe - m.aIm*m.dIm - ( m.bRe*m.cRe - m.bIm*m.cIm );
	double detIm = m.aRe*m.dIm + m.aIm*m.dRe - ( m.bRe*m.cIm + m.bIm*m.cRe );

	double nn = detRe * detRe + detIm * detIm;

	if( nn == 0 )
	    throw new RuntimeException( m + " matrix is not invertable" );

	assignAdjugate( m );
	assignTimes( detRe / nn, - detIm / nn );
    }

    /** assigns this with its own inverse */
    protected void assignInvert() {
	assignInvert( this );
    }

    /** assigns this with the adjoined of <code>a</code> */
    protected void assignAdjoined( final AbstractComplex2By2 a ) {

	double tmp;

	tmp =  a.bRe; bRe =  a.cRe; cRe = tmp;
	tmp = -a.bIm; bIm = -a.cIm; cIm = tmp;

	aRe = a.aRe; aIm = -a.aIm;
	dRe = a.dRe; dIm = -a.dIm;
    }

    /** assigns this with its own adjoined */
    protected void assignAdjoined() {
	assignAdjoined( this );
    }

    /** assigns this with <code>m</code> and normalizes (determinant equals 1) it */
    protected void assignNormalizeDeterminant( final AbstractComplex2By2 m ) {

	final double detRe = m.aRe*m.dRe - m.aIm*m.dIm  - ( m.bRe*m.cRe - m.bIm*m.cIm );
	final double detIm = m.aRe*m.dIm + m.aIm*m.dRe  - ( m.bRe*m.cIm + m.bIm*m.cRe );

	if( Math.abs( detIm ) > EPS || Math.abs( detRe - 1.0 ) > EPS ) {

	    final double rr = Math.sqrt( Complex.abs( detRe, detIm ) );
	    final double ii = Complex.arg( detRe, detIm ) / 2;

	    final double sqrtOfDetRe = rr * Math.cos( ii );
	    final double sqrtOfDetIm = rr * Math.sin( ii );

	    final double nn = sqrtOfDetRe*sqrtOfDetRe + sqrtOfDetIm*sqrtOfDetIm;

	    if( nn == 0 )
		throw new ArithmeticException("Matrix determinant is zero" );

	    if( this != m )
		assign(m);

	    assignTimes( sqrtOfDetRe, -sqrtOfDetIm );
	    assignDivide( nn );
	}
    }

    /** normalizes (determinant equals 1) this */
    protected void assignNormalizeDeterminant() {
	assignNormalizeDeterminant( this );
    }

    /**
     * Assign <code>this</code> by column.
     *
     * @param column1 a {@link Complex2}: the first column. Must not be
     * <code>null</code>
     * @param column2 a {@link Complex2}: the second column. Must not be
     * <code>null</code>
     */
    protected void assignByColumn(final Complex2 column1,
                                  final Complex2 column2) {
        assign(column1.aRe, column1.aIm, column2.aRe, column2.aIm,
               column1.bRe, column1.bIm, column2.bRe, column2.bIm);
    }

    /**
     * Assign <code>this</code> with <code>t this t<sup>-1</sup>.
     *
     * @param t an {@link AbstractComplex2By2}: the matrix with which to
     * conjugate
     * @throws IllegalArgumentException if determinant of parameter
     * <code>t</code> is <code>0.0</code>.
     */
    protected void assignConjugateWith(final AbstractComplex2By2 t)
        throws IllegalArgumentException {
        t.determinant(dummyComplex);
        if (dummyComplex.re == 0.0 && dummyComplex.im == 0.0) {
            throw new IllegalArgumentException("Det t == 0.0. ");
        }
        this.assignAdjugate();
        this.assignTimes(t, this);
        this.assignAdjugate();
        this.assignTimes(t, this);
        this.assignDivide(dummyComplex);
    }

    /**
     * Assign <code>this</code> with <code>t this t<sup>*</sup></code>.
     *
     * @param t an {@link AbstractComplex2By2}: the matrix to adjoin
     */
    protected void assignAdjoinedWith(final AbstractComplex2By2 t) {
        this.assignTimes(t, this);
        this.assignAdjoined();
        this.assignTimes(t, this);
        this.assignAdjoined();
    }

    /**
     * <p>Assign <code>this</code> by prescribing eigenvectors and
     * eigenvalues. </p>
     *
     * <p> The given eigenvectors must be linearly independent. </p>
     *
     * @param eigenvalue1Re a <code>double</code>: real part of the first
     * eigenvalue
     * @param eigenvalue1Im a <code>double</code>: imaginary part of the
     * first eigenvalue
     * @param eigenvector1 a {@link Complex2}: first eigenvector
     * @param eigenvalue2Re a <code>double</code>: real part of the second
     * eigenvalue
     * @param eigenvalue2Im a <code>double</code>: imaginary part of the
     * second eigenvalue
     * @param eigenvector2 a {@link Complex2}: second eigenvector
     */
    protected void assignByEigenvectors(final double eigenvalue1Re,
                                        final double eigenvalue1Im,
                                        final Complex2 eigenvector1,
                                        final double eigenvalue2Re,
                                        final double eigenvalue2Im,
                                        final Complex2 eigenvector2) {

        /* Let this be a diagonal matrix with the eigenvalues as entries. */
        this.assign(eigenvalue1Re, eigenvalue1Im, 0.0, 0.0,
                    0.0, 0.0, eigenvalue2Re, eigenvalue2Im);

        /* Make a matrix whose columns are the eigenvectors. */
        dummyComplex2By2.assignByColumn(eigenvector1, eigenvector2);

        /* Transform this. */
        this.assignConjugateWith(dummyComplex2By2);
    }

    /**
     * <p>Assign <code>this</code> by prescribing eigenvectors and
     * eigenvalues. </p>
     *
     * <p> The given eigenvectors must be linearly independent. </p>
     *
     * @param eigenvalue1 a {@link Field.Complex}: first eigenvalue
     * @param eigenvector1 a {@link Complex2}: first eigenvector
     * @param eigenvalue2 a {@link Field.Complex}: second eigenvalue
     * @param eigenvector2 a {@link Complex2}: second eigenvector
     */
    protected void assignByEigenvectors(final Complex eigenvalue1,
                                        final Complex2 eigenvector1,
                                        final Complex eigenvalue2,
                                        final Complex2 eigenvector2) {
        assignByEigenvectors(eigenvalue1.re, eigenvalue1.im, eigenvector1,
                             eigenvalue2.re, eigenvalue2.im, eigenvector2);
    }

    /** returns square of eucleaden norm */
    final public double absSqr() {

	return aRe*aRe + aIm*aIm + bRe*bRe + bIm*bIm + cRe*cRe + cIm*cIm +dRe*dRe + dIm*dIm;
    }

    /** returns eucleaden norm */
    final public double abs() {
	return Math.sqrt( absSqr() );
    }

    /** returns square of eucleaden distance of this to <code>s</code> */
    final public double distSqr( final AbstractComplex2By2 s ) {

	return
	    (s.aRe-aRe)*(s.aRe-aRe) + (s.aIm-aIm)*(s.aIm-aIm) +
	    (s.bRe-bRe)*(s.bRe-bRe) + (s.bIm-bIm)*(s.bIm-bIm) +
	    (s.cRe-cRe)*(s.cRe-cRe) + (s.cIm-cIm)*(s.cIm-cIm) +
	    (s.dRe-dRe)*(s.dRe-dRe) + (s.dIm-dIm)*(s.dIm-dIm);
    }

    /** returns eucleaden distance of this to <code>s</code> */
    final public double dist( final AbstractComplex2By2 s ) {
	return Math.sqrt( distSqr( s ) );
    }

    /** returns determinant of this */
    final public Complex determinant() {
	Complex det = new Complex();
	determinant( det );
	return det;
    }

    /** assigns <code>det</code> with the determinant of this
     ** @param det complex value which is assigned with the determinant */
    final public void determinant( final Complex det ) {

	final double detRe = aRe*dRe - aIm*dIm - bRe*cRe + bIm*cIm;
	final double detIm = aRe*dIm + aIm*dRe - bRe*cIm - bIm*cRe;

	det.assign( detRe, detIm );
    }


    /** returns trace of this */
    final public Complex trace() {
	Complex trace = new Complex();
	trace( trace );
	return trace;
    }

    /** assigns <code>trace</code> with the trace of this
     ** @param trace complex value which is assigned with the trace */
    final public void trace( final Complex trace ) {

	trace.assign( aRe + dRe, aIm + dIm );
    }

    /** returns the square of absulate value of the bigger eigenvalue.
     ** It creates a temporarily used instance of <code>Field.Complex</code> */
    final public double normSqr() {
	final Complex bigEV = new Complex();
	getEigenValues( null, bigEV );
	return bigEV.absSqr();
    }

    /** returns the absulate value of the bigger eigenvalue.
     ** It creates one temporarily used instance of <code>Field.Complex</code> */
    final public double norm() {
	final Complex bigEV = new Complex();
	getEigenValues( null, bigEV );
	return bigEV.abs();
    }

    /** assigns <code>smallEV</code> with the smaller eigenvalue of this
     ** and <code>bigEV</code> with the bigger;
     ** both, <code>smallEV</code> and <code>bigEV</code>, may be <code>null</code>. */
    final public void getEigenValues( final Complex smallEV, final Complex bigEV ) {

	final double detRe = aRe*dRe - aIm*dIm - bRe*cRe + bIm*cIm;
	final double detIm = aRe*dIm + aIm*dRe - bRe*cIm - bIm*cRe;

	final double aPlusDRe = aRe + dRe;
	final double aPlusDIm = aIm + dIm;

	// dis = (a+d)^2 - 4*det
	final double disRe = aPlusDRe * aPlusDRe - aPlusDIm * aPlusDIm - 4*detRe;
	final double disIm = aPlusDRe * aPlusDIm + aPlusDIm * aPlusDRe - 4*detIm;

	// compute sqrt of dis
	final double rr = Math.sqrt( Complex.abs( disRe, disIm ) );
	final double ii = Complex.arg( disRe, disIm ) / 2;

	double sqrtOfDisRe = rr * Math.cos( ii );
	double sqrtOfDisIm = rr * Math.sin( ii );

	// we want < a+b, dis^0.5 > to be non negative
	if( aPlusDRe * sqrtOfDisRe + aPlusDIm * sqrtOfDisIm < 0 ) {
	    sqrtOfDisRe = -sqrtOfDisRe;
	    sqrtOfDisIm = -sqrtOfDisIm;
	}

	if( smallEV != null )
	    smallEV.assign( ( aPlusDRe - sqrtOfDisRe ) / 2,
			    ( aPlusDIm - sqrtOfDisIm ) / 2 );

	if( bigEV != null )
	    bigEV.  assign( ( aPlusDRe + sqrtOfDisRe ) / 2,
			    ( aPlusDIm + sqrtOfDisIm ) / 2 );
    }

    /** assigns the first two entries of array <code>ev</code> to the
     ** eigenvalues of the this;
     ** the first entry contains the eigenvalue with the smaller euclidean norm.
     ** if an entry of the array is <code>null</code> an instance of <code>Field.Complex</code> is created. */
    final public void getEigenValues( final Complex[] ev ) {

	if( ev[0] == null ) ev[0] = new Complex();
	if( ev[1] == null ) ev[1] = new Complex();

	getEigenValues( ev[0], ev[1] );
    }

    /** returns array containing the two eigenvalues;
     ** the first entry contains the eigenvalue with the smaller euclidean norm. */
    final public Complex[] getEigenValues() {
	final Complex[] ev = new Complex[2];
	getEigenValues( ev );
	return ev;
    }

    /** this is considered to equal <code>o</code> if
     ** their euclidean distance is smaller than <code>{@link #EPS}</code> */
    public boolean equals(Object o) {
	try {
	    return distSqr( (AbstractComplex2By2)o ) < EPS * EPS;

	} catch(ClassCastException ex) {

	    return false;
	}
    }

    public String toString() {
	StringBuffer sb=new StringBuffer().append('(').append('(').append('(').append(aRe);
	if(aIm>=0.0)
	    sb.append('+');
	sb.append(aIm).append('i').append(')').append(',').append('(').append(bRe);
	if(bIm>=0.0)
	    sb.append('+');
	sb.append(bIm).append('i').append(')').append(')').append(',');
	sb.append('(').append('(').append(cRe);
	if(cIm>=0.0)
	    sb.append('+');
	sb.append(cIm).append('i').append(')').append(',').append('(').append(dRe);
	if(dIm>=0.0)
	    sb.append('+');
	sb.append(dIm).append('i').append(')').append(')').append(')');

	return sb.toString();
    }

    /**
     * Multiplies <code>m</code> with <code>v</code> and stores the result
     * in <code>w</code>.
     */
    public static void times( AbstractComplex2By2 m, Complex2 v, Complex2 w ) {

	if ( v == w ) {
	    double wARe = m.aRe*v.aRe-m.aIm*v.aIm + m.bRe*v.bRe-m.bIm*v.bIm;
	    double wAIm = m.aRe*v.aIm+m.aIm*v.aRe + m.bRe*v.bIm+m.bIm*v.bRe;
	    double wBRe = m.cRe*v.aRe-m.cIm*v.aIm + m.dRe*v.bRe-m.dIm*v.bIm;
	    double wBIm = m.cRe*v.aIm+m.cIm*v.aRe + m.dRe*v.bIm+m.dIm*v.bRe;

	    w.aRe=wARe; w.aIm=wAIm; w.bRe=wBRe; w.bIm=wBIm;
	} else {
	    w.aRe = m.aRe*v.aRe-m.aIm*v.aIm + m.bRe*v.bRe-m.bIm*v.bIm;
	    w.aIm = m.aRe*v.aIm+m.aIm*v.aRe + m.bRe*v.bIm+m.bIm*v.bRe;
	    w.bRe = m.cRe*v.aRe-m.cIm*v.aIm + m.dRe*v.bRe-m.dIm*v.bIm;
	    w.bIm = m.cRe*v.aIm+m.cIm*v.aRe + m.dRe*v.bIm+m.dIm*v.bRe;
	}
    }
}














