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

package de.jtem.mfc.field;

import java.io.Serializable;
import java.util.LinkedList;
import java.util.ListIterator;

/** This class provides a comprehensive library for complex numbers.
  * The underlying ideas are the following:
  * <p>
  * <em>Method names</em> refer to the action performed:
  * <br>
  *    For instances <code>a,b</code> of class <code>Field.Complex</code>, <code>a.plus(b)</code>
  *    produces a new instance of <code>Field.Complex</code> containing the sum <code>a+b</code>.
  * <p>
  * <em>Assign-methods</em> do NOT create new instances:
  * <br>
  *    For instances <code>a,b,c</code> of class <code>Field.Complex</code>,
  *    <code>c.assignPlus(a,b)</code> stores the sum <code>a+b</code> in <code>c</code>.
  *
  *    <code>a.assignPlus(b)</code> stores the sum <code>a+b</code> in <code>a</code>,
  *    like the operator "+=" does for simple types.
  * <p>
  * Field.Complex numbers INTERACT WITH DOUBLES:
  * <br>
  *   For an instance <code>a</code> of class <code>Field.Complex</code> and a double <code>r</code>,
  *   <code>a.times(r)</code> creates a new Field.Complex containing the product <code>a*r</code>.
  *
  * @author schmies
  */

public class ComplexConstant extends ComplexValue implements Field.Complex, Serializable, Cloneable {

    private static final long serialVersionUID = 1L;

    /** 0 */
    public static final ComplexConstant ZERO       =  new ComplexConstant( 0.0, 0.0 );
    /** 1 */
    public static final ComplexConstant ONE        =  new ComplexConstant( 1.0, 0.0 );
    /** i */
    public static final ComplexConstant I          =  new ComplexConstant( 0.0, 1.0 );
    /** -1 */
    public static final ComplexConstant NEG_ONE    =  new ComplexConstant(-1.0, 0.0 );
    /** -i */
    public static final ComplexConstant NEG_I      =  new ComplexConstant( 0.0,-1.0 );

    /** real part of complex number */
    double re;
    /** imaginary part of complex number */
    double im;

    /** Creates complex zero. */
    public ComplexConstant () {
      re = 0.;
      im = 0.;
    }
    /** Creates <code>this</code> with <code>re = aReal</code>.
     * @param aReal is a real part of <code>this</code>
     */
    public ComplexConstant (final double aReal) {
      re = aReal;
      im = 0.;
    }

    /** Creates <code>this</code> which is equal to u */
    public ComplexConstant(final ComplexConstant u) {
      re = u.re;
      im = u.im;
    }

    /** Creates <code>this</code> which is equal to u */
    public ComplexConstant(final Field.Complex u) {
      re = u.getRe();
      im = u.getIm();
    }

    /** Creates <code>this</code> with real and imag entry.
     * @param real and imag. part of <code>this</code>
     */
    public ComplexConstant (final double aReal,final double aImag ) {
      re =  aReal;
      im =  aImag;
    }
    /** Return real part of <code>this</code>.
     * @return a real part
     */
    public final double getRe() {
      return re;
    }
   /** Return imag. part of <code>this</code>.
     * @return a imag. part
     */
    public final double getIm() {
      return im;
    }

    /** Assigns <code>this</code> with real and imag entry.
     * @param real and imag. part
     */
    final void assign( final double aReal, final double aImag) {
      re =  aReal;
      im =  aImag;
    }

   /** Assigns <code>this</code> with <code>v</code>.
    * Assigns and returns instance of <code>this</code>.
    * @return <code>this</code>
    */
    final void assign(ComplexConstant v) {
      re =  v.re;
      im =  v.im;
    }

    /** Assigns <code>this</code> with <code>v</code>.
     * Assigns and returns instance of <code>this</code>.
     * @return <code>this</code>
     */
    final void assign(Field.Complex v) {
      re = v.getRe();
      im = v.getIm();
    }

   /** Assigns <code>this</code> with <code>re</code>.
    * Assigns and returns instance of <code>this</code>.
    * @return <code>this</code>
    */
    final void assign(double re) {
      this.re =  re;
      im =  0.;
    }

    /** Creates a new <code>Field.Complex</code> that is copy of <code>this</code>.
     * @return a new Field.Complex
     */
    public final ComplexConstant copy() {
      return new ComplexConstant(re,im);
    }
    /** Clones <code>this</code> to new <code>Field.Complex</code>.
     * @see #copy()
     * @return copy of this as <code>Object</code>
     */
    public final Object clone() {
      return copy();
    }
    /**  Returns <tt>true</tt> if <code>this</code> has any infinite entries.
     * @return <tt>true</tt> if infinite, <tt>false</tt> else
     */
    public final boolean isInfinite () {
      return (Double.isInfinite(re) || Double.isInfinite(im));
    }
    /** Returns <tt>true</tt> if <code>this</code> has any Not-a-Number(NaN) entries
     * @return <tt>true</tt> if is NaN, <tt>false</tt> otherwise
     */
    public final boolean isNaN () {
      return (Double.isNaN(re) || Double.isNaN(im));
    }

    /** Returns <tt>true</tt> if <code>this</code> has only zero entries.
     * @return <tt>true</tt> if <code>this</code> is zero, <tt>false</tt> otherwise
     */
    public boolean isZero() {
      return equals(0,0);
    }

    /** Returns sum of this and <code>v</code>
     * Creates new <code>Field.Complex</code> w = this + v .
     */
    public final ComplexConstant plus(final Field.Complex v) {
      return new ComplexConstant(re + v.getRe(), im + v.getIm());
    }

    /** Creates a new <code>Field.Complex</code> C  which is sum of <code>this</code> and <code>r</code>.
     * @return C = this + r.
     * @param summand <code>r</code> as a double
     */
    public final ComplexConstant plus( final double r ) {
      return new ComplexConstant(re + r, im);
    }
    /** Assigns <code>this</code> with sum of <code>this</code> and <code>x,y</code>.
     *	re += x
     *	im += y
     * @param summands <code>x</code> and <code>y</code> as a double
     */
    final void assignPlus( final double x, final double y) {
	re += x;
	im += y;
    }

    /** Assigns <code>this</code> with  sum of <code>this</code> and <code>v</code>.
     * @param summand <code>v</code> as a Field.Complex
     */
    final void assignPlus( final ComplexConstant v ) {
      re += v.re;
      im += v.im;
    }

    /** Assigns <code>this</code> with  sum of <code>this</code> and <code>v</code>.
     * @param summand <code>v</code> as a Field.Complex
     */
    final void assignPlus( final Field.Complex v ) {
      re += v.getRe();
      im += v.getIm();
    }

    /**  Assigns <code>this</code> with sum of <code>u</code> and <code>r</code>.
     * @param summands <code>u</code> as Field.Complex and <code>r</code> as double
     */
    final void assignPlus( final ComplexConstant u, final double r ) {
      re = u.re + r;
      im = u.im;
    }

    /**  Assigns <code>this</code> with sum of <code>u</code> and <code>r</code>.
     * @param summands <code>u</code> as Field.Complex and <code>r</code> as double
     */
    final void assignPlus( final Field.Complex u, final double r ) {
      re = u.getRe() + r;
      im = u.getIm();
    }

    /** Assigns <code>this</code> with  sum of <code>this</code> and <code>r</code>.
     * @param summand <code>r</code> as a double
     */
    final void assignPlus( final double r ) {
      re +=r;
    }

    /**  Assigns <code>this</code> with sum of <code>u</code> and <code>v</code>.
     * @param summands <code>u</code> and <code>v</code> as Field.Complex
     */
    final void assignPlus( final ComplexConstant u, final ComplexConstant v) {
      re = u.re + v.re;
      im = u.im + v.im;
    }

    /**  Assigns <code>this</code> with sum of <code>u</code> and <code>v</code>.
     * @param summands <code>u</code> and <code>v</code> as Field.Complex
     */
    final void assignPlus( final Field.Complex u, final Field.Complex v) {
      re = u.getRe() + v.getRe();
      im = u.getIm() + v.getIm();
    }

    /** Creates a new <code>Field.Complex</code> C which is product of substraction <code>v</code> from <code>this</code>.
     * @param  Field.Complex <code>v</code> as subtrahend
     * @return C = this - v.
     */
    public final ComplexConstant minus(final Field.Complex v) {
      return new ComplexConstant(re - v.getRe(), im - v.getIm());
    }

    /** creates a new <code>Field.Complex</code> which is product of substraction <code>r</code> from <code>this</code>.
     * @param  double <code>r</code> as subtrahend
     * @return C = this - r.
     */
    public final ComplexConstant minus( final double r ) {
      return new ComplexConstant(re - r, im);
    }

    /**  Assigns <code>this</code> with product of substraction between <code>this</code> and <code>x,y</code>.
     * this = this - (x+iy).
     * @param <code>x</code> and <code>y</code> as double
     */
    final void assignMinus( final double x, final double y) {
	re -= x;
	im -= y;
    }

    /**  Assigns <code>this</code> with product of substraction between <code>this</code> and <code>v</code>.
     * this -= v.
     * @param Field.Complex <code>v</code> as subtrahend
     */
    final void assignMinus( final ComplexConstant v ) {
      re -= v.re;
      im -= v.im;
    }

    /**  Assigns <code>this</code> with product of substraction between <code>this</code> and <code>v</code>.
     * this -= v.
     * @param Field.Complex <code>v</code> as subtrahend
     */
    final void assignMinus( final Field.Complex v ) {
      re -= v.getRe();
      im -= v.getIm();
    }

    /**  Assigns <code>this</code> with product of substraction between <code>this</code> and <code>r</code>.
     *  this.re -= r.
     * @param double <code>r</code> as subtrahend
     */
    final void assignMinus( final double r ) {
      re -=r;
    }

    /**  Assigns <code>this</code> with product of substraction between t<code>u</code> and <code>v</code>.
     * this = u - v
     * @param <code>u</code> and subtrahend <code>v</code>
     */
    final void assignMinus( final ComplexConstant u, final ComplexConstant v ) {
      re = u.re - v.re;
      im = u.im - v.im;
    }

    /**  Assigns <code>this</code> with product of substraction between t<code>u</code> and <code>v</code>.
     * this = u - v
     * @param <code>u</code> and subtrahend <code>v</code>
     */
    final void assignMinus( final Field.Complex u, final Field.Complex v ) {
      re = u.getRe() - v.getRe();
      im = u.getIm() - v.getIm();
    }

    /**  Assigns <code>this</code> with product of substraction between <code>u</code> and <code>r</code>.
     * this = u - r
     * @param Field.Complex <code>u</code> and subtrahend <code>r</code> as double
     */
    final void assignMinus( final ComplexConstant u, final double r ) {
      re = u.re - r;
      im = u.im;
    }

    /**  Assigns <code>this</code> with product of substraction between <code>u</code> and <code>r</code>.
     * this = u - r
     * @param Field.Complex <code>u</code> and subtrahend <code>r</code> as double
     */
    final void assignMinus( final Field.Complex u, final double r ) {
      re = u.getRe() - r;
      im = u.getIm();
    }

    /**  Assigns <code>this</code> with product of substraction between <code>r</code> and <code>u</code>.
     * this =  r - u, this.im = - u.im.
     * @param double <code>r</code> and subtrahend <code>u</code> as double.
     */
    final void assignMinus( final double r, final ComplexConstant u ) {
      re = r - u.re;
      im = -u.im;
    }

    /**  Assigns <code>this</code> with product of substraction between <code>r</code> and <code>u</code>.
     * this =  r - u, this.im = - u.im.
     * @param double <code>r</code> and subtrahend <code>u</code> as double.
     */
    final void assignMinus( final double r, final Field.Complex u ) {
      re = r - u.getRe();
      im = -u.getIm();
    }

    /** Creates a new <code>Field.Complex</code> C which is product of mltiplication between <code>this</code> with <code>v</code>.
     * @param Field.Complex <code>v</code> as multiplicand
     * @return Field.Complex <code>w</code> = <code>this</code> * <code>v</code>
     */
    public final ComplexConstant times( final Field.Complex v ) {
        return new ComplexConstant(re*v.getRe() - im*v.getIm(), re*v.getIm() + im*v.getRe());
    }

    /** Creates a new <code>Field.Complex</code> C which is product of mltiplication between <code>this</code> with <code>r</code>.
     * @param double <code>r</code> as multiplicand
     * @return Field.Complex <code>w</code> = <code>this</code> * <code>r</code>. */
    public final ComplexConstant times( final double r ) {
	return new ComplexConstant(re * r, im * r);
    }

    /**  Assigns <code>this</code> with  product of mltiplication between <code>this</code> with <code>x,y</code>.
     * @param double <code>x</code> and <code>y</code> as multiplicand
     */
    final void assignTimes( final double x, final double y) {
	final double rr = re*x - im*y;
	final double ii = re*y + im*x;
	re = rr;
	im = ii;
    }

    /**  Assigns <code>this</code> with product of mltiplication between <code>this</code> with <code>v</code>
     * @param Field.Complex <code>v</code> as multiplicand
     */
    final void assignTimes( final ComplexConstant v ) {
      double rr = re*v.re - im*v.im;
      double ii = re*v.im + im*v.re;
      re = rr;
      im = ii;
    }

    /**  Assigns <code>this</code> with product of mltiplication between <code>this</code> with <code>v</code>
     * @param Field.Complex <code>v</code> as multiplicand
     */
    final void assignTimes( final Field.Complex v ) {
      double rr = re*v.getRe() - im*v.getIm();
      double ii = re*v.getIm() + im*v.getRe();
      re = rr;
      im = ii;
    }

    /**  Assigns <code>this</code> with product of mltiplication between <code>this</code> with <code>r</code>
     * @param double <code>r</code> as multiplicand
     */
    final void assignTimes( final double r ) {
      re *=r;
      im *=r;
    }

    /**  Assigns <code>this</code> with product of between <code>x+iy</code>
     *   with <code>u+iv</code>.
     * this = (x+iy)(u+iv)
     * @param <code>u</code> as multiplicator and <code>v</code> as multiplicand
     */
    final void assignTimes( final double x, final double y, final double u, final double v ) {
      re = x*u - y*v;
      im = x*v + y*u;
    }

    /**  Assigns <code>this</code> with product of mltiplication between <code>u</code> with <code>v</code>.
     * this = u * v
     * @param <code>u</code> as multiplicator and <code>v</code> as multiplicand
     */
    final void assignTimes( final ComplexConstant u, final ComplexConstant v ) {
	double rr = u.re*v.re - u.im*v.im;
	double ii = u.re*v.im + u.im*v.re;
	re = rr;
	im = ii;
    }

    /**  Assigns <code>this</code> with product of mltiplication between <code>u</code> with <code>v</code>.
     * this = u * v
     * @param <code>u</code> as multiplicator and <code>v</code> as multiplicand
     */
    final void assignTimes( final Field.Complex u, final Field.Complex v ) {
      assignTimes( u.getRe(), u.getIm(), v.getRe(), v.getIm() );
    }

    /**  Assigns <code>this</code> with product of mltiplication between <code>u</code> with <code>r</code>.
     * this = u * r.
     * @param <code>u</code> as multiplicator and <code>r</code> as multiplicand
     */
    final void assignTimes( final ComplexConstant u, final double r ) {
	re = u.re*r;
	im = u.im*r;
    }

    /**  Assigns <code>this</code> with product of mltiplication between <code>u</code> with <code>r</code>.
     * this = u * r.
     * @param <code>u</code> as multiplicator and <code>r</code> as multiplicand
     */
    final void assignTimes( final Field.Complex u, final double r ) {
        re = u.getRe()*r;
        im = u.getIm()*r;
    }

    /** Returns new <code>Field.Complex</code> which is assigned with the
     * the product of <code>this</code> and <code>i</code>.
     * @return <code>i * this</code>
     */
    public final ComplexConstant timesI() {
	return new ComplexConstant( -im, re );
    }

    /** Multiplies <code>this</code> with <code>i</code>. */
    final void assignTimesI() {
	double rr = re;
	re = -im;
	im =  rr;
    }

    /** Assigns <code>this</code> with the product of <code>v</code> and <code>i</code>. */
    final void assignTimesI( final ComplexConstant v ) {
	double rr = v.re;
	re = -v.im;
	im =  rr;
    }

    /** Assigns <code>this</code> with the product of <code>v</code> and <code>i</code>. */
    final void assignTimesI( final Field.Complex v ) {
        double rr = v.getRe();
        re = -v.getIm();
        im =  rr;
    }

    /** Returns new <code>Field.Complex</code> which is assigned with the
     * the quotient of <code>this</code> and <code>i</code>.
     * @return <code>i * this</code>
     */
    public final ComplexConstant divideI() {
	return new ComplexConstant( im, -re );
    }

    /** Divides <code>this</code> with <code>i</code>. */
    final void assignDivideI() {
	double rr = re;
	re =  im;
	im = -rr;
    }

    /** Assigns <code>this</code> with the quotient of <code>v</code> and <code>i</code>. */
    final void assignDivideI( final ComplexConstant v ) {
	double rr = v.re;
	re = v.im;
	im =  -rr;
    }

    /** Assigns <code>this</code> with the quotient of <code>v</code> and <code>i</code>. */
    final void assignDivideI( final Field.Complex v ) {
        double rr = v.getRe();
        re = v.getIm();
        im =  -rr;
    }

    /**
     * Returns ratio of <code>this</code> and <code>v</code>.
     * Creates new <code>Field.Complex</code> for the result.
     * If <code>v</code> is zero then the components of the
     * result equal plus or minus infinity or zero depending on
     * the sign of the related component of <code>this</code>.
     * @param  <code>v</code> the divisor
     * @return <code>this/v</code> */
    public final ComplexConstant divide( final Field.Complex v ) {
        final double vRe=v.getRe();
        final double vIm=v.getIm();
        double nn = vRe*vRe + vIm*vIm;
        if (nn == 0.) {
            return new ComplexConstant( re/nn,im/nn );
        }
        return new ComplexConstant( (re*vRe + im*vIm)/nn, (im*vRe - re*vIm)/nn );
    }

    /** Returns ratio of <code>this</code> and <code>r</code>.
     * Creates new <code>Field.Complex</code> for the result.
     * If <code>r</code> is zero then the components of the
     * result equal plus or minus infinity or zero depending on
     * the sign of the related component of <code>this</code>.
     * @param <code>r</code> as divisor
     * @return <code>this/r</code>.
     */
    public final ComplexConstant divide( final double r) {
      return new ComplexConstant(re/r, im/r);
    }

    /** Assigns <code>this</code> with product of division between <code>this</code> and <code>x,y</code>: this /= x+i*y.
     * If <code>x,y</code> is zero then the components of the
     * result equal plus or minus infinity or zero depending on
     * the sign of the related component of <code>this</code>.
     * @param <code>x</code> and <code>y</code> as divisors
     */
    final void assignDivide( final double x, final double y) {
	final double nn = x*x + y*y;
	if (nn == 0.) {
	    re = re/nn;
	    im = im/nn;
	} else {
	    final double rr = re;
	    final double ii = im;
	    re = (rr*x + ii*y)/nn;
	    im = (ii*x - rr*y)/nn;
	}
    }

    /** Assigns <code>this</code> with product of division between <code>this</code> and <code>v</code>.
     *	<code>this = this / v</code>
     * If <code>v</code> is zero then the components of the
     * result equal plus or minus infinity or zero depending on
     * the sign of the related component of <code>this</code>.
     *  @param <code>v</code> as divisor
     */
    final void assignDivide( final ComplexConstant v ) {
	assignDivide( v.re, v.im );
    }

    /** Assigns <code>this</code> with product of division between <code>this</code> and <code>v</code>.
     *	<code>this = this / v</code>
     * If <code>v</code> is zero then the components of the
     * result equal plus or minus infinity or zero depending on
     * the sign of the related component of <code>this</code>.
     *  @param <code>v</code> as divisor
     */
    final void assignDivide( final Field.Complex v ) {
        assignDivide( v.getRe(), v.getIm() );
    }

    /**  Assigns <code>this</code> with product of division between <code>this</code> and <code>r</code>:this /= r.
     * If <code>r</code> is zero then the components of the
     * result equal plus or minus infinity or zero depending on
     * the sign of the related component of <code>this</code>.
     * @param <code>r</code> as divisor
     */
    final void assignDivide( final double r ) {
      re /= r;
      im /= r;
    }


    /**  assigns <code>this</code> with product of division between
     * <code>x+iy</code> and <code>u+iv</code>: this = (x+iy) / (u+iv).
     * If <code>v</code> is zero then the components of the
     * result equal plus or minus infinity or zero depending on
     * the sign of the related component of <code>u</code>.
     * @param <code>u</code> as dividend and <code>v</code> as divisor
     */
    final void assignDivide(final double x, final double y, final double u, final double v) {
      final double nn = u*u + v*v;
      if (nn == 0.) {
        re = re / nn;
        im = im / nn;
      }
      else {
        double rr = (x * u + y * v) / nn;
        double ii = (y * u - x * v) / nn;
        re = rr;
        im = ii;
      }
    }

    /**  assigns <code>this</code> with product of division between <code>u</code> and <code>v</code>: this = u / v.
     * If <code>v</code> is zero then the components of the
     * result equal plus or minus infinity or zero depending on
     * the sign of the related component of <code>u</code>.
     * @param <code>u</code> as dividend and <code>v</code> as divisor
     */
    final void assignDivide(final ComplexConstant u,
                                   final ComplexConstant v) {
      final double nn = v.re * v.re + v.im * v.im;
      if (nn == 0.) {
        re = re / nn;
        im = im / nn;
      }
      else {
        double rr = (u.re * v.re + u.im * v.im) / nn;
        double ii = (u.im * v.re - u.re * v.im) / nn;
        re = rr;
        im = ii;
      }
    }

    /**  assigns <code>this</code> with product of division between <code>u</code> and <code>v</code>: this = u / v.
     * If <code>v</code> is zero then the components of the
     * result equal plus or minus infinity or zero depending on
     * the sign of the related component of <code>u</code>.
     * @param <code>u</code> as dividend and <code>v</code> as divisor
     */
    final void assignDivide( final Field.Complex u, final Field.Complex v ) {
      assignDivide( u.getRe(), u.getIm(), v.getRe(), v.getIm() );
    }

    /**  assigns <code>this</code> with product of division between <code>u</code> and <code>r</code>: this = u / r.
     * If <code>r</code> is zero then the components of the
     * result equal plus or minus infinity or zero depending on
     * the sign of the related component of <code>u</code>.
     * @param <code>u</code> as dividend and <code>r</code> as divisor
     */
    final void assignDivide( final ComplexConstant u, final double r ) {
      	re = u.re / r;
	im = u.im / r;
    }

    /**  assigns <code>this</code> with product of division between <code>u</code> and <code>r</code>: this = u / r.
     * If <code>r</code> is zero then the components of the
     * result equal plus or minus infinity or zero depending on
     * the sign of the related component of <code>u</code>.
     * @param <code>u</code> as dividend and <code>r</code> as divisor
     */
    final void assignDivide(final Field.Complex u, final double r) {
      re = u.getRe() / r;
      im = u.getIm() / r;
    }

    /**  assigns <code>this</code> with product of division between <code>r</code> and <code>u</code>: this = r / u.
     * If <code>u</code> is zero then the components of the
     * result equal plus or minus infinity or zero depending on
     * the sign of the related component of <code>r</code>.
     * @param <code>r</code> as dividend and <code>z</code> as divisor
     */
    final void assignDivide(final double r, final ComplexConstant z) {
      final double nn = absSqr(z);
      if (nn == 0.) {
        re = r / nn;
        im = 0;
      }
      else {
        re = r * z.re / nn;
        im = -r * z.im / nn;
      }
    }

    /**  assigns <code>this</code> with product of division between <code>r</code> and <code>u</code>: this = r / u.
     * If <code>u</code> is zero then the components of the
     * result equal plus or minus infinity or zero depending on
     * the sign of the related component of <code>r</code>.
     * @param <code>r</code> as dividend and <code>z</code> as divisor
     */
    final void assignDivide(final double r, final Field.Complex z) {
      final double u = z.getRe();
      final double v = z.getIm();
      final double nn = u * u + v * v;
      if (nn == 0.) {
        re = r / nn;
        im = 0;
      }
      else {
        re = r * u / nn;
        im = -r * v / nn;
      }
    }

    /** Returns invert of <code>this</code>.
     * Creates new <code>Field.Complex</code> for the result.
     * In case of |this|==0 the following result occurs:
     * this.re = Double.POSITIVE_INFINITY, this.im = 0.
     */
    public final ComplexConstant invert() {
      double nn = re*re + im*im;
      if (nn == 0.) {
        return new ComplexConstant(Double.POSITIVE_INFINITY, 0);
      }
      return new ComplexConstant(re/nn, -im/nn);
    }

    /** Returns invert of <code>u</code>.
     * Creates new <code>Field.Complex</code> that is invert of <code>u</code>.
     * In case of |u|==0 the following result occurs:
     * u.re = Double.POSITIVE_INFINITY, u.im = 0.
     */
    public final static ComplexConstant invert( final Field.Complex u ) {
      final double uRe=u.getRe();
      final double uIm=u.getIm();
      double nn = uRe*uRe + uIm*uIm;
      if (nn == 0.) {
        return new ComplexConstant(Double.POSITIVE_INFINITY, 0);
      }
      return new ComplexConstant(uRe/nn, -uIm/nn);
    }

    /** u = 1./u.
     * In case of |u|==0 the following result occurs:
     * u.re = Double.POSITIVE_INFINITY, u.im = 0.
     */
    final void assignInvert() {
      double nn = re*re +im*im;
      if (nn == 0.) {
          re = Double.POSITIVE_INFINITY;
          im = 0;
          return;
      }
      re /= nn;
      im /= -nn;
    }

     /** Assigns <code>this</code> with inverse of <code>u</code>.
      * this=1/u.
      * In case of |u|==0 the following result occurs:
      * this.re = Double.POSITIVE_INFINITY, this.im = 0
      */
    final void assignInvert(ComplexConstant u) {
      double nn = u.re*u.re + u.im*u.im;
      if (nn == 0.) {
          re = Double.POSITIVE_INFINITY;
          im = 0;
      } else {
	  re = u.re/nn;
	  im = -u.im/nn;
      }
    }

    /** Assigns <code>this</code> with inverse of <code>u</code>.
     * this=1/u.
     * In case of |u|==0 the following result occurs:
     * this.re = Double.POSITIVE_INFINITY, this.im = 0
     */
   final void assignInvert(Field.Complex u) {
     final double uRe=u.getRe();
     final double uIm=u.getIm();
     double nn = uRe*uRe + uIm*uIm;
     if (nn == 0.) {
         re = Double.POSITIVE_INFINITY;
         im = 0;
     } else {
         re = uRe/nn;
         im = -uIm/nn;
     }
   }

    /** Returns conjugation of <code>u</code>.
     *	Creates new <code>Field.Complex</code> that conjugates <code>u</code>.
     * @param <code>u</code> as Field.Complex
     */
    public final static ComplexConstant conjugate( final Field.Complex u ) {
        return new ComplexConstant( u.getRe(), -u.getIm() );
    }

    /** Returns conjugation of <code>this</code>.
     *	Creates new <code>Field.Complex</code> that conjugates <code>this</code>.
     */
    public final ComplexConstant conjugate() {
	return new ComplexConstant( re, -im );
    }

    /** Conjugates <code>this</code>.
     * <code>this.im = -this.im</code>
     */
    final void assignConjugate() {
	im = -im;
    }

   /** Assigns <code>this</code> with conjugation of <code>u</code>.
    *	<code>this.re =  u.re</code>
    *	<code>this.im = -u.im</code>
    */
    final void assignConjugate( final ComplexConstant u ) {
	re =  u.re;
	im = -u.im;
    }

    /** Assigns <code>this</code> with conjugation of <code>u</code>.
     *	<code>this.re =  u.re</code>
     *	<code>this.im = -u.im</code>
     */
     final void assignConjugate( final Field.Complex u ) {
         re =  u.getRe();
         im = -u.getIm();
     }

    /** Returns negative of <code>u</code>.
     *	Creates new <code>Field.Complex</code> that equals -<code>u</code>.
     */
    public final static ComplexConstant neg( final Field.Complex u ) {
        return new ComplexConstant( -u.getRe(), -u.getIm() );
    }

    /** Returns negative of <code>this</code>.
     *	Creates new <code>Field.Complex</code> that equals -<code>this</code>.
     */
    public final ComplexConstant neg() {
	return new ComplexConstant( -re, -im );
    }
    /** Assigns <code>this</code> with negative sign.
     * Assigns <code>this</code> with negative sign: this = -this .
     */
    final void assignNeg() {
	re = -re;
	im = -im;
    }

    /** Assigns <code>this</code> with negative <code>u</code>.
     * <code>this = -u</code>.
     */
    final void assignNeg( final ComplexConstant u ) {
	re = -u.re;
	im = -u.im;
    }

    /** Assigns <code>this</code> with negative <code>u</code>.
     * <code>this = -u</code>.
     */
    final void assignNeg( final Field.Complex u ) {
        re = -u.getRe();
        im = -u.getIm();
    }

    /** Returns calculation root of <code>u</code>.
     * Creates a new <code>Field.Complex</code> that equals square root of <code>u</code>.
     * @return w = sqr(u)
     */
    public final static ComplexConstant sqr( final Field.Complex u ) {
      final double uRe=u.getRe();
      final double uIm=u.getIm();
      ComplexConstant w = new ComplexConstant();
      w.re = uRe * uRe - uIm * uIm;
      w.im = 2. * uRe * uIm;
      return w;
    }

    /** Returns calculation root of <code>this</code>.
     * Creates a new <code>Field.Complex</code> that equals square root of <code>this</code>
     * @return w = sqr(this)
     */
    public final ComplexConstant sqr() {
	return sqr( this );
    }
    /** Assigns <code>this</code> with root of <code>this</code>.
     * this = sqr(this).
     */
    final void assignSqr() {
	double rr = re;
	double ii = im;
	re = rr*rr - ii*ii;
	im = 2.*rr*ii;
    }

    /** Assigns <code>this</code> with root of <code>u</code>.
     * this = sqr(u).
     */
    final void assignSqr( ComplexConstant u) {
	if( this == u ) {
	    assignSqr();
	} else {
	    re = u.re*u.re - u.im*u.im;
	    im = 2.*u.re*u.im;
	}
    }

    /** Assigns <code>this</code> with root of <code>u</code>.
     * this = sqr(u).
     */
    final void assignSqr( Field.Complex u) {
      if( this == u ) {
          assignSqr();
      } else {
        final double uRe = u.getRe();
        final double uIm = u.getIm();
        re = uRe * uRe - uIm * uIm;
        im = 2. * uRe * uIm;
      }
    }

    /** Returns cube of <code>u</code>.
     * Creates new <code>Field.Complex</code> <code>w = u^3</code>.
     */
    public final static ComplexConstant cube( final Field.Complex u) {
        ComplexConstant w = new ComplexConstant();
        final double uRe = u.getRe();
        final double uIm = u.getIm();
        w.re = uRe * (uRe*uRe - 3.*uIm*uIm);
        w.im = uIm * (3.*uRe*uRe - uIm*uIm);
        return w;
    }

    /** Returns cube of <code>this</code>.
     * Creates new <code>Field.Complex</code> <code>w = this^3</code>.
     */
    public final ComplexConstant cube() {
	return cube( this );
    }
    /** Assigns <code>this</code> with cube of <code>this</code>
     *	<code>this = this^3</code>.
     */
    final void assignCube() {
	double rr = re;
	double ii = im;
	re = rr * (rr*rr - 3.*ii*ii);
	im = ii * (3.*rr*rr - ii*ii);
    }

    /** Assigns this cube of <code>u</code>
     * <code>this = u^3</code>.
     */
    final void assignCube( final ComplexConstant u ) {
	if( this == u ) {
	    assignCube();
	} else {
	    re = u.re * (u.re*u.re - 3.*u.im*u.im);
	    im = u.im * (3.*u.re*u.re - u.im*u.im);
	}
    }

    /** Assigns this cube of <code>u</code>
     * <code>this = u^3</code>.
     */
    final void assignCube( final Field.Complex u ) {
      if (this == u) {
        assignCube();
      }
      else {
        final double uRe = u.getRe();
        final double uIm = u.getIm();
        re = uRe * (uRe * uRe - 3. * uIm * uIm);
        im = uIm * (3. * uRe * uRe - uIm * uIm);
      }
    }

    /** Returns real part of <code>u</code>.
     * @return <code>u.re</code>.
     */
    public final static double re( final ComplexConstant u) {
	return u.re;
    }

    /** Returns imag. part of <code>u</code>.
     * @return <code>u.im</code>.
     */
    public final static double im(ComplexConstant u) {
	return u.im;
    }

    /** Returns exponent of <code>u</code>.
     * Creates new <code>Field.Complex</code> w = exp(u) */
    public final static ComplexConstant exp( final Field.Complex u) {
        ComplexConstant w = new ComplexConstant();
        double r = Math.exp(u.getRe());
        w.re = r * Math.cos(u.getIm());
        w.im = r * Math.sin(u.getIm());
        return w;
    }

    /** Returns exponent of <code>this</code>.
     * Creates new <code>Field.Complex</code> w = exp(this) */
    public final ComplexConstant exp() {
	return exp( this );
    }

    /** Assigns <code>this</code> with exponent of <code>this</code>.
     * <code>this = exp(this)</code>.
     */
    final void assignExp() {
	final double r = Math.exp(re);
	re = r * Math.cos( im );
	im = r * Math.sin( im );
    }

    /** Assigns <code>this</code> with exponent of <code>x+iy</code>.
     * <code>this = exp(x+iy)</code>.
     */
    final void assignExp( final double x, final double y ) {
	final double r = Math.exp(x);
	re = r * Math.cos(y);
	im = r * Math.sin(y);
    }

    /** Assigns <code>this</code> with exponent of <code>u</code>.
     * <code>this = exp(u)</code>.
     */
    final void assignExp( final ComplexConstant u ) {
	final double r = Math.exp(u.re);
	re = r * Math.cos(u.im);
	im = r * Math.sin(u.im);
    }

    /** Assigns <code>this</code> with exponent of <code>u</code>.
     * <code>this = exp(u)</code>.
     */
    final void assignExp( final Field.Complex u ) {
        final double r = Math.exp(u.getRe());
        re = r * Math.cos(u.getIm());
        im = r * Math.sin(u.getIm());
    }

    /** Returns product of log <code>u</code>.
     * Creates new <code>Field.Complex</code> that equals log(<code>u</code>).
     */
    public final static ComplexConstant log( final Field.Complex u) {
      ComplexConstant v=new ComplexConstant(u);
      v.assignLog();
      return v;
    }

    /** Returns product of log <code>this</code>.
     * Creates new <code>Field.Complex</code> that equals log(<code>this</code>).
     */
    public final ComplexConstant log() {
	return log( this );
     }

    /** Assigns <code>this</code> with log(<code>this</code>).
     * <code>this = log(this)</code>.
     */
    final void assignLog() {
	assignLog( this );
    }

    /** Assigns <code>this</code> with log(<code>u</code>).
     * <code>this = log(u)</code>.
     */
    final void assignLog( final ComplexConstant u ) {
	final double r = u.abs();

	if ( r == 0.) {
	    re = Double.NEGATIVE_INFINITY;
	    im = 0;
	} else {

	    final double arg = u.arg();

	    re = Math.log( r );
	    im = arg;
	}
    }

    /** Assigns <code>this</code> with log(<code>u</code>).
     * <code>this = log(u)</code>.
     */
    final void assignLog( final Field.Complex u ) {
      final double uRe=u.getRe();
      final double uIm=u.getIm();
      final double r = abs(uRe, uIm);

      if (r == 0.) {
        re = Double.NEGATIVE_INFINITY;
        im = 0;
      }
      else {

        final double arg = arg(uRe, uIm);

        re = Math.log(r);
        im = arg;
      }
    }

    /** Returns <code>w</code> = sqrt(<code>u</code>).
     * Creates new <code>Field.Complex</code> : w = sqrt(u).
     */
    public final static ComplexConstant sqrt( final Field.Complex u ) {
      final double uRe=u.getRe();
      final double uIm=u.getIm();
      double rr = Math.sqrt(abs(uRe, uIm));
      double ii = arg(uRe, uIm)/2.;
      return new ComplexConstant( rr*Math.cos(ii), rr*Math.sin(ii) );
    }

    /** Returns <code>w</code> = sqrt(<code>this</code>).
     * Creates new <code>Field.Complex</code> : w = sqrt(u).
     */
    public final ComplexConstant sqrt() {
	return sqrt( this );
    }

    /** Assigns <code>this</code> with sqrt(<code>this</code>).
     * <code>this</code> = sqrt(<code>this</code>).
     */
    final void assignSqrt() {
	double rr = Math.sqrt(abs());
	double ii = arg()/2.;
	re = rr*Math.cos(ii);
	im = rr*Math.sin(ii);
    }

    /** Assigns <code>this</code> with sqrt(<code>u</code>).
     * <code>this</code> = sqrt(<code>u</code>).
     */
    final void assignSqrt( final ComplexConstant u ) {
	final double rr = Math.sqrt(abs(u));
	final double ii = arg(u)/2.;
	re = rr*Math.cos(ii);
	im = rr*Math.sin(ii);
    }

    /** Assigns <code>this</code> with sqrt(<code>u</code>).
     * <code>this</code> = sqrt(<code>u</code>).
     */
    final void assignSqrt( final Field.Complex u ) {
      final double uRe=u.getRe();
      final double uIm=u.getIm();
      final double rr = Math.sqrt(abs(uRe, uIm));
      final double ii = arg(uRe, uIm)/2.;
      re = rr*Math.cos(ii);
      im = rr*Math.sin(ii);
    }

    /** Returns <code>u</code> to the pow of <code>v</code>.
     * Creates new <code>Field.Complex</code> w = u^v
     */
    public final static ComplexConstant pow( final Field.Complex u, final Field.Complex v ) {
         ComplexConstant w = new ComplexConstant(u);
         w.assignPow(v);
         return w;
    }

    /** Returns <code>this</code> to the pow of <code>v</code>.
     * Creates new <code>Field.Complex</code> w = this^v.
     */
    public final ComplexConstant pow(final Field.Complex v) {
      return pow(this, v);
    }

    /** Assigns <code>this</code> with product of <code>this^v</code>.*/
    final void assignPow(final ComplexConstant v) {
      assignPow(this, v);
    }

    /** Assigns <code>this</code> with product of <code>this^v</code>.*/
    final void assignPow(final Field.Complex v) {
      assignPow(this, v);
    }

    /** Assigns <code>this</code> with product of <code>(x=iy)^(u+iv)</code>.*/
    final void assignPow( final double x, final double y, final double u, final double v ) {
        //set to 1
        if (absSqr(u,v) == 0.) {
            re=1.;
            im=0.;
            return;
        }
       //setlog
        double r = abs( x, y );

         if ( r == 0. ) {
          re = 0;
          im = 0;
          return;
        }

        final double arg = arg(x,y);

        final double tmpRe = Math.log( r );
        final double tmpIm = arg;

        //setTimes
        double rr = tmpRe*u - tmpIm*v;
        double ii = tmpRe*v + tmpIm*u;
       //exp
        r = Math.exp(rr);
        re = r * Math.cos(ii);
        im = r * Math.sin(ii);
    }

    /** Assigns <code>this</code> with product of <code>u^v</code>.*/
       final void assignPow( final ComplexConstant u, final ComplexConstant v ) {
         assignPow( u.re, u.im, v.re, v.im );
       }

       /** Assigns <code>this</code> with product of <code>u^v</code>.*/
       final void assignPow(final Field.Complex u, final Field.Complex v) {
         assignPow(u.getRe(), u.getIm(), v.getRe(), v.getIm() );
       }

    /** Returns <code>u</code> to the pow of <code>r</code>.
     * Creates new <code>Field.Complex</code> w = u^r
     */
    public final static ComplexConstant pow( final Field.Complex u, final double r) {
      ComplexConstant c =new ComplexConstant();
      c.assignPow(u,r);
      return c;
    }

   /** Returns <code>this</code> to the pow of <code>r</code>.
     * Creates new <code>Field.Complex</code> w = this^r
     */
    public final ComplexConstant pow( final double r ) {
	return pow( this, r );
    }

    /** Assigns <code>this</code> with <code>this</code> to the pow of r.
     * <code>this = this^r</code>
     */
    final void assignPow( final double r ) {
         assignPow( this , r );
    }

    /** Assigns <code>this</code> with <code>u</code> to the pow of r.
     * <code>this = (x+iy)^r</code>
     */
    final void assignPow( final double x, final double y, final double r ) {

	final double rad = abs( x, y );
	double tmpRe,tmpIm;

	if ( rad == 0.) {

	    tmpRe = Double.NEGATIVE_INFINITY;
	    tmpIm = 0;

	} else {
	    final double arg = arg(x,y);

	    tmpRe = Math.log( rad );
	    tmpIm = arg;
	}

	tmpRe *= r;
	tmpIm *= r;

	final double rr = Math.exp(tmpRe);
	re = rr * Math.cos( tmpIm );
	im = rr * Math.sin( tmpIm );
    }

    /** Assigns <code>this</code> with <code>u</code> to the pow of r.
     * <code>this = u^r</code>
     */
    final void assignPow(final ComplexConstant u, final double r) {
      assignPow(u.re, u.im, r);
    }

    /** Assigns <code>this</code> with <code>u</code> to the pow of r.
     * <code>this = u^r</code>
     */
    final void assignPow(final Field.Complex u, final double r) {
      assignPow(u.getRe(), u.getIm(), r);
    }

    /** Assigns <code>this</code> with <code>r</code> to the pow of u.
     * <code>this = r^u</code>
     */
    final void assignPow(final double r, final ComplexConstant u) {
      assignPow(r, 0, u.re, u.im);
    }

    /** Assigns <code>this</code> with <code>r</code> to the pow of u.
     * <code>this = r^u</code>
     */
    final void assignPow(final double r, final Field.Complex u) {
      assignPow(r, 0, u.getRe(), u.getIm());
    }

    /** Assigns <code>this</code> with <code>u</code> to the pow of r.
     * <code>this = (x+iy)^r</code
     */
    final void assignPow(final double x, final double y, int r) {
      ComplexConstant u2n = new ComplexConstant(x,y); // u^{2^n}
      if (r < 0) {
        r = -r;
        u2n.assignInvert();
      }

      assign(1., 0.);
      while (r != 0) {
        if ( (r & 1) == 1)
          assignTimes(u2n);
        u2n.assignSqr(); // u^{2^{n+1}}
        r >>= 1; // r / 2
      }
    }

    /** Assigns <code>this</code> with <code>u</code> to the pow of r.
     * <code>this = u^r</code
     */
    final void assignPow(final ComplexConstant u, int r) {
      assignPow( u.re, u.im, r );
    }

    /** Assigns <code>this</code> with <code>u</code> to the pow of r.
     * <code>this = u^r</code
     */
    final void assignPow(final Field.Complex u, int r) {
      assignPow(u.getRe(), u.getIm(), r);
    }


    /** Returns <code>u</code> to the pow of <code>r</code>.
     * Creates new <code>Field.Complex</code> w = u^r
     */
    public final static ComplexConstant pow( final Field.Complex u, final int r) {
    ComplexConstant c = new ComplexConstant();
    c.assignPow(u,r);
    return c;
    }

    /** Returns <code>this</code> to the pow of <code>r</code>.
     * Creates new <code>Field.Complex</code> w = this^r
     */
    public final ComplexConstant pow( final int r ) {
      ComplexConstant c = new ComplexConstant();
      c.assignPow( re, im, r );
      return c;
    }

    /** Assigns <code>this</code> with <code>this</code> to the pow of r.
     * <code>this = this^r</code>
     */
    final void assignPow( final int r ) {
         assignPow( this, r );
    }

    /** Returns sinus of <code>u</code>.
     *	Creates a new <code>Field.Complex<code> w which equals sin(u)
     */
    public final static ComplexConstant sin( final Field.Complex u ) {
	ComplexConstant w = new ComplexConstant();
        w.assignSin( u.getRe(), u.getIm() );
	return w;
    }

    /** Returns sinus of <code>this</code>.
     *	Creates a new <code>Field.Complex<code> w which equals sin(this)
     */
    public final ComplexConstant sin() {
	return sin( this );
    }

    /** Assigns <code>this</code> with sinus <code>this</code>.
     *	this = sin(this)
     */
    final void assignSin() {
	assignSin( re, im );
    }

    /** Assigns <code>this</code> with sinus <code>x+iy</code>.
     *	this = sin(x+iy)
     */
    final void assignSin(final double x, final double y) {
      final double sinRe = Math.sin(x);
      final double cosRe = Math.cos(x);
      final double expIm = Math.exp(y);
      re = sinRe * (expIm + 1 / expIm) / 2;
      im = cosRe * (expIm - 1 / expIm) / 2;
    }

    /** Assigns <code>this</code> with sinus <code>u</code>.
     *	this = sin(u)
     */
    final void assignSin( final ComplexConstant u ) {
	assignSin( u.re, u.im );
    }

    /** Assigns <code>this</code> with sinus <code>u</code>.
     *	this = sin(u)
     */
    final void assignSin(final Field.Complex u) {
      assignSin( u.getRe(), u.getIm() );
    }

    /** Returns cosinus of <code>u</code>.
     *	Creates a new <code>Field.Complex<code> w which equals cos(x+iy)
     */
    public final static ComplexConstant cos( final Field.Complex u ) {
      ComplexConstant w = new ComplexConstant();
      w.assignCos( u.getRe(), u.getIm() );
      return w;
    }

    /** Returns cosinus of <code>this</code>.
     *	Creates a new <code>Field.Complex<code> w which equals cos(this)
     */
    public final ComplexConstant cos() {
	return cos( this );
    }

    /** Assigns <code>this</code> with cosinus <code>this</code>.
     *	this = cos(this)
     */
    final void assignCos() {
	assignCos( re, im );
    }

    /** Assigns <code>this</code> with cosinus <code>x+iy</code>.
     *	this = cos(x+iy)
     */
    final void assignCos(final double x, final double y) {
      final double cosRe = Math.cos(x);
      final double sinRe = Math.sin(x);
      final double expIm = Math.exp(y);
      re = cosRe * (expIm + 1 / expIm) / 2;
      im = -sinRe * (expIm - 1 / expIm) / 2;
    }

    /** Assigns <code>this</code> with cosinus <code>u</code>.
     *	this = cos(u)
     */
    final void assignCos(final ComplexConstant u) {
      assignCos( u.re, u.im );
    }

    /** Assigns <code>this</code> with cosinus <code>u</code>.
     *	this = cos(u)
     */
    final void assignCos(final Field.Complex u) {
      assignCos( u.getRe(), u.getIm() );
    }

    
    /** Returns Tangent of <code>u</code>.
     *	Creates a new <code>Field.Complex<code> w which equals Tan(x+iy)
     */
    public final static ComplexConstant tan( final Field.Complex u ) {
    	ComplexConstant w = new ComplexConstant();
    	w.assignTan( u.getRe(), u.getIm() );
    	return w;
    }

    /** Returns Tangent of <code>this</code>.
     *	Creates a new <code>Field.Complex<code> w which equals Tan(this)
     */
    public final ComplexConstant tan() {
    	return tan( this );
    }

    /** Assigns <code>this</code> with Tangent <code>this</code>.
     *	this = Tan(this)
     */
    final void assignTan() {
    	assignTan( re, im );
    }

    /** Assigns <code>this</code> with Tangent <code>x+iy</code>.
     *	this = Tan(x+iy)
     */
    final void assignTan(final double x, final double y) {
    	final double cosRe = Math.cos(x);
    	final double sinRe = Math.sin(x);
    	final double expIm = Math.exp(y);
    	final double reCos = cosRe * (expIm + 1 / expIm) / 2;
    	final double imCos = -sinRe * (expIm - 1 / expIm) / 2;
    	re = sinRe * (expIm + 1 / expIm) / 2;
    	im = cosRe * (expIm - 1 / expIm) / 2;
    	assignDivide( reCos, imCos );
    }

    /** Assigns <code>this</code> with Tangent <code>u</code>.
     *	this = Tan(u)
     */
    final void assignTan(final ComplexConstant u) {
    	assignTan( u.re, u.im );
    }

    /** Assigns <code>this</code> with Tangent <code>u</code>.
     *	this = Tan(u)
     */
    final void assignTan(final Field.Complex u) {
    	assignTan( u.getRe(), u.getIm() );
    }
    
    /** Returns sinushyperbolicus of <code>u</code>.
     *	Creates a new <code>Field.Complex<code> w which equals sinh(u)
     */
    public final static ComplexConstant sinh( final Field.Complex u ) {
	ComplexConstant w = new ComplexConstant();
        w.assignSinh( u.getRe(), u.getIm() );
	return w;
    }

    /** Returns sinushyperbolicus of <code>this</code>.
     *	Creates a new <code>Field.Complex<code> w which equals sinh(this)
     */
    public final ComplexConstant sinh() {
	return sinh( this );
    }
    /** Assigns <code>this</code> with sinushyperbolicus <code>this</code>.
     *	this = sinh(this)
     */
    final void assignSinh() {
	assignSinh( re, im );
    }

    /** Assigns <code>this</code> with sinushyperbolicus <code>x+iy</code>.
     *	this = sinh(x+iy)
     */
    final void assignSinh(final double x, final double y) {
      final double sinIm = Math.sin(y);
      final double cosIm = Math.cos(y);
      final double expRe = Math.exp(x);
      re = cosIm * (expRe - 1 / expRe) / 2;
      im = sinIm * (expRe + 1 / expRe) / 2;
    }

    /** Assigns <code>this</code> with sinushyperbolicus <code>u</code>.
     *	this = sinh(u)
     */
    final void assignSinh( final ComplexConstant u ) {
	assignSinh( u.re, u.im );
    }

    /** Assigns <code>this</code> with sinushyperbolicus <code>u</code>.
     *	this = sinh(u)
     */
    final void assignSinh(final Field.Complex u) {
      assignSinh(u.getRe(), u.getIm() );
    }

    /** Returns cosinushyperbolicus of <code>u</code>.
     *	Creates a new <code>Field.Complex<code> w which equals cosh(u)
     */
    public final static ComplexConstant cosh( final Field.Complex u ) {
	ComplexConstant w = new ComplexConstant();
	w.assignCosh( u.getRe(), u.getIm() );
	return w;
    }

    /** Returns cosinushyperbolicus of <code>this</code>.
     *	Creates a new <code>Field.Complex<code> w which equals cosh(this)
     */
    public final ComplexConstant cosh() {
	return cosh( this );
    }

    /** Assigns <code>this</code> with cosinushyperbolicus <code>this</code>.
     *	this = cosh(this)
     */
    final void assignCosh() {
      assignCosh( re, im );
    }

    /** Assigns <code>this</code> with cosinushyperbolicus <code>x+iy</code>.
     *	this = cosh(x+iy)
     */
    final void assignCosh( final double x, final double y ) {
	final double sinIm = Math.sin(y);
	final double cosIm = Math.cos(y);
	final double expRe = Math.exp(x);
	re = cosIm * ( expRe + 1/expRe ) / 2;
	im = sinIm * ( expRe - 1/expRe ) / 2;
    }

    /** Assigns <code>this</code> with cosinushyperbolicus <code>u</code>.
     *	this = cosh(u)
     */
    final void assignCosh(final ComplexConstant u) {
      assignCosh(u.re, u.im);
    }

    /** Assigns <code>this</code> with cosinushyperbolicus <code>u</code>.
     *	this = cosh(u)
     */
    final void assignCosh(final Field.Complex u) {
      assignCosh( u.getRe(), u.getIm() );
    }

    /** Returns <code>w = r*exp(i*phi)</code>.
     * Creates a new <code>Field.Complex</code> w = r*exp(i*phi)
     */
    public final static ComplexConstant fromPolar( final double r, final double f ) {
	ComplexConstant w = new ComplexConstant();
	w.re = r * Math.cos(f);
	w.im = r * Math.sin(f);
	return w;
    }

    /** Assigns <code>this = r*exp(i*phi)</code>. */
    final void assignFromPolar( final double r, final double f) {
	re = r * Math.cos(f);
	im = r * Math.sin(f);
    }

    /** Returns cross ratio: (c-b)*(a-d)/(b-a)/(d-c). */
   public final void assignCrossRatio(final double aRe, final double aIm,
                                      final double bRe, final double bIm,
                                      final double cRe, final double cIm,
                                      final double dRe, final double dIm ) {

           assignDivide( cRe-bRe, cIm-bIm , bRe-aRe, bIm-aIm );
           assignTimes(  aRe-dRe, aIm-dIm );
           assignDivide( dRe-cRe, dIm-cIm );
   }

   /** Returns cross ratio: (c-d)*(a-d)/(b-a)/(d-c). */
   public final static ComplexConstant crossRatio(final Field.Complex a,
                                                  final Field.Complex b,
                                                  final Field.Complex c,
                                                  final Field.Complex d) {
     ComplexConstant cr = new ComplexConstant();
     cr.assignCrossRatio(a, b, c, d);
     return cr;
           }

   /** Returns cross ratio: (c-d)*(a-d)/(b-a)/(d-c). */
   public final void assignCrossRatio( final ComplexConstant a,final ComplexConstant b,
                                       final ComplexConstant c,final ComplexConstant d ) {
     assignCrossRatio( a.re, a.im, b.re, b.im, c.re, c.im, d.re, d.im );
   }

   /** Returns cross ratio: (c-d)*(a-d)/(b-a)/(d-c). */
   public final void assignCrossRatio( final Field.Complex a,final Field.Complex b,
                                       final Field.Complex c,final Field.Complex d ) {
       assignCrossRatio( a.getRe(), a.getIm(), b.getRe(), b.getIm(),
                         c.getRe(), c.getIm(), d.getRe(), d.getIm() );
   }

    /** Assigns <code>this=random()</code>. */
    final void assignRandom() {
      re = Math.random();
      im = Math.random();
    }

    /** Returns <code>r = <u,v>_R^2</code>.
     * Creates double with sum of multiplication between real and imag. parts of <code>u,v</code>
     */
    public final static double dot(final ComplexConstant u,
                                   final ComplexConstant v) {
      return u.re * v.re + u.im * v.im;
    }

    /** Returns <code>r = <this,v>_R^2</code>.
     * Creates double with sum of multiplication between real and imag. parts of <code>this,v</code>
     */
    public final double dot( final ComplexConstant v ) {
	return re * v.re + im * v.im;
    }

    /** Returns <code>r = <this,v>_R^2</code>.
     * Creates double with sum of multiplication between real and imag. parts of <code>this,v</code>
     */
    public final double dot( final Field.Complex v ) {
      return re * v.getRe() + im * v.getIm();
    }

    /** Returns r = det(u,v).
     * Creates double with vlaue of determinant <code>u,v</code>.
     */
    public final static double det( final ComplexConstant u, final ComplexConstant v ) {
	return u.re * v.im - u.im * v.re;
    }

    /** Returns r = det(this,v).
     * Creates double with vlaue of determinant <code>this,v</code>.
     */
    public final double det( final ComplexConstant v ) {
	return re * v.im - im * v.re;
    }

    /** Returns r = det(this,v).
     * Creates double with vlaue of determinant <code>this,v</code>.
     */
    public final double det(final Field.Complex v) {
      return re * v.getIm() - im * v.getRe();
    }

    /** Returns r = |u-v|^2 */
    public final static double distSqr( final ComplexConstant u, final ComplexConstant v ) {
	double subRe = u.re - v.re;
	double subIm = u.im - v.im;

	return subRe*subRe + subIm*subIm;
    }

    /** Returns r = |this-v|^2 */
    public final double distSqr( final ComplexConstant v ) {
	double subRe = re - v.re;
	double subIm = im - v.im;

	return subRe*subRe + subIm*subIm;
    }

    /** Returns r = |this-v|^2 */
    public final double distSqr( final Field.Complex v ) {
      double subRe = re - v.getRe();
      double subIm = im - v.getIm();

      return subRe*subRe + subIm*subIm;
    }

    /** Returns r = |u-v| */
    public final static double dist( final ComplexConstant u, final ComplexConstant v ) {
	return Math.sqrt( distSqr( u, v ) );
    }

    /** Returns r = |this-v| */
    public final double dist( final ComplexConstant v ) {
	return Math.sqrt( this.distSqr( v ) );
    }

    /** Returns r = |this-v| */
    public final double dist(final Field.Complex v) {
      return Math.sqrt(this.distSqr(v));
    }

    /** Returns r = |u| */
    public final static double abs( final ComplexConstant u ) {
	return u.abs();
    }

    /** Retuurns r = |u| */
    public final double abs() {
	return Math.sqrt(re*re + im*im);
    }

    /** Returns |u.re+iu.im|^2 */
    public final static double absSqr( final ComplexConstant u ) {
	return u.re*u.re + u.im*u.im;
    }

    /** Returns r = |u|^2 */
    public final double absSqr() {
	return re*re + im*im;
    }

    /** Returns arg( u ) */
    public final static double arg( final ComplexConstant u ) {
	return arg( u.re, u.im );
    }

    /** Returns arg( this ) */
    public final double arg() {
	return arg( re, im );
    }

    /**Indicates whether <code>c</code> equal to <code>this</code>*/
    public final boolean equals( Field.Complex c, final double eps){
	return (Math.abs(re-c.getRe()) < eps) && (Math.abs(im-c.getIm()) < eps);
    }

    /**Indicates whether some other Field.Complex is "equal to" this one.*/
    public final boolean equals( final Field.Complex c )
    {
	    return this == c ? true : equals(c.getRe(), c.getIm());
    }

    /**Indicates whether some other Field.Complex is "equal to" this one.*/
    public final boolean equals( final ComplexConstant c )
    {
      return this == c ? true : equals(c.re, c.im);
    }


    /**Indicates whether x + iy is "equal to" this one.*/
    public final boolean equals(final double x, final double y, final double eps) {
      return (Math.abs(re - x) < eps) && (Math.abs(im - y) < eps);
    }

    /**Indicates whether x + iy is "equal to" this one.*/
    public final boolean equals( final double x, final double y )
    {
        return (Math.abs(re-x) < EPS) && (Math.abs(im-y) < EPS);
    }


  /**Indicates whether <code>c</code> equal to <code>this</code>*/
    public final boolean equals( ComplexConstant c, final double eps){
      return (Math.abs(re-c.re) < eps) && (Math.abs(im-c.im) < eps);
    }

  /**Indicates whether <code>c1</code> equal to <code>c2</code>*/
    public static final boolean arraysCoincide( final ComplexConstant[] c1, final ComplexConstant[] c2, final double eps ){
	if( c1==null && c2==null ) return true;
	if(( c1== null && c2!=null )|| ( c2== null && c1!=null )||( c2.length!=c1.length))
	    return false;
	final int l = c1.length;
	for(int i=0; i<l;i++){
	    if( !c1[i].equals(c2[i], eps)) return false;
	}
	return true;
    }


  /**Indicates whether <code>c1</code> equal to <code>c2</code>*/
    public static final boolean arraysCoincide( final ComplexConstant[] c1, final ComplexConstant[] c2){
      return ComplexConstant.arraysCoincide(c1, c2, ComplexConstant.EPS );
    }

     /**Indicates whether <code>c1</code> equal to <code>c2</code> as sets.*/
        public static final boolean setsCoincide( final ComplexConstant[] c1, final ComplexConstant[] c2 ){
          return ComplexConstant.setsCoincide( c1, c2, ComplexConstant.EPS );
        }

 /**Indicates whether <code>c1</code> equal to <code>c2</code> as sets.*/
    public static final boolean setsCoincide( final ComplexConstant[] c1, final ComplexConstant[] c2,
					   final double eps ){
   	if( c1==null && c2==null ) return true;
	if(( c1== null && c2!=null )|| ( c2== null && c1!=null )||( c2.length!=c1.length))
	    return false;
	final int l = c1.length;
	LinkedList liste = new LinkedList();
	for( int i = 0; i< l; i++)
	    liste.add( c2[i] );
	ListIterator it;
	for( int i=0;i<l;i++){
	    it = liste.listIterator( 0 );
	    while( it.hasNext()){
		if( c1[i].equals( (ComplexConstant)it.next(),eps )) {
		    it.remove();// remove aus der "vergleichliste "
		    break;
		}else if( !it.hasNext()) return false;
	    }
	}
	return liste.size()==0;
 }



    /** Copies the values (not references) of source array <code>s</code>
	starting from position <code>is</code> into the destination array
	<code>d</code> starting form position <code>ds</code>.
	The methods creates instances of <code>ComplexConstant</code>
	and handles further more the case that source and desination
	overlap.
    */
    public static final void arraycopy( final Field.Complex [] s, int is,
					final ComplexConstant [] d, int id, final int size ) {

	if( !( is+size <= s.length && id+size <= d.length ) )
	    throw new IndexOutOfBoundsException();

	if( s == d && id == is )
	    return;

	if( id <= is ) {

	    for( int i=0; i<size; i++, is++, id++ ) {

		if( s[is] == null )
		    d[id] = null;
		else
		    d[id] = new ComplexConstant( s[is] );
	    }

	} else {

	    is += size-1;
	    id += size-1;

	    for( int i=0; i<size; i++, is--, id-- ) {

		if( s[is] == null )
		    d[id] = null;
		else
		    d[id] = new ComplexConstant( s[is] );
	    }
	}
    }

    /** Returns copy of array <code>s</code>.
	The methods copies the values of <code>s</code>
	not only the references, therefore it differs
	from <code>clone()</code>.
    */
    public static final ComplexConstant[] copy( final Field.Complex [] s ){

	final ComplexConstant [] d = new ComplexConstant[s.length];

	arraycopy( s, 0, d, 0, s.length );

	return d;
    }

}