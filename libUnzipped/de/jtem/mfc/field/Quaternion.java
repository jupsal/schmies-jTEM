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

import de.jtem.mfc.vector.Real3;
/**
 * This class provides a comprehensive library for quaternions.
 * The quaternions are represented thru four real values:
 * q = re = ix + iy = iz.
 * <p>
 * The four fields <code>re,x,y,z</code> can be accessed publicly.
 * <p>
  * We impose the following nameing conventions for methods:
  * <p>
  * <em>Method names</em> refer to the action performed:
  * <br>
  *    For instances <code>a,b</code> of class <code>Quaternion</code>, <code>a.plus(b)</code>
  *    produces a new instance of "Quaternion" containing the sum <code>a+b</code>.
  * <p>
  * <em>Assign-methods</em> do NOT create new instances:
  * <br>
  *    For instances <code>a,b,c</code> of class <code>Quaternion</code>,
  *    <code>c.assignPlus(a,b)</code> stores the sum <code>a+b</code> in <code>c</code>.
  *
  *    <code>a.assignPlus(b)</code> stores the sum <code>a+b</code> in a,
  *    like the operator "+=" does for simple types.
  * <p>
  * All non-commutatitve operations are RIGHT OPERATIONS !
  *<br>
  *   For instances <code>a,b</code> of class <code>Quaternion</code>, <code>a.assignTimes(b)</code> stores
  *   the product <code>a*b</code> in <code>a</code>. To store <code>b*a</code> in <code>a</code>,
  *   you would have to call <code>a.assignTimes(b,a)</code>.
  * <p>
  *
  * Quaternions INTERACT WITH DOUBLES:
  * <br>
  *   For an instance <code>a</code> of class <b>Quaternion</b> and a double <code>r</code>,
  *   <code>a.times(r)</code> creates a new Quaternion containing the product <code>a*r</code>.
  *
  *
  * @author	marina
  */
 public class Quaternion implements Serializable, Cloneable {

   final static double EPS    = 1.e-14;
   final static double EPSSQR = EPS * EPS;

    private static final long serialVersionUID = 1L;
     /** @deprecated */
    public static final Quaternion ZERO       =  new Quaternion( 0.0, 0.0, 0.0, 0.0);
     /** @deprecated */
    public static final Quaternion ONE        =  new Quaternion( 1.0, 0.0, 0.0, 0.0);
     /** @deprecated */
    public static final Quaternion I          =  new Quaternion( 0.0, 1.0, 0.0, 0.0);
     /** @deprecated */
    public static final Quaternion J          =  new Quaternion( 0.0, 0.0, 1.0, 0.0);
     /** @deprecated */
    public static final Quaternion K          =  new Quaternion( 0.0, 0.0, 0.0, 1.0);
     /** @deprecated */
     public static final Quaternion NEG_ONE   =  new Quaternion(-1.0, 0.0, 0.0, 0.0);
	/** @deprecated */
    public static final Quaternion NEG_I      =  new Quaternion( 0.0,-1.0, 0.0, 0.0);
     /** @deprecated */
    public static final Quaternion NEG_J      =  new Quaternion( 0.0, 0.0,-1.0, 0.0);
     /** @deprecated */
    public static final Quaternion NEG_K      =  new Quaternion( 0.0, 0.0, 0.0,-1.0);
     /** real part of quaternion.*/
    public double re;
     /** component of imaginary part of quaternion.*/
    public double x, y, z;

    /** Creates zero quaternion. */
    public Quaternion () {
    }

    /** Creates purely real quaternion. */
    public Quaternion ( final double re ) {
      this.re = re;
    }

    /**
     * Creates purely imaginary quaternion.
     * @param im imaginary part of quaternion.
     */
    public Quaternion ( final Real3 im ) {
      x = im.x;
      y = im.y;
      z = im.z;
    }

    /** Creates quaternion u + j*v. */
    public Quaternion ( final Complex u, final Complex v) {
      re =  u.re;
      x  =  u.im;
      y  =  v.re;
      z  = -v.im;
    }

    /** Creates quaternion re + ix + jy + kz. */
    public Quaternion ( final double re, final double x,  final double y,  final double z ) {
      this.re = re;
      this.x  = x;
      this.y  = y;
      this.z  = z;
    }

    /** Creates quaternion re + ix + jy + kz with im = (x,y,z). */
    public Quaternion ( final double re, final Real3 im ) {
     this.re = re;
     this.x  = im.x;
     this.y  = im.y;
     this.z  = im.z;
   }

    /** Creates <copy of <code>q</code>.*/
    public Quaternion ( final Quaternion q ) {
      assign(q);
    }

    /** Returns real part. */
    public final double getRe() {
      return re;
    }

    /** Returns imaginary part.*/
    public final Real3 getIm() {
      return new Real3( x,y,z);
    }

    /** Returns x component of imaginary part. */
    public final double getX() {
      return x;
    }

    /** Return y component of imaginary part.*/
    public final double getY() {
      return y;
    }

    /** Return z component of imaginary part.*/
    public final double getZ() {
      return z;
    }

    /** Sets real part.
     * @param re real part
     */
    public final void setRe( final double re) {
      this.re = re;
    }

    /** Sets imaginary part.
     * @param im imaginary part
     */
    public final void setIm( final Real3 im ) {
      this.x = im.x;
      this.y = im.y;
      this.z = im.z;
    }

    /** Sets x component of imaginary part.
     * @param double x
     */
    public final void setX( final double x) {
      this.x = x;
    }

    /** Sets y component of imaginary part.
     * @param double y
     */
    public final void setY( final double y ) {
      this.y = y;
    }

     /** Sets z component of imaginary part.
      * @param double z
      */
    public final void setZ( final double z) {
      this.z = z;
    }

    /** @deprecated
     * @see #assign(double)
     * @return <code>this</code>
     */
    public final Quaternion set( final double aR) {
      re = aR; x = 0.0; y = 0.0; z = 0.0;
      return this;
    }

    /**
     * Assigns this with re.
     */
    public final void assign( final double re ) {
      this.re = re;
      x = 0.0;
      y = 0.0;
      z = 0.0;
    }

    /**  @deprecated quaternion = complex1 + j*complex2
     * @see #assign(Field.Complex,Field.Complex)
     * @return <code>this</code>
     */
    public final Quaternion set( final Complex u, final Complex v) {
      re = u.re; x = u.im; y = v.re; z = -v.im;
      return this;
    }

    /**
     * Assigns this with ix + iy + iz with im=(x,y,z).
     */
    public final void assign( final Real3 im ) {
      this.re = 0;
      x = im.x;
      y = im.y;
      z = im.z;
    }

    /**
     * Assigns this with <code>u + i v</code>.
     */
    public final void assign( final Complex u, final Complex v) {
      re =  u.re;
      x  =  u.im;
      y  =  v.re;
      z  = -v.im;
    }

    /** @deprecated
     * @see #assign( double aR, double aI, double aJ, double aK )
     * @return <code>this</code>
     */
    public final Quaternion set( final double aR,  final double aI,  final double aJ,  final double aK) {
      re = aR; x = aI; y = aJ; z = aK;
      return this;
    }

    /**
     * Assigns this with re + ix + iy + iz.
     */
    public final void assign( final double re,  final double x,  final double y,  final double z ) {
      this.re = re;
      this.x  = x;
      this.y  = y;
      this.z  = z;

    }

    /** @deprecated
     * @see #assign( Quaternion q )
     * @return <code>this</code>
     */
    public final Quaternion set( final Quaternion q) {
      re = q.re; x = q.x; y = q.y; z = q.z;
      return this;
    }

    /**
     * Assigns <code>this</code> with  q.
     */
    public final Quaternion assign( final Quaternion q) {
      re = q.re; x = q.x; y = q.y; z = q.z;
      return this;
    }

    /**
     * Creates copy of this.
     * @return copy of this
     */
    public final Quaternion copy() {
      return new Quaternion(re,x,y,z);
    }

    /** Clones this one
     * @return copy of <code>this</code> Object
     * @see #copy()
     */
    public final Object clone() {
      return copy();
    }

    /**
     * Checks whether this has any infinite entries.
     * @return <tt>true</tt> if infinite, <tt>false</tt> otherwise
     */
    public final boolean isInfinite () {
      return ((Double.isInfinite(re) || Double.isInfinite(x)) || (Double.isInfinite(y) || Double.isInfinite(z)));
    }

   /**
    * Checks whether this has any Not-a-Number(NaN) entries.
    * @return <tt>true</tt> if is NaN, <tt>false</tt> otherwise
    */
    public final boolean isNaN () {
      return ((Double.isNaN(re) || Double.isNaN(x)) || (Double.isNaN(y) || Double.isNaN(z)));
    }


    /**
     * Creates sum of this and <code>q</code>.
     */
    public final Quaternion plus( final Quaternion q) {
      return new Quaternion(re+q.re, x+q.x, y+q.y, z+q.z);
    }

    /**
     * Creates sum of this and purley imaginary quaternion given by <code>im</code>.
     */
    public final Quaternion plus( final Real3 im) {
      return new Quaternion(re, x+im.x, y+im.y, z+im.z);
    }

    /**
     * Creates sum of this and real number <code>re</code>.
     */
    public final Quaternion plus( final double re) {
      return new Quaternion(this.re+re, x, y, z);
    }

    /**  @deprecated this += q
     * @see #assignPlus( Quaternion q)
     * @return this
     */
    public final Quaternion setPlus( final Quaternion q) {
      re += q.re; x += q.x; y += q.y; z += q.z;
      return this;
    }

    /**
     * Adds <code>q</code> to this.
     */
    public final void assignPlus( final Quaternion q) {
      re += q.re;
      x += q.x;
      y += q.y;
      z += q.z;
    }

    /**
     * Adds purely imagrinary quaternion given by <code>im</code> to this.
     */
    public final void assignPlus( final Real3 im ) {
      x += im.x;
      y += im.y;
      z += im.z;
    }

    /**
     * Adds real number re to this.
     */
    public final void assignPlus( final double re) {
      this.re += re;
    }

    /**  @deprecated this = q + qq
     * @see #assignPlus( Quaternion q, Quaternion qq )
     * @return this
     */
    public final Quaternion setPlus( final Quaternion q, final Quaternion qq) {
      re = q.re + qq.re; x = q.x + qq.x; y = q.y + qq.y; z = q.z + qq.z;
      return this;
    }

    /**
     * Assigns <code>this</code> with product of sum between <code>q</code> and <code>qq</code>
     * @param <code>q</code> and <code>qq</code> as summands
     */

    public final void assignPlus( final Quaternion q, final Quaternion qq) {
      re = q.re + qq.re; x = q.x + qq.x; y = q.y + qq.y; z = q.z + qq.z;
    }

    /** @deprecated this = q + aR
     * @see #assignPlus( Quaternion q, double aR )
     * @return this
     */

    public final Quaternion setPlus( final Quaternion q, final double aR) {
      re = q.re + aR; x = q.x; y = q.y; z = q.z;
      return this;
    }

    /** Assigns <code>this</code> with product of sum between <code>q</code> and <code>aR</code>
     * @param <code>q</code> and <code>aR</code> as summands
     */
    public final void assignPlus( final Quaternion q, final double aR) {
      re = q.re + aR; x = q.x; y = q.y; z = q.z;
    }

    /** Creates a new <code>Quaternion</code> C which is product of substraction <code>q</code> from <code>this</code>
     * @param <code>q</code> as subtrahend
     * @return C = this - q
     */
    public final Quaternion minus( final Quaternion q) {
      return new Quaternion(re-q.re, x-q.x, y-q.y, z-q.z);
    }
    /** Creates a new <code>Quaternion</code> C which is product of substraction <code>aR</code> from <code>this</code>
     * @param <code>aR</code> as subtrahend
     * @return C = this - aR
     */
    public final Quaternion minus( final double aR) {
      return new Quaternion(re-aR, x, y, z);
    }
    /** @deprecated  this -= q
     * @see #assignMinus( Quaternion q )
     * @return this
     */
    public final Quaternion setMinus( final Quaternion q) {
      re -= q.re; x -= q.x; y -= q.y; z -= q.z;
      return this;
    }
     /** Assigns <code>this<code> with a product  of substraction <code>q</code> from <code>this</code>
     * @param <code>q</code> as subtrahend
     */
    public final void assignMinus( final Quaternion q) {
      re -= q.re; x -= q.x; y -= q.y; z -= q.z;
    }
    /**  @deprecated this -= aR
     * @see #assignMinus( double aR )
     * @return this
     */
    public final Quaternion setMinus( final double aR) {
      re -= aR;
      return this;
    }
    /** Assigns <code>this<code> with a product  of substraction <code>aR</code> from <code>this</code>
     * @param <code>aR</code> as subtrahend
     */
    public final void assignMinus( final double aR) {
      re -= aR;
    }
    /** @deprecated this = q - qq
     * @see #assignMinus( Quaternion q, Quaternion qq )
     * @return this
     */
    public final Quaternion setMinus( final Quaternion q, final Quaternion qq) {
      re = q.re - qq.re; x = q.x - qq.x; y = q.y - qq.y; z = q.z - qq.z;
      return this;
    }
    /** Assigns <code>this<code> with a product  of substraction <code>qq</code> from <code>q</code>
     * @param <code>qq</code> as subtrahend and <code>q</code>
     */
    public final void assignMinus( final Quaternion q, final Quaternion qq) {
      re = q.re - qq.re; x = q.x - qq.x; y = q.y - qq.y; z = q.z - qq.z;
    }
    /** @deprecated this = q - aR
     * @see #assignMinus( Quaternion q, double aR )
     * @return this
     */
    public final Quaternion setMinus( final Quaternion q, final double aR) {
      re = q.re - aR; x = q.x; y = q.y; z = q.z;
      return this;
    }
    /**  Assigns <code>this<code> with a product  of substraction <code>aR</code> from <code>q</code>
     * @param <code>aR</code> as subtrahend and <code>q</code>
     */
    public final void assignMinus( final Quaternion q, final double aR) {
      re = q.re - aR; x = q.x; y = q.y; z = q.z;
    }

    /**
     * Creates product of this and <code>q</code>.
     */
    public final Quaternion times( final Quaternion q) {
	return new Quaternion(this.re*q.re-this.x*q.x-this.y*q.y-this.z*q.z,
			      this.re*q.x+this.x*q.re+this.y*q.z-this.z*q.y,
			      this.re*q.y-this.x*q.z+this.y*q.re+this.z*q.x,
			      this.re*q.z+this.x*q.y-this.y*q.x+this.z*q.re);
    }

    /**
    * Creates product of this and purely imaginary quaternion given by <code>im</code>.
    */
   public final Quaternion times( final Real3 im) {
       return new Quaternion(-this.x*im.x-this.y*im.y-this.z*im.z,
                             this.re*im.x+this.y*im.z-this.z*im.y,
                             this.re*im.y-this.x*im.z+this.z*im.x,
                             this.re*im.z+this.x*im.y-this.y*im.x);
   }


   /**
    * Creates product of this and real number.
    */
   public final Quaternion times( final double re) {
     return new Quaternion(this.re*re, x*re, y*re, z*re);
   }

   /** @deprecated this *= q
    * @see #assignTimes( Quaternion q)
    */
   public final Quaternion setTimes( final Quaternion q) {
      double rr=re,ii=x,jj=y,kk=z;
      if (q==this) {
	re = rr*rr-ii*ii-jj*jj-kk*kk; x = rr*ii+ii*rr+jj*kk-kk*jj;
	y = rr*jj-ii*kk+jj*rr+kk*ii; z = rr*kk+ii*jj-jj*ii+kk*rr;
      } else {
	re = rr*q.re-ii*q.x-jj*q.y-kk*q.z; x = rr*q.x+ii*q.re+jj*q.z-kk*q.y;
	y = rr*q.y-ii*q.z+jj*q.re+kk*q.x; z = rr*q.z+ii*q.y-jj*q.x+kk*q.re;
      }
      return this;
    }

    /**
     * Multiplies <code>q</code> from the right to this.
     * this*=q
     */
    public final void assignTimes( final Quaternion q) {
      double rr=re,ii=x,jj=y,kk=z;
      if (q==this) {
	re = rr*rr-ii*ii-jj*jj-kk*kk; x = rr*ii+ii*rr+jj*kk-kk*jj;
	y = rr*jj-ii*kk+jj*rr+kk*ii; z = rr*kk+ii*jj-jj*ii+kk*rr;
      } else {
	re = rr*q.re-ii*q.x-jj*q.y-kk*q.z; x = rr*q.x+ii*q.re+jj*q.z-kk*q.y;
	y = rr*q.y-ii*q.z+jj*q.re+kk*q.x; z = rr*q.z+ii*q.y-jj*q.x+kk*q.re;
      }
    }

    /**
     * Multiplies purely imaginar quaternion given by <code>im</code> from the right to this.
     * this*=im
     */
    public final void assignTimes(final Real3 im) {
      double rr = re, ii = x, jj = y, kk = z;
      re = - ii * im.x - jj * im.y - kk * im.z;
      x = rr * im.x + jj * im.z - kk * im.y;
      y = rr * im.y - ii * im.z + kk * im.x;
      z = rr * im.z + ii * im.y - jj * im.x;
    }



     /** @deprecated this *= aR
     * @see #assignTimes( double aR )
     * @return this
     */
    public final Quaternion setTimes( final double aR) {
      re *= aR; x *= aR; y *= aR; z *= aR;
      return this;
    }
     /**
      * Multiplies this with real number.
      */
    public final void assignTimes( final double aR) {
      re *= aR; x *= aR; y *= aR; z *= aR;
    }

     /** @deprecated this = q * qq
     * @see #assignTimes( Quaternion q, Quaternion qq )
     * @return this
     */
    public final Quaternion setTimes( final Quaternion q,  final Quaternion qq) {
      if (q==this || qq==this) {
	double rr = q.re*qq.re-q.x*qq.x-q.y*qq.y-q.z*qq.z;
	double ii = q.re*qq.x+q.x*qq.re+q.y*qq.z-q.z*qq.y;
	double jj = q.re*qq.y-q.x*qq.z+q.y*qq.re+q.z*qq.x;
	z = q.re*qq.z+q.x*qq.y-q.y*qq.x+q.z*qq.re; re=rr; x=ii; y=jj;
      } else {
	re = q.re*qq.re-q.x*qq.x-q.y*qq.y-q.z*qq.z;
	x = q.re*qq.x+q.x*qq.re+q.y*qq.z-q.z*qq.y;
	y = q.re*qq.y-q.x*qq.z+q.y*qq.re+q.z*qq.x;
	z = q.re*qq.z+q.x*qq.y-q.y*qq.x+q.z*qq.re;
      }

      return this;
    }

    /**
     * Assigns this with product of <code>a</code> and <code>b</code>.
     */
    public final void assignTimes( final Quaternion a, final Quaternion b) {
      if (a==this || b==this) {
	double rr = a.re*b.re-a.x*b.x-a.y*b.y-a.z*b.z;
	double ii = a.re*b.x+a.x*b.re+a.y*b.z-a.z*b.y;
	double jj = a.re*b.y-a.x*b.z+a.y*b.re+a.z*b.x;
	z = a.re*b.z+a.x*b.y-a.y*b.x+a.z*b.re; re=rr; x=ii; y=jj;
      } else {
	re = a.re*b.re-a.x*b.x-a.y*b.y-a.z*b.z;
	x = a.re*b.x+a.x*b.re+a.y*b.z-a.z*b.y;
	y = a.re*b.y-a.x*b.z+a.y*b.re+a.z*b.x;
	z = a.re*b.z+a.x*b.y-a.y*b.x+a.z*b.re;
      }
    }


    /**
     * Assigns this with product of <code>a</code> and <code>b</code>,
     * where <code>b</code> is purley imagiary.
     */
    public final void assignTimes(final Quaternion a, final Real3 b) {
      if (a == this ) {
        double rr = - a.x * b.x - a.y * b.y - a.z * b.z;
        double ii = a.re * b.x + a.y * b.z - a.z * b.y;
        double jj = a.re * b.y + a.z * b.x;
        z = a.re * b.z + a.x * b.y - a.y * b.x;
        re = rr;
        x = ii;
        y = jj;
      }
      else {
        re = - a.x * b.x - a.y * b.y - a.z * b.z;
        x = a.re * b.x + a.y * b.z - a.z * b.y;
        y = a.re * b.y - a.x * b.z + a.z * b.x;
        z = a.re * b.z + a.x * b.y - a.y * b.x;
      }
    }


    /**
     * Assigns this with product of <code>a</code> and <code>b</code>,
     * where <code>a</code> is purley imagiary.
     */
    public final void assignTimes(final Real3 a, final Quaternion b) {
      if (b == this) {
        double rr = -a.x * b.x - a.y * b.y - a.z * b.z;
        double ii = a.x * b.re + a.y * b.z - a.z * b.y;
        double jj = -a.x * b.z + a.y * b.re + a.z * b.x;
        z = a.x * b.y - a.y * b.x + a.z * b.re;
        re = rr;
        x = ii;
        y = jj;
      }
      else {
        re = -a.x * b.x - a.y * b.y - a.z * b.z;
        x = a.x * b.re + a.y * b.z - a.z * b.y;
        y = -a.x * b.z + a.y * b.re + a.z * b.x;
        z = a.x * b.y - a.y * b.x + a.z * b.re;
      }
    }

     /** @deprecated this = q * aR
     * @see #assignTimes( Quaternion q, double aR )
     * @return this
     */
    public final Quaternion setTimes( final Quaternion q,  final double aR) {
      re = q.re*aR; x = q.x*aR; y = q.y*aR; z = q.z*aR;
      return this;
    }

     /**
      * Assigns this with product of <code>q</code> and real number <code>re</code>,.
      */
    public final void assignTimes( final Quaternion q, final double re) {
      this.re = q.re*re; x = q.x*re; y = q.y*re; z = q.z*re;
    }

    /** Returns the lengthSqr of <code>this</code>
     * @return lengthSqr as double
     */
    public final double lengthSqr() {
      return re*re + x*x + y*y + z*z;
    }
    /** Returns the lengthSqr of <code>q</code>
     * @return lengthSqr as double
     */
    public final static double lengthSqr( final Quaternion q) {
      return q.re*q.re + q.x*q.x + q.y*q.y + q.z*q.z;
    }
    /** Returns the length of <code>this</code>
     * @return lengthSqr as double
     */
    public final double length() {
      return Math.sqrt(lengthSqr());
    }
    /** Returns the length of <code>q</code>
     * @return lengthSqr as double
     */
    public final static double length( final Quaternion q) {
      return Math.sqrt(lengthSqr(q));
    }
    /** ABS = LENGTH */
    public final double absSqr() {
      return re*re + x*x + y*y + z*z;
    }
    /** ABS = LENGTH */
    public final static double absSqr( final Quaternion q) {
      return q.re*q.re + q.x*q.x + q.y*q.y + q.z*q.z;
    }
        /** ABS = LENGTH */
    public final double abs() {
      return Math.sqrt(absSqr());
    }
        /** ABS = LENGTH */
    public final static double abs( final Quaternion q) {
      return Math.sqrt(absSqr(q));
    }


    /** In case of |q|==0 we set q^-1 = (Double.POSITIVE_INFINITY,0,0,0) */
    /** creates a new <code>Quaternion</code> C = this^-1 */
    public final Quaternion invert() {
      double l=lengthSqr();
      if (l==0.0) return new Quaternion(Double.POSITIVE_INFINITY,0.0,0.0,0.0);
      return new Quaternion(re/l,-x/l,-y/l,-z/l);
    }
    /** Creates a new <code>Quaternion</code> that is invert of <code>q</code>.
     * C = q^-1
     */
    public final static Quaternion invert( final Quaternion q) {
      double l=lengthSqr(q);
      if (l==0.0) return new Quaternion(Double.POSITIVE_INFINITY,0.0,0.0,0.0);
      return new Quaternion(q.re/l,-q.x/l,-q.y/l,-q.z/l);
    }
    /** @deprecated this = this^-1
     * @see #assignInvert()
     * @return this
     */
    public final Quaternion setInvert() {
      double l=lengthSqr();
      if (l==0.0) { re = Double.POSITIVE_INFINITY; x = 0.0; y = 0.0; z = 0.0; }
      else { re = re/l; x = -x/l; y = -y/l; z = -z/l; }
      return this;
    }
    /** Assigns <code>this</code> with <ode>this</code>^-1.
     * this = this^-1
     */
    public final void assignInvert() {
      double l=lengthSqr();
      if (l==0.0) { re = Double.POSITIVE_INFINITY; x = 0.0; y = 0.0; z = 0.0; }
      else { re = re/l; x = -x/l; y = -y/l; z = -z/l; }
    }
     /** @deprecated this = q^-1
     * @see #assignInvert( Quaternion q)
     * @return this
     */
    public final Quaternion setInvert( final Quaternion q) {
      double l=lengthSqr(q);
      if (l==0.0) { re = Double.POSITIVE_INFINITY; x = 0.0; y = 0.0; z = 0.0; }
      else { re = q.re/l; x = -q.x/l; y = -q.y/l; z = -q.z/l; }
      return this;
    }
     /**  Assigns <code>this</code> with <ode>q</code>^-1.
      * this = q^-1
      * @param <code>q</code>
      */
    public final void assignInvert( final Quaternion q) {
      double l=lengthSqr(q);
      if (l==0.0) { re = Double.POSITIVE_INFINITY; x = 0.0; y = 0.0; z = 0.0; }
      else { re = q.re/l; x = -q.x/l; y = -q.y/l; z = -q.z/l; }
    }

    /** Creates a new <code>Quaternion</code> C which is productof division between  <code>this</code> and <code>q</code>
     * @param <code>q</code> as divisor
     * @return C = this * q^-1
     */
    public final Quaternion divide( final Quaternion q) {
      double l=lengthSqr(q);

      if (l == 0.0) return new Quaternion( re/l , x/l , y/l, z/l);

      double rTmp,iTmp,jTmp,kTmp; // temp variables

      if (l==0.0) { rTmp = Double.POSITIVE_INFINITY; iTmp = 0.0; jTmp = 0.0; kTmp = 0.0; }
      else { rTmp = q.re/l; iTmp = -q.x/l; jTmp = -q.y/l; kTmp = -q.z/l; }

	return new Quaternion(this.re*rTmp-this.x*iTmp-this.y*jTmp-this.z*kTmp,
			      this.re*iTmp+this.x*rTmp+this.y*kTmp-this.z*jTmp,
			      this.re*jTmp-this.x*kTmp+this.y*rTmp+this.z*iTmp,
			      this.re*kTmp+this.x*jTmp-this.y*iTmp+this.z*rTmp);


    }

    /** Creates a new <code>Quaternion</code> C which is productof division between  <code>this</code> and <code>q</code>
     * @param <code>aR</code> as divisor
     * @return C = q/aR
     */
    public final Quaternion divide( final double aR) {
      return new Quaternion(re/aR, x/aR, y/aR, z/aR);
    }

    /** @deprecated this *= q^-1
     * @see #assignDivide( Quaternion q)
     * @return this
     */
    public final Quaternion setDivide( final Quaternion q) {
	assignDivide(q);
	return this;
    }

    /**  Assigns <code>this</code> with a product of division between <code>this</code> and <code>q</code>: this *= q^-1
     * @param <code>q</code> as divisor
     */
    public final void assignDivide( final Quaternion q) {
      double rTmp,iTmp,jTmp,kTmp; // temp variables
      double l=lengthSqr(q);
      if (l == 0.0)
      {
	  re /=l; y /= l;x /=l; z /=l; //       r /=0.0;i /= 0.0;j /=0.0;k /=0.0;

       rTmp = Double.POSITIVE_INFINITY; iTmp = 0.0; jTmp = 0.0; kTmp = 0.0;
       }
      else { rTmp = q.re/l; iTmp = -q.x/l; jTmp = -q.y/l; kTmp = -q.z/l; }
      assignTimes(new Quaternion(rTmp,iTmp,jTmp,kTmp));
    }

    /** @deprecated this /= aR
     * @see #assignDivide( double aR )
     * @return this
     */
    public final Quaternion setDivide( final double aR) {
      re /= aR; x /= aR; y /= aR; z /= aR;
      return this;
    }
    /**  Assigns <code>this</code> with a product of division between <code>this</code> and <code>aR</code>.
     * this /=aR
     * @param aR as divisor
     */
    public final void assignDivide( final double aR) {
      re /= aR; x /= aR; y /= aR; z /= aR;
    }

    /** @deprecated this = q * qq^-1
     * @see #assignDivide( Quaternion q, Quaternion qq )
     * @return this
    */
    public final Quaternion setDivide( final Quaternion q, final Quaternion qq) {
	assignDivide(q,qq);
       return this;
    }

    /**  Assigns <code>this</code> with a product of division between <code>q</code> and <code>qq</code>
     * @param <code>q</code> as dividend
     * @param <code>qq</code> as divisor
     */
    public final void assignDivide( final Quaternion q, final Quaternion qq) {
      double l=lengthSqr(qq);
      double rTmp,iTmp,jTmp,kTmp;
      if (l == 0.0)
       {
       re = q.re/l; x = q.x/l; y = q.y/l; z = q.z/l;
       rTmp = Double.POSITIVE_INFINITY; iTmp = 0.0; jTmp = 0.0; kTmp = 0.0;
       }
      else { rTmp = qq.re/l; iTmp = -qq.x/l; jTmp = -qq.y/l; kTmp = -qq.z/l; }
      assignTimes(q,new Quaternion(rTmp,iTmp,jTmp,kTmp));
    }
    /** @deprecated this = q/aR
     * @see #assignDivide( Quaternion q, double aR )
     * @return this
    */
    public final Quaternion setDivide( final Quaternion q, final double aR) {
      re = q.re/aR; x = q.x/aR; y = q.y/aR; z = q.z/aR;
      return this;
    }
    /**  Assigns <code>this</code> with a product of division between <code>q</code> and <code>aR</code>
     * @param <code>q</code> as dividend
     * @param <code>aR</code> as divisor
     */
    public final void assignDivide( final Quaternion q, final double aR) {
      re = q.re/aR; x = q.x/aR; y = q.y/aR; z = q.z/aR;
    }


    /** Returns C which is conjugated <code>this</code>.
     * creates a new <code>Quaternion</code> C which is conjugated <code>this</code>.
     * C = bar(this)
     */
    public final Quaternion conjugate() {
      return new Quaternion(re,-x,-y,-z);
    }
    /** Returns C which is conjugated <code>q</code>.
     * creates a new <code>Quaternion</code> C which is conjugated <code>q</code>.
     * C = bar(q)
     */
    public final static Quaternion conjugate( final Quaternion q) {
      return new Quaternion(q.re,-q.x,-q.y,-q.z);
    }
    /** @deprecated this = bar(this)
     * @see #assignConjugate()
     * @return this
     */
    public final Quaternion setConjugate () {
      x = -x; y = -y; z = -z;
      return this;
    }
    /** conjugates <code>this</code>: this = bar(this).
     * assigns <code>this</code> with conjugated <code>this</code>: this = bar(this)
     */
    public final void assignConjugate () {
      x = -x; y = -y; z = -z;
    }
     /** @deprecated this = bar(q)
     * @see #assignConjugate( Quaternion q)
     * @return this
     */
    public final Quaternion setConjugate( final Quaternion q) {
      re = q.re; x = -q.x; y = -q.y; z = -q.z;
      return this;
    }
     /**  Assigns <code>this</code> to conjugated <code>q</code>:
      * this = bar(q)
      */
    public final void assignConjugate( final Quaternion q) {
      re = q.re; x = -q.x; y = -q.y; z = -q.z;
    }
    /** BAR = CONJUGATE Returns C = bar(this).
     ** creates a new <code>Quaternion</code> C = bar(this) */
    public final Quaternion bar() {
      return new Quaternion(re,-x,-y,-z);
    }
    /** Returns C = bar(q)
     * Creates a new <code>Quaternion</code> C = bar(q) */
    public final static Quaternion bar( final Quaternion q) {
      return new Quaternion(q.re,-q.x,-q.y,-q.z);
    }
    /** @deprecated this = bar(this) */
    public final Quaternion setBar() {
       x = -x; y = -y; z = -z;
      return this;
    }
    /** Assigns <code>this = bar(this)</code> */
    public final void assignBar() {
       x = -x; y = -y; z = -z;
    }
     /**  @deprecated this = bar(q) */
    public final Quaternion setBar( final Quaternion q) {
      re = q.re; x = -q.x; y = -q.y; z = -q.z;
      return this;
    }
     /** Assigns <code>this = bar(q)</code> */
    public final void assignBar( final Quaternion q) {
      re = q.re; x = -q.x; y = -q.y; z = -q.z;
    }


    /** Returns C = -this.
     *	creates a new <code>Quaternion</code> C which is nedative to <code>this</code>
     * @return C = -this
     */
    public final Quaternion neg() {
      return new Quaternion(-re,-x,-y,-z);
    }
     /** Returns C = -q.
     *	creates a new <code>Quaternion</code> C which is nedative to <code>q</code>
     * @return C = -this
     */
    public final static Quaternion neg( final Quaternion q) {
      return new Quaternion(-q.re,-q.x,-q.y,-q.z);
    }
    /** @deprecated this = -this
     * @see #assignNeg()
     * @return this
     */
    public final Quaternion setNeg() {
      re = -re; x = -x; y = -y; z = -z;
      return this;
    }
    /** assigns <code>this</code> with negative sign.
     * this = -this
     */
    public final void assignNeg() {
      re = -re; x = -x; y = -y; z = -z;
    }
     /** @deprecated  this = -q
      * @see #assignNeg( Quaternion q)
      * @return this
      */
    public final Quaternion setNeg( final Quaternion q) {
      re = -q.re; x = -q.x; y = -q.y; z = -q.z;
      return this;
    }
     /**  Assigns <code>this</code> with negative <code>q</code>.
      * this = -q.
      */
    public final void assignNeg( final Quaternion q) {
      re = -q.re; x = -q.x; y = -q.y; z = -q.z;
    }


    /** In case of |q|==0 we set norm(q) = (0,0,0,0)
     ** creates a new Quaternion C = this/length(this) */
    public final Quaternion norm() {
      double l=length();
      if(l==0.0) l=1.0;
      return new Quaternion(re/l,x/l,y/l,z/l);
    }
     /** Returns C = q/length(q).
      * creates a new <code>Quaternion</code> C = q/length(q)
      */
    public final static Quaternion norm( final Quaternion q) {
      double l=length(q);
      if(l==0.0) l=1.0;
      return new Quaternion(q.re/l,q.x/l,q.y/l,q.z/l);
    }
    /** @deprecated this = this/length(this)
     * @see #assignNorm()
     * @return this
     */
    public final Quaternion setNorm() {
      double l=length();
      if(l==0.0) l=1.0;
      re /= l; x /= l; y /= l; z /= l;
      return this;
    }
    /** Assigns this = this/length(this) */
    public final void assignNorm() {
      double l=length();
      if(l==0.0) l=1.0;
      re /= l; x /= l; y /= l; z /= l;
    }
     /** @deprecated  this = q/length(q)
     * @see #assignNorm( Quaternion q)
     * @return this
     */
    public final Quaternion setNorm( final Quaternion q) {
      double l=length(q);
      if(l==0.0) l=1.0;
      re = q.re/l; x = q.x/l; y = q.y/l; z = q.z/l;
      return this;
    }
     /**  Assigns this = q/length(q) */
    public final void assignNorm(Quaternion q) {
      double l=length(q);
      if(l==0.0) l=1.0;
      re = q.re/l; x = q.x/l; y = q.y/l; z = q.z/l;
    }


     /** Returns C = (p+q)/2.
      * creates a new <code>Quaternion</code> C = (p+q)/2 */
    public final static Quaternion mean( final Quaternion p, final Quaternion q) {
      return new Quaternion(0.5*(p.re+q.re),0.5*(p.x+q.x),0.5*(p.y+q.y),0.5*(p.z+q.z));
    }
     /**  @deprecated this = (p+q)/2
     * @see #assignMean( Quaternion p, Quaternion q)
     * @return this
     */
    public final Quaternion setMean( final Quaternion p, final Quaternion q) {
      re = 0.5*(p.re+q.re); x = 0.5*(p.x+q.x);
      y = 0.5*(p.y+q.y); z = 0.5*(p.z+q.z);
      return this;
    }
     /**  Assigns this = (p+q)/2 */
    public final void assignMean( final Quaternion p, final Quaternion q) {
      re = 0.5*(p.re+q.re); x = 0.5*(p.x+q.x);
      y = 0.5*(p.y+q.y); z = 0.5*(p.z+q.z);
    }



     /**  Returns <p,q>_R^4 */
    public final double dot( final Quaternion p, final Quaternion q) {
      return p.re*q.re + p.x*q.x + p.y*q.y * p.z+q.z;
    }
     /**  Returns <this,p>_R^4 */
    public final double dot( final Quaternion p) {
      return re*p.re + x*p.x + y*p.y * z+p.z;
    }



    /** Returns <code>r</code> of <code>q</code>: r=q.r */
    public final static double r( final Quaternion q) {
	return q.re;
    }
    /** Returns <code>i</code> of <code>q</code>: i = q.i */
    public final static double i( final Quaternion q) {
	return q.x;
    }
    /** Returns <code>j</code> of <code>q</code>: j = q.j */
    public final static double j( final Quaternion q) {
	return q.y;
    }
    /** Returns <code>k</code> of <code>q</code>: k = q.k */
    public final static double k( final Quaternion q) {
	return q.z;
    }

    /** Returns <code>this</code> as string */
    public final String toString() {
      StringBuffer sb=new StringBuffer().append('(').append(re);
      if(x>=0.0) sb.append('+');
      sb.append(x).append('i');
      if(y>=0.0) sb.append('+');
      sb.append(y).append('j');
      if(z>=0.0) sb.append('+');
      return sb.append(z).append('k').append(')').toString();
    }
    /**Indicates whether some other object is "equal to" this one.*/
    public final boolean equals( final Object o) {
      try {
        Quaternion c=(Quaternion)o;
        return abs( minus( c ) ) < Double.MIN_VALUE;
      } catch(ClassCastException ex) {
        return false;
      }
    }

    /** @deprecated
      * @see #assignRandom()
      * @return this
      */
    public final Quaternion setRandom()
    {
       re = Math.random();
       x = Math.random();
       y = Math.random();
       z = Math.random();

       return this;
    }
    /** Assigns <code>this = random()</code> */
    public final void assignRandom()
    {
       re = Math.random();
       x = Math.random();
       y = Math.random();
       z = Math.random();

    }
    /** Assigns this = (d[offset], ... , d[offset+3])
      * returns <code>this</code>                     */
    public final Quaternion set( final double[] d, final int offset ) {
      re=d[offset]; x=d[1+offset]; y=d[2+offset]; z=d[3+offset]; return this;
    }
    /** Assigns this = (d[offset], ... , d[offset+3]) */
    public final void assign( final double[] d, final int offset ) {
      re=d[offset]; x=d[1+offset]; y=d[2+offset]; z=d[3+offset];
    }
    /** Assigns (d[offset], ... , d[offset+3]) = this */
    public final Quaternion get( final double[] d, final int offset ) {
      d[offset]=re; d[1+offset]=x; d[2+offset]=y; d[3+offset]=z; return this;
    }

    /** Returns C = -bar(this)<sup>-1</sup>.
      * creates a new <code>Quaternion</code> C = -bar(this)<sup>-1</sup> */
    public final Quaternion dual() {
      return dual(this);
    }
    /** Returns C = -bar(q)<sup>-1</sup>.
      * Creates a new <code>Quaternion</code> C = -bar(q)<sup>-1</sup> */
    public final static Quaternion dual( final Quaternion q ) {
      double l=-q.lengthSqr();
      if (l==0.0) return new Quaternion(Double.POSITIVE_INFINITY,0.0,0.0,0.0);
      return new Quaternion(q.re/l,q.x/l,q.y/l,q.z/l);
    }
    /** @deprecated this = -bar(this)<sup>-1</sup> */
    public final Quaternion setDual() {
      assignDual();
      return this;
    }
    /** Assigns <code>this = -bar(this)<sup>-1</sup></code> */
    public final void assignDual() {
      double l=-lengthSqr();
      if (l==0.0) { re=Double.POSITIVE_INFINITY; x=0.0; y=0.0; z=0.0; }
      else { re/=l; x/=l; y/=l; z/=l; }
    }
     /**  @deprecated this = -bar(q)<sup>-1</sup> */
    public final Quaternion setDual( final Quaternion q) {
      assignDual(q);
      return this;
    }
     /** Assigns <code>this = -bar(q)<sup>-1</sup></code> */
    public final void assignDual( final Quaternion q) {
      double l=-q.lengthSqr();
      if (l==0.0) { re=Double.POSITIVE_INFINITY; x=0.0; y=0.0; z=0.0; }
      else { re=q.re/l; x=q.x/l; y=q.y/l; z=q.z/l; }
    }

    /** Returns r = |u-v|^2 */
    public final static double distSqr( final Quaternion u, final Quaternion v ) {
	double subR = u.re - v.re;
	double subI = u.x - v.x;
	double subJ = u.y - v.y;
	double subK = u.z - v.z;

	return subR*subR + subI*subI + subJ*subJ + subK*subK;
    }

    /** Returns r = |this-v|^2 */
    public final double distSqr( final Quaternion v ) {
	double subR = re - v.re;
	double subI = x - v.x;
	double subJ = y - v.y;
	double subK = z - v.z;

	return subR*subR + subI*subI + subJ*subJ + subK*subK;
    }

    /** Returns r = |u-v| */
    public final static double dist( final Quaternion u, final Quaternion v ) {
	return Math.sqrt( distSqr( u, v ) );
    }

    /** Returns r = |this-v| */
    public final double dist( final Quaternion v ) {
	return Math.sqrt( this.distSqr( v ) );
    }



    /**Indicates whether some other quaternion is "equal to" this one.*/
    public final boolean equals( final Quaternion c )
    {
      return this == c ? true : equals(c.re, c.x, c.y, c.z);
    }

    /**Indicates whether re + ix + iy + iz is "equal to" this one according to eps.*/
    public final boolean equals(final double re, final double x, final double y,
                               final double z, final double eps ) {
     return (Math.abs(this.re - re) < eps)
         && (Math.abs(this.x - x) < eps) && (Math.abs(this.y - y) < eps) &&
         (Math.abs(this.z - z) < eps);
   }

    /**Indicates whether re + ix + iy + iz is is "equal to" this one.*/
    public final boolean equals(final double re, final double x, final double y,
                                final double z) {
      return equals( re, x, y, z, EPS );
    }

  /**Indicates whether <code>q</code> equal to <code>this</code>*/
    public final boolean equals( Quaternion q, final double eps){
      return equals( q.re, q.x, q.y, q.z );
    }

  /**Indicates whether <code>c1</code> equal to <code>c2</code>*/
    public static final boolean arraysCoincide( final Complex[] c1, final Complex[] c2, final double eps ){
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
    public static final boolean arraysCoincide( final Complex[] c1, final Complex[] c2){
      return Complex.arraysCoincide(c1, c2, Complex.EPS );
    }

     /**Indicates whether <code>c1</code> equal to <code>c2</code> as sets.*/
        public static final boolean setsCoincide( final Quaternion[] c1, final Quaternion[] c2 ){
          return Quaternion.setsCoincide( c1, c2, Complex.EPS );
        }

 /**Indicates whether <code>c1</code> equal to <code>c2</code> as sets.*/
    public static final boolean setsCoincide( final Quaternion[] c1, final Quaternion[] c2,
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
		if( c1[i].equals( (Quaternion)it.next(),eps )) {
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
	The methods creates instances of <code>Field.Complex</code> if necessary
	and handles further more the case that source and desination
	overlap.
	Attention! The case that instances appear in both arrays <code>s</code>
	and <code>d</code> yields to unpredictable results.
    */
    public static final void arraycopy( final Quaternion [] s, int is,
					final Quaternion [] d, int id, final int size ) {

	if( !( is+size <= s.length && id+size <= d.length ) )
	    throw new IndexOutOfBoundsException();

	if( s == d && id == is )
	    return;

	if( id <= is ) {

	    for( int i=0; i<size; i++, is++, id++ ) {

		if( s[is] == null )
		    d[id] = null;
		else if( d[id] != null )
		    d[id].assign( s[is] );
		else
		    d[id] = new Quaternion( s[is] );
	    }

	} else {

	    is += size-1;
	    id += size-1;

	    for( int i=0; i<size; i++, is--, id-- ) {

		if( s[is] == null )
		    d[id] = null;
		else if( d[id] != null )
		    d[id].assign( s[is] );
		else
		    d[id] = new Quaternion( s[is] );
	    }
	}
    }

    /** Returns copy of array <code>s</code>.
	The methods copies the values of <code>s</code>
	not only the references, therefore it differs
	from <code>clone()</code>.
    */
    public static final Quaternion[] copy( final Quaternion [] s ){

	final Quaternion [] d = new Quaternion[s.length];

	arraycopy( s, 0, d, 0, s.length );

	return d;
    }

}
