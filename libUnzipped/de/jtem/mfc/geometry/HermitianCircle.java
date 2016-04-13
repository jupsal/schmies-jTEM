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

package de.jtem.mfc.geometry;

import java.io.Serializable;

import de.jtem.mfc.field.Complex;
import de.jtem.mfc.group.Moebius;
import de.jtem.mfc.matrix.Complex2By2;
import de.jtem.mfc.matrix.HermitianComplex2By2;
import de.jtem.mfc.vector.Real3;

/**
 * <p>A class to represent <strong>oriented circles in 2-dimensional Moebius
 * geometry</strong>.</p>
 *
 * <p>An oriented circle is represented by a hermitian complex
 * <code>2&times;2</code>-matrix with non-positive determinant.</p>
 *
 * <p>Matrices which differ by a positive factor represent the same
 * oriented circle. </p>
 *
 * <p>Matrices which differ by a negative factor represent the same
 * circle with opposite orientation. </p>
 */
public class HermitianCircle extends HermitianComplex2By2 
    implements Serializable, Cloneable {

    private static final long serialVersionUID = 1L;

    private static final Real3 dummyReal3 = new Real3();
    private static final Complex dummyComplexA = new Complex();
    private static final Complex dummyComplexB = new Complex();
    private static final Moebius dummyMoebius = new Moebius();
    private static final Complex2By2 dummyComplex2By2 = new Complex2By2();

    /**
     * Create unit circle in plane (= equator in Rieman sphere).
     */ 
    public HermitianCircle() {
        super(1.0, 0.0, 0.0, -1.0);
    }

    public HermitianCircle(HermitianCircle c){
        super(c);
    }

    /**
     * Set the value of aRe.
     * @param v  value to assign to aRe
     * @throws IllegalArgumentException if determinant would become positive
     */
    public void setARe(double v) throws IllegalArgumentException {
        if (v * dRe - bRe * bRe - bIm * bIm <= 0.0) {
            super.setARe(v);
        } else {
            throw new IllegalArgumentException("Determinant would become " +
                                               "positive. ");
        }
    }
    
    /**
     * Set the value of the upper right element. Also sets the lower left
     * element to the complex conjugate.
     *
     * @param b a {@link Field.Complex}
     * @throws IllegalArgumentException if determinant would become positive
     */
    public void setB(final Complex  b) throws IllegalArgumentException {
        if (aRe * dRe - b.re * b.re - b.im * b.im <= 0.0) {
            super.setB(b);
        } else {
            throw new IllegalArgumentException("Determinant would become " +
                                               "positive. ");
        }
    }
    
    /**
     * Set the value of the lower left element. Also sets the upper right
     * element to the complex conjugate.
     *
     * @param c a {@link Field.Complex}.
     * @throws IllegalArgumentException if determinant would become positive
     */
    public void setC(final Complex c) throws IllegalArgumentException {

        if (aRe * dRe - c.re * c.re - c.im * c.im <= 0.0) {
            super.setC(c);
        } else {            
            throw new IllegalArgumentException("Determinant would become " +
                                               "positive. ");
        }
    }

    /**
     * Set the value of dRe.
     * @param v  Value to assign to dRe.
     * @throws IllegalArgumentException if determinant would become positive
     */
    public void setDRe(double v) throws IllegalArgumentException {
        if (aRe * v - bRe * bRe - bIm * bIm <= 0.0) {
            super.setDRe(v);
        } else {
            throw new IllegalArgumentException("Determinant would become " +
                                               "positive. ");
        }        
    }

    /** <p> Assign by prescribing 3 points on the circle. </p>
     *
     *  
     *
     *
     */
    public void assignThrough(ComplexProjective1 z1, ComplexProjective1 z2, ComplexProjective1 z3){

        assign(0, 0, 1, 0);

        dummyMoebius.assign(z1, z2, z3);
        apply(dummyMoebius);
    }  

    /** 
     * <p>Assign by Lorentz vector.</p> 
     *
     * <p>The arguments must satisfy &emsp;<code>x0<sup>2</sup> -
     * x1<sup>2</sup> - x2<sup>2</sup> - x3<sup>2</sup> <=
     * 0.0</sup></code></p>
     *
     * @throws IllegalArgumentException if arguments do not satisfy above
     * condition
     */
    public void assignLorentzVector(final double x0,
                                    final double x1,
                                    final double x2,
                                    final double x3) 
        throws IllegalArgumentException{
        super.assignLorentzVector(x0, x1, x2, x3);
        if (realDeterminant() > 0.0) {
            throw new IllegalArgumentException("Lorentz vector does not " +
                                               "represent a circle. ");
        }
    }

    /**
     * Test for equality.
     *
     * @param circ the circle to compare <code>this</code> with 
     * @return <code>true</code> if <code>this</code> and <code>circ</code>
     * represent the same oriented circle, <code>false</code>
     * otherwise.
     */
    public boolean equals(HermitianCircle circ) {
        final double x0 = this.lorentzVectorX0();
        final double x1 = this.lorentzVectorX1();
        final double x2 = this.lorentzVectorX2();
        final double x3 = this.lorentzVectorX3();
        final double y0 = circ.lorentzVectorX0();
        final double y1 = circ.lorentzVectorX1();
        final double y2 = circ.lorentzVectorX2();
        final double y3 = circ.lorentzVectorX3();

        final double xx = x0 * x0 + x1 * x1 + x2 * x2 + x3 * x3;
        final double yy = y0 * y0 + y1 * y1 + y2 * y2 + y3 * y3;
        final double xy = x0 * y0 + x1 * y1 + x2 * y2 + x3 * y3;

        /* The following expression is homogenous in x and y. It is
         * non-negative, and equal to 0 if and only if x and y are linearly
         * dependent (Schwarz inequality).
         */
        final double h = 1.0 - xy * xy /(xx * yy);
            
        if (h > EPS) {
            return false;
        } else {
            // Check whether circles are equally oriented.

            // Find x which is largest in abs-value and corresponding y.
            double x;
            double y;
            if (Math.abs(x1) > Math.abs(x0)) {
                x = x1;
                y = y1;
            } else {
                x = x0;
                y = y0;
            }
            if (Math.abs(x2) > Math.abs(x)) {
                x = x2;
                y = y2;
            }
            if (Math.abs(x3) > Math.abs(x)) {
                x = x3;
                y = y3;
            }
            
            if (x * y > 0.0) {
                return true;
            } else {
                return false;
            }
        }
    }

    /**
     * Assign by spherical center and spherical radius.
     *
     * @param center a {@link ComplexProjective1}: the spherical center
     * @param radius a <code>double</code>: the spherical radius
     */ 
    public void assignSpherical(final ComplexProjective1 center, 
                                final double radius) {
        final double x0;
        final double x1;
        final double x2;
        final double x3;
        center.projectTo(dummyReal3);
        x0 = Math.cos(radius);
        x1 = dummyReal3.x;
        x2 = dummyReal3.y;
        x3 = dummyReal3.z;
        assignLorentzVector(x0, x1, x2, x3);
    }

    /** 
     * Return the spherical radius.
     *
     * @return a <code>double</code> between <code>0</code> and
     * <code>&pi;</code>: the spherical radius
     */
    public double sphericalRadius() {
        final double x0 = lorentzVectorX0();
        final double x1 = lorentzVectorX1();
        final double x2 = lorentzVectorX2();
        final double x3 = lorentzVectorX3();
        final double n = Math.sqrt(x1 * x1 + x2 * x2 + x3 * x3);
        return Math.acos(x0 / n);
    }

    /**
     * Calculate the spherical center in real 3-space. Since the circles are
     * oriented, the spherical center is unique.
     *
     * @param center a {@link Real3} used to hand over the result
     */
    public void sphericalCenter(final Real3 center) {
        final double x1 = lorentzVectorX1();
        final double x2 = lorentzVectorX2();
        final double x3 = lorentzVectorX3();
        center.assign(x1, x2, x3);
        center.normalize();
    }

    /**
     * Calculate the spherical center in 1-dimensional complex projective
     * space. Since the circles are oriented, the spherical center is unique.
     *
     * @param center a {@link ComplexProjective1} used to hand over the
     * result
     */
    public void sphericalCenter(final ComplexProjective1 center) {
        sphericalCenter(dummyReal3);
        dummyReal3.projectTo(center);
    }

    /**
     * Assign by euclidean center and euclidean radius.
     *
     * @param center a {@link ComplexProjective1}: the Euclidean center
     *
     * @param r a <code>double</code>: the euclidean radius
     */ 
    public void assignEuclidean(final ComplexProjective1 center,
                                final double radius) {

        final double newA;
        final double newD;

        center.getA(dummyComplexA);
        center.getB(dummyComplexB);

        newA = dummyComplexB.absSqr();
        newD = dummyComplexA.absSqr() - newA * radius * radius;

        dummyComplexB.assignConjugate();

        dummyComplexB.assignTimes(dummyComplexA);
        
        dummyComplexB.assignTimes(-1.0);

        this.assign(newA, dummyComplexB, newD);

        if (radius < 0.0) {
            this.assignTimes(-1.0);
        }

    }

    /**
     * Assign by euclidean center and euclidean radius.
     *
     * <p>TODO: maybe one should rather implement this using {@link #assignEuclidean(ComplexProjective1, double)}?
     *
     * @param m1 a <code>double</code>: the x coordinate of the euclidean center
     * @param m2 a <code>double</code>: the y coordinate of the euclidean center
     * @param r a <code>double</code>: the euclidean radius
     */ 
    public void assignEuclidean(final double m1, final double m2, 
                                final double r) {
        //         final double x0;
        //         final double x1;
        //         final double x2;
        //         final double x3;
        
        assign(1.0, -m1, -m2, (m1 * m1 + m2 * m2 - r * r));

        // normalize determinant if possible
        if (r != 0.0) {
            assignDivide(r);
        }
        //     if (r != 0.0) {
        //         assign(1/r, -m1/r, -m2/r, (m1 * m1 + m2 * m2 - r * r)/r);
        //     } 
        //     else {
        //         assign(1.0, -m1, -m2, (m1 * m1 + m2 * m2 - r * r));
        //     }

        // hier ist`s noch nicht richtig
        //       if (Double.isInfinite (m1) && Double.isInfinite (m2)) {
        //           x0 = 1;
        //           x1 = -1;
        //           x2 = 1;
        //           x3 = 1;
        //       } 
        //       else {
        //           x0 = (1 + (m1 * m1 + m2 * m2));
        //           x1 = (1 - (m1 * m1 + m2 * m2));
        //           x2 = 2 * m1;
        //           x3 = 2 * m2;
        //       } 
        //        assignLorentzVector(x0, x3, x2, x1);
    }

    /** 
     * Return the euclidean radius. 
     * The result may be negative according to the orientation of the circle
     *
     * @return a <code>double</code>
     */
    public double euclideanRadius() {
        
        return Math.sqrt(-realDeterminant())/aRe;

        //          if (realDeterminant() == 0){
        //              return 0;
        //          }
        //          else{
        //              return Math.sqrt(-realDeterminant()/aRe);
        //          }
        //         return Math.sqrt(-realDeterminant())/aRe;
    }

    /**
     * Calculate the euclidean center in the complex plane.
     *
     * @param center a {@link Field.Complex} used to hand over the result
     */
    public void euclideanCenter(final Complex center){

        center.assign(-bRe, -bIm);

        //         if (realDeterminant() == 0){
        //             center.assignTimes(realDeterminant()* aRe);
        //         }
        //         else{
        //             center.assignDivide(aRe);
        //         }
        
        center.assignDivide(aRe);

    }

    /**
     * Calculate the euclidean center in the complex plane.
     *
     * @param center a {@link ComplexProjective1} used to hand over the result
     */
    public void euclideanCenter(final ComplexProjective1 center){

        euclideanCenter(dummyComplexA);
        dummyComplexB.assign(1.0, 0.0);

        //         if (realDeterminant() == 0){
        //             center.assignTimes(realDeterminant()* aRe);
        //         }
        //         else{
        //             center.assignDivide(aRe);
        //         }
        
        center.assign(dummyComplexA, dummyComplexB);

    }

    /**
     * Assign by hyperbolic center and hyperbolic radius.
     *
     * @param m1 a <code>double</code>: the x coordinate of the hyperbolic
     * center (in the Poincare unit disc model)
     * @param m2 a <code>double</code>: the y coordinate of the hyperbolic
     * center (in the Poincare unit disc model)
     * @param r a <code>double</code>: the hyperbolic radius
     * @throws {@link IllegalArgumentException} if m1 * m1 + m2 * m2 >= 1
     */ 
    public void assignHyperbolic(final double m1, final double m2, 
                                 final double r) {
        final double nsq = m1 * m1 + m2 * m2;
        if (nsq >= 1.0) {
            throw new IllegalArgumentException("The point (" + m1 + ", " + m2 
                                               + ") does not lie inside the "
                                               + "unit circle. ");
        }
        final double er = Math.exp(r);
        final double s = 0.25 * (er + 1.0/er - 2.0);  // s = sinh^2 (r/2)
        final double ddd = 1.0 / (1.0 - nsq);
        final double newA = ddd + s;
        final double newD = nsq * ddd - s;
        final double newBre = - m1 * ddd;
        final double newBim = - m2 * ddd;
        assign(newA, newBre, newBim, newD);
    }

    /**
     * Assign by hyperbolic center and hyperbolic radius.
     *
     * @param z a {@link de.jtem.mfc.field.Field.Complex}: the hyperbolic center (in the
     * Poincare unit disc model)
     * @param r a <code>double</code>: the hyperbolic radius
     * @throws {@link IllegalArgumentException} if |z| >= 1
     */ 
    public void assignHyperbolic(final Complex z, final double r) {
        assignHyperbolic(z.re, z.im, r);
    }

    /**
     * Assign by hyperbolic center and hyperbolic radius.
     *
     * @param z a {@link ComplexProjective1}: the hyperbolic center (in the
     * Poincare unit disc model)
     * @param r a <code>double</code>: the hyperbolic radius
     * @throws {@link IllegalArgumentException} if <code>z</code> does not
     * lie inside the unit circle
     */ 
    public void assignHyperbolic(final ComplexProjective1 z, final double r) {
        z.projectTo(dummyComplexA);
        assignHyperbolic(dummyComplexA, r);
    }

    

    /**
     * Assign <code>this</code> with <code>m<sup>-1</sup> this m<sup>-1 *</sup></code>.
     * 
     * @param m an {@link Moebius}: the Moebius transform to apply 
     */
    public void apply(Moebius m){
        
        dummyMoebius.assignInvert(m);
        dummyMoebius.assignAdjoined();
        this.assignAdjoinedWith(dummyMoebius);

    } 
    /**
     * For test purposes:
     * Assign <code>this</code> with <code>m<sup>-1</sup> this m<sup>-1 *</sup></code>.
     * 
     * @param m an {@link Complex2By2}: the Moebius transform to apply 
     */
    public void apply(Complex2By2 m){
        
        dummyComplex2By2.assignInvert(m);
        dummyComplex2By2.assignAdjoined();
        this.assignAdjoinedWith(dummyComplex2By2);

    } 
}
