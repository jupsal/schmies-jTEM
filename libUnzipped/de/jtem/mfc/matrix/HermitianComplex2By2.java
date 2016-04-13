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

/**
 * <p>A class to represent Hermitian complex
 * <code>2&times;2</code>-matrices.</p>
 *
 * <p>A complex <code>2&times;2</code>-matrix is Hermitian if the diagonal
 * entries are real and the off diagonal entries are complex conjugate to
 * each other.</p>
 * 
 */
public class HermitianComplex2By2 extends AbstractComplex2By2 
                                  implements Serializable, Cloneable {
    
    private static final long serialVersionUID = 1L;

    /**
     * Create the identity matrix.
     */
    public HermitianComplex2By2() {
        super();
    }

    /**
     * Create a matrix with the same data as m.
     */
    public HermitianComplex2By2(HermitianComplex2By2 m) {
        super(m);
    }

    /**
     * Create a hermitian matrix with given diagonal elements and given
     * upper right element.
     *
     * @param theA   the upper left entry
     * @param theBRe the real part of the upper right entry
     * @param theBIm the imaginary part of the upper right left entry
     * @param theD   the lower right entry
     */
    public HermitianComplex2By2(final double theA, 
                                final double theBRe, final double theBIm, 
                                final double theD) {
        super(theA, 0.0, theBRe, theBIm, theBRe, -theBIm, theD, 0.0);
    }

    /**
     * Create a hermitian matrix with given diagonal elements and given
     * upper right element.
     *
     * @param a a <code>double</code>: the real upper left entry
     * @param b a {@link Field.Complex}: the complex upper right entry
     * @param d a <code>double</code>: the real lower right entry
     */
    public HermitianComplex2By2(final double a, 
                                final Complex b, 
                                final double d) {
        this(a, b.re, b.im, d);
    }

    /**
     * Get the value of aRe.
     * @return value of aRe.
     */
    public double getARe() {
        return aRe;
    }
    
    /**
     * Set the value of aRe.
     * @param v  Value to assign to aRe.
     */
    public void setARe(double  v) {
        this.aRe = v;
    }
    
    /**
     * Set the value of the upper right element. Also sets the lower left
     * element to the complex conjugate.
     *
     * @param b a {@link Field.Complex}
     */
    public void setB(final Complex  b) {
        this.bRe = b.re;
        this.bIm = b.im;
        this.cRe = b.re;
        this.cIm = -b.im;
    }
    
    /**
     * Set the value of the lower left element. Also sets the upper right
     * element to the complex conjugate.
     *
     * @param c a {@link Field.Complex}.
     */
    public void setC(final Complex c) {
        this.cRe = c.re;
        this.cIm = c.im;
        this.bRe = c.re;
        this.bIm = -c.im;
    }

    /**
     * Get the value of dRe.
     * @return value of dRe.
     */
    public double getDRe() {
        return dRe;
    }
    
    /**
     * Set the value of dRe.
     * @param v  Value to assign to dRe.
     */
    public void setDRe(double  v) {
        this.dRe = v;
    }

    /** 
     * Assign hermitian matrix by diagonal and upper right elements.
     *
     * @param newA   the upper left entry
     * @param newBRe the real part of the upper right entry
     * @param newBIm the imaginary part of the upper right entry
     * @param newD   the lower right entry
     */
    public void assign(final double newA, 
                       final double newBRe, final double newBIm, 
                       final double newD) {


        /* this.aIm and this.dIm need not be assigned because they are zero
         * anyway. */
        this.aRe = newA;
        this.aIm = 0.0;
        this.bRe = newBRe;
        this.bIm = newBIm;
        this.cRe = newBRe;
        this.cIm = -newBIm;
        this.dRe = newD;
        this.dIm = 0.0;
    }

    /** 
     * Assign hermitian matrix by diagonal and upper right elements.
     *
     * @param a a <code>double</code>: the real upper left entry
     * @param b a {@link Field.Complex}: the complex upper right entry
     * @param d a <code>double</code>: the real lower right entry
     */
    public void assign(final double a, 
                       final Complex b, 
                       final double d) {
        assign(a, b.re, b.im, d);
    }
  
  /** assigns this with the prescribed matrix */
    public void assign( final HermitianComplex2By2 s ) { 
	assign( s.aRe, s.bRe, s.bIm, s.dRe); 
    }

     /**
     * <p>Assign <code>this</code> with <code>t this t<sup>-1</sup>.</p>
     *
     * <p>Simply calls {@link
     * AbstractComplex2By2#assignConjugateWith(AbstractComplex2By2)}.
     * 
     * @param t an {@link AbstractComplex2By2}: the matrix with which to
     * conjugate
     */
    public void assignConjugateWith(final AbstractComplex2By2 t) {
        super.assignConjugateWith(t);
    }

     /**
     * <p>Assign <code>this</code> with <code>t this t<sup>*</sup>.</p>
     *
     * <p>Simply calls {@link
     * AbstractComplex2By2#assignAdjoinedWith(AbstractComplex2By2)}.
     * 
     * @param t an {@link AbstractComplex2By2}: the matrix to adjoin with
     */
    public void assignAdjoinedWith(final AbstractComplex2By2 t) {
        super.assignAdjoinedWith(t);
    }

   /**
     * <p>Assign hermitian matrix by a Lorentz vector.</p>
     *
     * <p>Sets <code>this</code> to <br>
     *
     * <pre>
     *     [ x0 - x3       -x1 - i x2 ]
     *     [                          ]
     *     [ -x1 + i x2     x0 + x3   ].
     * </pre></p>
     *
     * <p> Note that the determinant is &emsp;<code>x0<sup>2</sup> -
     * x1<sup>2</sup> - x2<sup>2</sup> - x3<sup>2</sup></code>.</p>
     */
    public void assignLorentzVector(final double x0, final double x1,
                                    final double x2, final double x3) {
        assign(x0 - x3, -x1, -x2, x0 + x3);
    }
 
    /**
     * Return the <code>x0</code>-component of the Lorentz vector
     * corresponding to <code>this</code>.
     * 
     * @return <code>x0</code>
     * @see #assignLorentzVector(double, double, double, double)
     */
    final public double lorentzVectorX0() {
        return 0.5 * (aRe + dRe);
    }

    /**
     * Return the <code>x1</code>-component of the Lorentz vector
     * corresponding to <code>this</code>.
     * 
     * @return <code>x1</code>
     * @see #assignLorentzVector(double, double, double, double)
     */
    final public double lorentzVectorX1() {
        return -bRe;
    }

    /**
     * Return the <code>x2</code>-component of the Lorentz vector
     * corresponding to <code>this</code>.
     * 
     * @return <code>x2</code>
  e  * @senL#assigVeorentzouctor(double, double, double, d  ble)
     */
    final public double lorentzVectorX2() {
        return -bIm;
    }

    /**
     * Return the <code>x3</code>-component of the Lorentz vector
     * corresponding to <code>this</code>.
     * 
     * @return <code>x3</code>
     * @see #assignLorentzVector(double, double, double, double)
     */
    final public double lorentzVectorX3() {
        return 0.5 * (dRe - aRe);
    }

    /**
     * <p>Return the real determinant of <code>this</code>.</p>
     *
     * <p>The determinant of a hermitian matrix is real.</p>
     *
     * @return the real determinant <code>ad - bc</code>
     */
    public double realDeterminant() {
        return aRe * dRe - bRe * bRe - bIm * bIm;
    }

}
