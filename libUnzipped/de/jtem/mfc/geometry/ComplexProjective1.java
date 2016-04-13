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
import de.jtem.mfc.vector.Complex2;
import de.jtem.mfc.vector.Real3;

/**
 * <p>A class that represents points in 1-dimensional complex projective
 * space, also known as the Riemann sphere.</p>
 *
 * <p>The elements of a {@link Complex2} are interpreted as homogeneous
 * coordinates. </p>
 */
public class ComplexProjective1 extends Complex2 implements Serializable, 
                                                            Cloneable {

    private static final long serialVersionUID = 1L;

    /**
     * Creates the point with homogeneous coordinates <code>(1,0)</code>.
     * This is the north pole of the Riemann sphere.
     */
    public ComplexProjective1() {
        super(1.0, 0.0, 0.0, 0.0);
    }

    /** 
     * Creates the point with homogenous coordinates <code>(a, b)</code>.
     */
    public ComplexProjective1(final Complex a, final Complex b) {
        super(a, b);
    }

    /**
     * Creates a copy of p.
     */
    public ComplexProjective1(ComplexProjective1 p) {
        super(p);
    }
    
    /** 
     * Creates the point with homogenous coordinates <code>(a, b)</code>,
     * where <code>a = theARe + i theAIm</code> and <code>b = theBRe + i
     * theBIm</code>.
     */
    public ComplexProjective1(final double theARe, final double theAIm, 
                              final double theBRe, final double theBIm) {
        super(theARe, theAIm, theBRe, theBIm);
    }
    
    /**
     * Apply a Moebius transformation.
     */
    public void apply(Moebius m) {
        m.applyTo(this);
    }

    /**
     * <p> Project to the complex plane. </p>
     *
     * <p> The point with coordinates <code>(1,0)</code> is sent to
     * infinity.<br> The point with coordinates <code>(0,1)</code> is
     * sent to zero.</p>
     *
     * <p> Simply assigns <code>a / b</code> to
     * <code>z</code>. </p>
     *
     * @param z is set to <code>a/b</code>.
     */
    final public void projectTo(final Complex z) {
        z.assign(aRe, aIm);
        z.assignDivide(bRe, bIm);
    }

    /**
     * <p> Project stereographically to the sphere in real 3-space. </p>
     *
     * <p> The point with coordinates <code>(1,0)</code> is sent to the
     * north pole <code>(0,0,1)</code>.<br> The point with coordinates
     * <code>(0,1)</code> is sent to the south pole <code>(0,0,-1)</code>.
     *
     * <p> The components of <code>v</code> are set to <br>
     * <code> 
     *    v.x = 2 Re(a conjugate(b)) / (|a|<sup>2</sup> + |b|<sup>2</sup>)
     *    <br>
     *    v.y = 2 Im(a conjugate(b)) / (|a|<sup>2</sup> + |b|<sup>2</sup>)
     *    <br>
     *    v.z = (|a|<sup>2</sup> - |b|<sup>2</sup>) / 
     *          (|a|<sup>2</sup> + |b|<sup>2</sup>)
     * </code>.
     * </p>
     *
     * @param v a {@link Real3} used to hand over the result. Must not be
     * <code>null</code>.
     */
    final public void projectTo(final Real3 v) {

        /* Calculate norm of a and b and their sum */
        final double aNormSq = aRe * aRe + aIm * aIm;
        final double bNormSq = bRe * bRe + bIm * bIm;
        final double normSqSum = aNormSq + bNormSq;

        /* Calculate real and imaginary parts of (2 times a times conjugate
         * of b) */
        final double twoAConjugateBRe = 2.0 * (aRe * bRe + aIm * bIm);
        final double twoAConjugateBIm = 2.0 * (aIm * bRe - aRe * bIm);

        v.x = twoAConjugateBRe / normSqSum;
        v.y = twoAConjugateBIm / normSqSum;
        v.z = (aNormSq - bNormSq) / normSqSum;
    }

    public static final boolean equals (final double aRe1, final double aIm1,
                                        final double bRe1, final double bIm1,
                                        final double aRe2, final double aIm2,
                                        final double bRe2, final double bIm2){

        return (Math.abs((aRe1 * bRe2 - aIm1 * bIm2)
                         - (aRe2 * bRe1 - aIm2 * bIm1)) < EPS &&
                Math.abs((aIm1 * bRe2 + aRe1 * bIm2)
                         - (aIm2 * bRe1 + aRe2 * bIm1)) < EPS );

    }  
    
    public static final boolean equals (final double aRe1, final double aIm1,
                                        final double bRe1, final double bIm1,
                                        final double aRe2, final double aIm2){

        return equals (aRe1, aIm1, bRe1, bIm1, aRe2, aIm2, 1, 0); 

    }
    
    public final boolean equals (final double aRe1, final double aIm1,
                                        final double bRe1, final double bIm1){

        return equals (aRe1, aIm1, bRe1, bIm1, aRe, aIm, bRe, bIm); 

    }

    /**Indicates whether some other object is "equal to" this one.*/
    public final boolean equals( final Object o )
    {
	try
	{
            return equals((ComplexProjective1) o);
	}
	catch(ClassCastException ex)
	{
	    return false;
	}
    }   

    /**Indicates whether some other ComplexProjective1 is "equal to" this one.*/
    public final boolean equals( final ComplexProjective1 c )
    {
	    return this == c ? true : equals(c.aRe, c.aIm, c.bRe, c.bIm, 
                                             aRe, aIm, bRe, bIm);
    }

}
