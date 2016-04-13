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

package de.jtem.mfc.group;

import java.io.Serializable;

import de.jtem.mfc.field.Complex;
import de.jtem.mfc.field.ComplexConstant;
import de.jtem.mfc.field.Field;
import de.jtem.mfc.geometry.ComplexProjective1;
import de.jtem.mfc.matrix.AbstractComplex2By2;

/** <p>This class represents linear fractional transformations or M&ouml;bius
 * transformations. </p>
 *
 * <p align=center>
 * <code>
 * <table>
 * <tr> <td>             </td> <td> a z + b             </td> </tr>
 * <tr> <td>z |-> &nbsp; </td> <td> <hr noshade size=1> </td> </tr>
 * <tr> <td>             </td> <td> c z + d             </td> </tr>
 * </table>
 * </code></p>
 *
 * <p>
 * All methods of <code>Complex2By2</code> try to avoid the
 * creation of temporarily used objects in their internal computations,
 * unless the opposite is stated in their describtions.

 * <p>This class has <code>final</code> member objects, which are
 * used as dummy variables. Thus, it is definitely not thread-safe.

 * @author schmies
 **/

public class Moebius extends AbstractComplex2By2
  implements Complex.FunctionOnComplex, Serializable {


    private static final long serialVersionUID = 1L;

    final ComplexProjective1 dummyComplexProjective1 =
        new ComplexProjective1();

    public Moebius() {
    }

    public Moebius(final double ar, final double ai,
                   final double br, final double bi,
                   final double cr, final double ci,
                   final double dr, final double di) {
      super.assign( ar, ai, br, bi, cr, ci, dr, di );
    }

    public Moebius(final Field.Complex a, final Field.Complex b,
                   final Field.Complex c, final Field.Complex d) {
      super.assign( a.getRe(), a.getIm(), b.getRe(), b.getIm(), c.getRe(), c.getIm(), d.getRe(), d.getIm() );
    }

    public Moebius( final AbstractComplex2By2 m ) {
      super( m );
    }

    /** creates loxodromic transformation s with the two different
     ** fix points A=(ar,ai) and B=(br,bi), such that: <t>(s(z)-B)/(s(z)-A)=m(z-B)/(z-A)
     ** <p>
     ** m=(mr,mi) is the square of one of the eigenvalues of the element of SL2C which
     ** represents s. */
    public Moebius( final double ar, final double ai,
		    final double br, final double bi,
		    final double mr, final double mi) {
	assign(ar,ai,br,bi,mr,mi);
    }

    /** creates loxodromic transformation s with the two different
     ** fix points A and B, such that: <t>(s(z)-B)/(s(z)-A)=m(z-B)/(z-A)
     ** <p>
     ** m is the square of one of the eigenvalues of the element of SL2C which
     ** represents s. */
    public Moebius( final Field.Complex A, final Field.Complex B, final Field.Complex m ) {
	assign( A.getRe(), A.getIm(), B.getRe(), B.getIm(), m.getRe(), m.getIm() );
    }

    /**
     * Assign this with prescribed entries.
     */
    public void assign( final double ar, final double ai,
			final double br, final double bi,
			final double cr, final double ci,
			final double dr, final double di) {
	super.assign( ar, ai, br, bi, cr, ci, dr, di );
    }

    /**
     * Assign this with prescribed entries.
     */
    public void assign(final Field.Complex a, final Field.Complex b,
		       final Field.Complex c, final Field.Complex d) {
	super.assign( a.getRe(), a.getIm(), b.getRe(), b.getIm(), c.getRe(), c.getIm(), d.getRe(), d.getIm() );
    }

    /**
     * Assign this with the given transformation.
     */
    public void assign( final AbstractComplex2By2 m) {
        super.assign(m);
    }

    /**
     * Assign this with the given transformation.
     */
    public void assign( final Moebius m) {
        super.assign(m);
    }

    /** set instance to the loxodromic transformation s with the two different
     ** fix points A=(ar,ai) and B=(br,bi), such that: <t>(s(z)-B)/(s(z)-A)=m(z-B)/(z-A)
     ** <p>
     ** m=(mr,mi) is the square of one of the eigenvalues of the element of SL2C which
     ** represents s. */
    final public void assign( final double ar, final double ai,
			      final double br, final double bi,
			      final double mr, final double mi) {

	final double abRe = ar*br - ai*bi;
	final double abIm = ar*bi + ai*br;

	aRe = ( ar*mr - ai*mi ) - br;
	aIm = ( ar*mi + ai*mr ) - bi;

	bRe = abRe - ( abRe*mr - abIm*mi );
	bIm = abIm - ( abRe*mi + abIm*mr );

	cRe = mr - 1;
	cIm = mi;

	dRe = ar - ( br*mr - bi*mi );
	dIm = ai - ( br*mi + bi*mr );

	try {
	assignNormalizeDeterminant();
	}
	catch( ArithmeticException e ) {
	    System.out.println( this );
	    throw e;
	}
    }




    /** Assign this with the identity transformation */
    public void assignIdentity() {
	super.assignIdentity();
    }

    /** set instance to the loxodromic transformation s with the two different
     ** fix points A and B, such that: <t>(s(z)-B)/(s(z)-A)=m(z-B)/(z-A)
     ** <p>
     ** m is the square of one of the eigenvalues of the element of SL2C which
     ** represents s. */
    final public void assign( final Field.Complex a, final Field.Complex b, final Field.Complex m ) {
	assign( a.getRe(), a.getIm(), b.getRe(), b.getIm(), m.getRe(), m.getIm() );
    }

    /**
     * <p>Assign <code>this</code> by prescribing its fixed points and the
     * eigenvalues <code>lambda</code>, <code>1/lambda</code> of the
     * corresponding <code>2&times;2</code>-matrix.</p>
     *
     * <p>If <code>fixPoint1</code> and <code>fixPoint2</code> correspond to
     * the complex numbers <code>z1</code> and <code>z2</code>, then
     * <code>this</code> is set to the transformation <code>T</code> with
     *
     * <pre>
     *      T(z) - z2           2  z - z2
     *      ---------  =  lambda   ------  .
     *      T(z) - z1              z - z1
     * </pre></p>
     *
     * @param fixPoint1 a {@link ComplexProjective1}: the first fix point
     * @param fixPoint2 a {@link ComplexProjective1}: the second fix point
     * @param lambdaRe a <code>double</code>: real part of <code>lambda</code>
     * @param lambdaIm a <code>double</code>: imaginary part of
     * <code>lambda</code>
     * @throws IllegalArgumentException if <code>lambdaRe == 0.0 && lambdaIm
     * == 0.0</code> or if the fix points coincide.
     */
    public void assign(final ComplexProjective1 fixPoint1,
                       final ComplexProjective1 fixPoint2,
                       final double lambdaRe,
                       final double lambdaIm) throws IllegalArgumentException {
        if (lambdaRe == 0.0 && lambdaIm == 0.0) {
            throw new IllegalArgumentException("lambda must not be zero. ");
        }
        final double lambdaAbsSq = lambdaRe * lambdaRe + lambdaIm * lambdaIm;
        final double lambdaInvRe = lambdaRe / lambdaAbsSq;
        final double lambdaInvIm = -lambdaIm / lambdaAbsSq;
        super.assignByEigenvectors(lambdaRe, lambdaIm, fixPoint1,
                                   lambdaInvRe, lambdaInvIm, fixPoint2);
    }

    /**
     * <p>Assign <code>this</code> with a Euclidean scale-rotation.</p>
     *
     * <p>One fix point is <code>center</code>, the other is at infinity.</p>
     *
     * @param center a {@link ComplexProjective1}: center of the
     * scale-rotation
     * @param logScale a <code>double</code>: the logarithmic scale
     * @param angle a <code>double</code>: the angle of rotation
     */
    public void assignEuclideanLogScaleRotation(final
                                                ComplexProjective1 center,
                                                final double logScale,
                                                final double angle) {

        /* Let dummyComplexProjective1 be the point at infinity*/
        dummyComplexProjective1.aRe = 1.0;
        dummyComplexProjective1.aIm = 0.0;
        dummyComplexProjective1.bRe = 0.0;
        dummyComplexProjective1.bIm = 0.0;

        /* assign by fixed points and eigenvector */
        final double halfAngle = 0.5 * angle;
        final double scale = Math.exp(0.5 * logScale);
        assign(dummyComplexProjective1, center,
//         assign(center, dummyComplexProjective1,
               scale * Math.cos(halfAngle),
               scale * Math.sin(halfAngle));
    }

    /**
     * <p>Assign <code>this</code> with a sperical scale-rotation.</p>
     *
     * <p>One fix point is <code>center</code>, the other is the point
     * diametrically opposite on the Riemann sphere.</p>
     *
     * @param center a {@link ComplexProjective1}: center of the
     * scale-rotation
     * @param logScale a <code>double</code>: the logarithmic scale
     * @param angle a <code>double</code>: the angle of rotation
     */
    public void assignSphericalLogScaleRotation(final
                                                ComplexProjective1 center,
                                                final double logScale,
                                                final double angle) {

        /* Let dummyComplexProjective1 be the point diametrically opposite
         * center on the Riemann sphere. */
        dummyComplexProjective1.aRe = -center.bRe;
        dummyComplexProjective1.aIm = center.bIm;
        dummyComplexProjective1.bRe = center.aRe;
        dummyComplexProjective1.bIm = -center.aIm;

        /* assign by fixed points and eigenvector */
        final double halfAngle = 0.5 * angle;
        final double scale = Math.exp(0.5 * logScale);
        assign(dummyComplexProjective1, center,
               scale * Math.cos(halfAngle),
               scale * Math.sin(halfAngle));
    }

    /**
     * Set this to moebius transform which maps a onto A, b onto B, and c
     * onto C.  this methods creates manny temporarily used objects and
     * therefore it needs to be optimized.
     */
    final public void assign( final Complex a, final Complex b, final Complex c,
			   final Complex A, final Complex B, final Complex C ) {

	Complex Alpha = B.minus(C); Alpha.assignDivide( B.minus( A ) );
	Complex Beta  = A.times( Alpha ); Beta.assignNeg();
	ComplexConstant Gamma = ComplexConstant.ONE;
	Complex Delta = C.times( -1 );

	Moebius M = new Moebius( Alpha, Beta, Gamma, Delta );

	Complex alpha = b.minus(c); alpha.assignDivide( b.minus( a ) );
	Complex beta  = a.times( alpha ); beta.assignNeg();
	ComplexConstant gamma = ComplexConstant.ONE;
	Complex delta = c.times( -1 );

	Moebius m = new Moebius( alpha, beta, gamma, delta );

	assignInvert( M );
	assignTimes( m );
    }

     /**
      * assigns this with the Moebius transform wich maps (0,1) to a, (1,1) to b and
      * (1,0) (=infinity) to c
      * uses tons of new Field.Complex numbers - optimize!!
      */
     final public void assign( final ComplexProjective1 a, final ComplexProjective1 b, final ComplexProjective1 c) {

         Complex a1 = a.getA();
         Complex a2 = a.getB();
         Complex b1 = b.getA();
         Complex b2 = b.getB();
         Complex c1 = c.getA();
         Complex c2 = c.getB();

         Complex num   = new Complex();
         Complex tmp   = new Complex();
         Complex denom = new Complex();

         num.assignTimes(c1, a2);
         tmp.assignTimes(a1, c2);
         num.assignMinus(tmp);

         denom.assignTimes(c1, b2);
         tmp.assignTimes(b1, c2);
         denom.assignMinus(tmp);

         denom.assignDivide(num);

         tmp.assignTimes(denom, a1);
         setB(tmp);

         tmp.assignTimes(denom, a2);
         setD(tmp);

         denom.assignTimes(b1, a2);
         tmp.assignTimes(b2, a1);
         denom.assignMinus(tmp);

         denom.assignDivide(num);

         tmp.assignTimes(denom, c1);
         setA(tmp);

         tmp.assignTimes(denom, c2);
         setC(tmp);

	 super.assignNormalizeDeterminant();
     }

    /**
     * Assign this with the inverse of <code>m</code>.
     */
    public void assignInvert(final Moebius m) {
        super.assignInvert(m);
    }

    /**
     * Assign this with its own inverse.
     */
    public void assignInvert() {
	super.assignInvert(this);
    }


    /**
     * Returns inverse of this.
     */
    public Moebius invert() {
	Moebius result = new Moebius();
	result.assignInvert( this );
	return result;
    }

    /**
     * Assign this with its own adjoined.
     */
    public void assignAdjoined() {
	super.assignAdjoined();
    }

    /**
     * Returns adjoined of this.
     */
    public Moebius adjoined() {
	Moebius result = new Moebius();
	result.assignAdjoined( this );
	return result;
    }

    /** Assign this with the product of itself and a. */
    public void assignTimes(final Moebius a) {
	super.assignTimes(this, a);
    }

    /**
     * Multiplies <code>a</code> and <code>b</code> and store
     * the result in <code>c</code>.
     */
    public static void times(final Moebius a,
                             final Moebius b,
                             final Moebius c) {
        AbstractComplex2By2.times(a, b, c);
    }

    /**
     * Assign this with the product of <code>a</code> and <code>b</code>.
     */
    public void assignTimes(final Moebius a, final Moebius b) {
	super.times( a, b, this );
    }

    /**
     * Returns product of this with <code>a</code>.
     */
    public Moebius times( final AbstractComplex2By2 a ) {
	Moebius result = new Moebius();
	super.times( this, a, result );
	return result;
    }

    /**
     * Assign this with the product of itself and the inverse of
     * <code>m</code>.
     */
    public void assignDivide(final Moebius m) {
	super.assignDivide(this, m);
    }

    /**
     * Assign this with the product of itself and the inverse of
     * <code>c</code>.
     */
    public void assignDivide(final Complex c) {
	super.assignDivide(this, c);
    }

    /**
     * Assign this with the product of itself and the inverse of
     * <code>c</code>.
     */
    public void assignDivide(final double c) {
        super.assignDivide(this, c);
    }


    /**
     * Assign this with the product of <code>a</code> and the inverse of
     * <code>b</code>.
     */
    public void assignDivide(final Moebius a, final Moebius b) {
        super.assignDivide(a, b);
    }


    /**
     * Returns ratio of this and <code>a</code>.
     */
    public Moebius divide( final AbstractComplex2By2 a ) {
	Moebius result = new Moebius();
	super.divide( this, a, result );
	return result;
    }

    /**
     * Apply <code>this</code> to a {@link ComplexProjective1}.
     *
     * @param z a {@link ComplexProjective1}: point to which
     * <code>this</code> is applied. Value is not changed (unless <code>v ==
     * w</code>)
     * @param w a {@link ComplexProjective1}: used to hand over the result
     */
    public void applyTo(final ComplexProjective1 v,
                        final ComplexProjective1 w) {
        AbstractComplex2By2.times(this, v, w);
    }

    /**
     * Apply <code>this</code> to a {@link ComplexProjective1}.
     *
     * @param z a {@link ComplexProjective1}: point to which
     * <code>this</code> is applied. Value is changed to the result.
     */
    public void applyTo(final ComplexProjective1 v) {
        applyTo(v, v);
    }

    /** applies instance to z and assigns r with the result */
    final public void applyTo( final Complex z, final Complex r ) {
	final double d1r = aRe*z.re - aIm*z.im + bRe;
	final double d1i = aRe*z.im + aIm*z.re + bIm;

	final double d2r = cRe*z.re - cIm*z.im + dRe;
	final double d2i = cRe*z.im + cIm*z.re + dIm;

	final double d = d2r*d2r + d2i*d2i;

	r.re = ( d1r*d2r + d1i*d2i ) / d;
	r.im = ( d1i*d2r - d1r*d2i ) / d;
    }


    /** applies instance to z and assigns r with the result */
    final public void applyTo( final Field.Complex z, final Complex r ) {
        Complex w = new Complex(z);
        applyTo(w,r);
    }

    /** applies instance to z and assigns r with the result */
    final public void eval( final Complex z, final Complex r ) {
    	applyTo(z,r);
    }
    
    /** applies instance to z and returns the result */
    final public Complex applyTo( final Complex z ) {
	final Complex r = new Complex();
	applyTo( z, r );
	return r;
    }

    final public Complex valueAt( final Complex z ) {
        final Complex r = new Complex();
        applyTo( z, r );
        return r;
    }


        /** applies inverse of instance to z and assigns r with the result */
        final public void applyInverseTo( final Complex z, final Complex r ) {


            final double d1r =  dRe*z.re - dIm*z.im - bRe;
            final double d1i =  dRe*z.im + dIm*z.re - bIm;

            final double d2r = -cRe*z.re + cIm*z.im + aRe;
            final double d2i = -cRe*z.im - cIm*z.re + aIm;

            final double d = d2r*d2r + d2i*d2i;

            r.re = ( d1r*d2r + d1i*d2i ) / d;
            r.im = ( d1i*d2r - d1r*d2i ) / d;
        }

        /** applies inverse of instance to z and assigns r with the result */
        final public void alt_applyInverseTo( final Complex z, final Complex r ) {
          final double d1r =  dRe*z.getRe() - dIm*z.getIm() - bRe;
            final double d1i =  dRe*z.getIm() + dIm*z.getRe() - bIm;

            final double d2r = -cRe*z.getRe() + cIm*z.getIm() + aRe;
            final double d2i = -cRe*z.getIm() - cIm*z.getRe() + aIm;

            final double d = d2r*d2r + d2i*d2i;

            r.re = ( d1r*d2r + d1i*d2i ) / d;
            r.im = ( d1i*d2r - d1r*d2i ) / d;
        }

    /** applies inverse of instance to z and returns the result */
    final public Complex applyInverseTo( final Complex z ) {
	final Complex r = new Complex();
	applyInverseTo( z, r );
	return r;
    }

    /** evaluates differential of instance at z and assigns r with the result */
    final public void applyDifferentialTo( final Complex z, final Complex r ) {
	final double dr = cRe*z.re - cIm*z.im + dRe;
	final double di = cRe*z.im + cIm*z.re + dIm;

	double d = dr*dr + di*di; d*=d;

	r.re = ( dr*dr - di*di ) / d;
	r.im = -2 * dr*di / d;
    }

    /** evaluates differential of instance at z and returns the result */
    final public Complex applyDifferentialTo( final Complex z ) {
	final Complex r = new Complex();
	applyDifferentialTo( z, r );
	return r;
    }

    /** evaluates differential of inverse of instance at z and assigns r with the result */
    final public void applyInverseDifferentialTo( final Complex z, final Complex r ) {
	double dr = -cRe*z.re + cIm*z.im + aRe;
	double di = -cRe*z.im - cIm*z.re + aIm;

	double d = dr*dr + di*di; d*=d;

	r.re = ( dr*dr - di*di ) / d;
	r.im = -2*dr*di / d;
    }

    /** evaluates differential of inverse of instance at z and returns the result */
    final public Complex applyInverseDifferentialTo( final Complex z ) {
	final Complex r = new Complex();
	applyInverseDifferentialTo( z, r );
	return r;
    }

    final public boolean isElliptic() {

	double detRe = aRe*dRe - aIm*dIm  - ( bRe*cRe - bIm*cIm );
	double detIm = aRe*dIm + aIm*dRe  - ( bRe*cIm + bIm*cRe );

	double nn = detRe * detRe + detIm * detIm;

	double normTraceRe = ( (aRe+dRe)*detRe + (aIm+dIm)*detIm ) / nn;
	double normTraceIm = ( (aRe+dRe)*detIm - (aIm+dIm)*detRe ) / nn;

	if( Math.abs( normTraceIm ) > EPS)
	    return false;

	return Math.abs( normTraceRe ) < 2 - EPS;
    }

    final public boolean isParabolic() {

	double detRe = aRe*dRe - aIm*dIm  - ( bRe*cRe - bIm*cIm );
	double detIm = aRe*dIm + aIm*dRe  - ( bRe*cIm + bIm*cRe );

	double nn = detRe * detRe + detIm * detIm;

	double normTraceRe = ( (aRe+dRe)*detRe + (aIm+dIm)*detIm ) / nn;
	double normTraceIm = ( (aRe+dRe)*detIm - (aIm+dIm)*detRe ) / nn;

	if( Math.abs( normTraceIm ) > EPS)
	    return false;

	return Math.abs( Math.abs( normTraceRe ) - 2 ) < EPS;
    }

    final public boolean isHyperbolic() {

	double detRe = aRe*dRe - aIm*dIm  - ( bRe*cRe - bIm*cIm );
	double detIm = aRe*dIm + aIm*dRe  - ( bRe*cIm + bIm*cRe );

	double nn = detRe * detRe + detIm * detIm;

	double normTraceRe = ( (aRe+dRe)*detRe + (aIm+dIm)*detIm ) / nn;
	double normTraceIm = ( (aRe+dRe)*detIm - (aIm+dIm)*detRe ) / nn;

	if( Math.abs( normTraceIm ) > EPS)
	    return false;

	return Math.abs( normTraceRe ) > 2 + EPS;
    }

    /** we call a moebius loxodromic if its norm is not 1, this is not standart. */
    final public boolean isLoxodromic() {
	return Math.abs( normSqr() - 1 ) > EPS;
    }

    /** assigns array fp with fix points of instance and returns the number of fixpoints,
     ** which can be either 1, 2, or infinity. In the last case 0 is returned.
     ** if a needed entry of the array fp is null an instance of Field.Complex is created. */
    final public int getFixPoints( final Complex[] fp ) {

	final double dMinusARe = dRe - aRe;
	final double dMinusAIm = dIm - aIm;

	if( Math.abs(cRe) + Math.abs(cIm) < 2*EPS ) {                     // c=0, affine case

	    if( Math.abs( dMinusARe ) + Math.abs( dMinusAIm ) < 2*EPS ) { // a=d

		if( Math.abs(bRe) + Math.abs(bIm) < 2*EPS ) {             // b=0, identity
		    return 0;
		} else {                                                      // b!=0
		    if( fp[0] == null )
			fp[0] = new Complex();
		    fp[0].assign( Double.POSITIVE_INFINITY, 0 );
		    return 1;
		}
	    } else {                                                          // a!=d
		if( fp[0] == null )
		    fp[0] = new Complex();

		if( Math.abs(bRe) + Math.abs(bIm) < 2*EPS ) {             // b=0, linear case
		    fp[0].assign( 0 );
		    if( fp[1] == null )
			fp[1] = new Complex();
		    fp[1].assign( Double.POSITIVE_INFINITY, 0 );
		    return 1;
		} else {                                                      // b!=0
		    fp[0].assign( bRe, bIm );
		    fp[0].assignDivide( dMinusARe, dMinusAIm );
		    return 1;
		}
	    }
	} else {
	    if( fp[0] == null )
		fp[0] = new Complex();

	    // dis = (d-a)^2 + 4cb
	    final double disRe = dMinusARe*dMinusARe - dMinusAIm*dMinusAIm + 4 * ( cRe*bRe - cIm*bIm );
	    final double disIm = dMinusARe*dMinusAIm + dMinusAIm*dMinusARe + 4 * ( cRe*bIm + cIm*bRe );

	    if( Math.abs( disRe ) + Math.abs( disIm ) < 2*EPS ) {        // dis=0
		fp[0].assign( -dMinusARe, -dMinusAIm );
		fp[0].assignDivide( 2*cRe, 2*cIm );
		return 1;
	    } else {                                                         // dis!=0
		if( fp[1] == null )
		    fp[1] = new Complex();

		final double rr = Math.sqrt( Complex.abs( disRe, disIm ) );
		final double ii = Complex.arg( disRe, disIm ) / 2;

		final double sqrtOfDisRe = rr * Math.cos( ii );
		final double sqrtOfDisIm = rr * Math.sin( ii );

		fp[0].assign( -dMinusARe, -dMinusAIm );
		fp[1].assign( -dMinusARe, -dMinusAIm );
		fp[0].assignMinus( sqrtOfDisRe, sqrtOfDisIm );
		fp[1].assignPlus ( sqrtOfDisRe, sqrtOfDisIm );
		fp[0].assignDivide( 2*cRe, 2*cIm );
		fp[1].assignDivide( 2*cRe, 2*cIm );
		return 2;
	    }
	}
    }

    /** returns array of fixpoints of instance, if its number is either 1 or 2.
     ** if instance is the identity the method returns null */
    final public Complex[] getFixPoints() {
	Complex [] tmpFP = new Complex[2];

	int numOfFixPoints = getFixPoints( tmpFP );

	if( numOfFixPoints == 0 )
	    return null;

	Complex [] fp = new Complex[ numOfFixPoints ];

	for( int i=0; i<numOfFixPoints; i++ )
	    fp[i] = tmpFP[i];

	return fp;
    }

    /** computes the image of a circle given by center and radius.
     ** @param centerOfCircle center of circle which is to map
     ** @param radiusOfCircle radius of circle which is to map
     ** @param centerOfMappedCircle center of mapped circle
     ** @return radius of mapped circle */
    final public double getRadiusOfMappedCircle( final Complex centerOfCircle,
						 final  double radiusOfCircle,
						 final Complex centerOfMappedCircle ) {
	double radiusOfMappedCircle;

	if( cRe == 0 && cIm == 0 ) {

	    // sigma(z) = a/d z + b/d

	    double dAbsSqr =   dRe*dRe + dIm*dIm;

	    double tmpARe  = ( aRe*dRe + aIm*dIm ) / dAbsSqr;
	    double tmpAIm  = (-aRe*dIm + aIm*dRe ) / dAbsSqr;

	    double tmpBRe  = ( bRe*dRe + bIm*dIm ) / dAbsSqr;
	    double tmpBIm  = (-bRe*dIm + bIm*dRe ) / dAbsSqr;

	    centerOfMappedCircle.re = centerOfCircle.re * tmpARe - centerOfCircle.im * tmpAIm + tmpBRe;
	    centerOfMappedCircle.im = centerOfCircle.re * tmpAIm + centerOfCircle.im * tmpARe + tmpBIm;

	    radiusOfMappedCircle = radiusOfCircle * Math.sqrt( tmpARe*tmpARe + tmpAIm*tmpAIm );

	} else {

	    // sigma(z) = h g f(z)  with f(z) = cz + d  , g(z) = 1/z , and h(z) = -z/c + a/c

	    centerOfMappedCircle.re = centerOfCircle.re * cRe - centerOfCircle.im * cIm + dRe;
	    centerOfMappedCircle.im = centerOfCircle.re * cIm + centerOfCircle.im * cRe + dIm;

	    double cAbsSqr =  cRe*cRe + cIm*cIm;
	    double cAbs    =  Math.sqrt( cAbsSqr );

	    radiusOfMappedCircle = radiusOfCircle * cAbs;

	    double factor = centerOfMappedCircle.absSqr() - radiusOfMappedCircle*radiusOfMappedCircle;

	    centerOfMappedCircle.assignDivide( factor );
            centerOfMappedCircle.assignConjugate();

	    radiusOfMappedCircle /= Math.abs( factor );

	    double tmpARe  = ( aRe*cRe + aIm*cIm ) / cAbsSqr;
	    double tmpAIm  = (-aRe*cIm + aIm*cRe ) / cAbsSqr;

	    double tmpRe  = ( centerOfMappedCircle.re*cRe + centerOfMappedCircle.im*cIm ) / cAbsSqr;
	    double tmpIm  = (-centerOfMappedCircle.re*cIm + centerOfMappedCircle.im*cRe ) / cAbsSqr;

	    centerOfMappedCircle.re = tmpARe - tmpRe;
	    centerOfMappedCircle.im = tmpAIm - tmpIm;

	    radiusOfMappedCircle /= cAbs;
	}

	return radiusOfMappedCircle;
    }
}


