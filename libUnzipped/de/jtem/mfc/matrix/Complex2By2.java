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
import de.jtem.mfc.vector.Complex2;

/**
 ** This class represents a complex 2-by-2 matrix.
 ** <pre>
 ** [  a   b  ]     [  aRe + i aIm      bRe + i bIm  ] 
 ** [         ]  =  [                                ]
 ** [  c   d  ]     [  dRe + i dIm      dRe + i dIm  ] 
 ** </pre>
 ** <p>
 ** To minimize the object count its 4 complex entries <code>a,b,c,d</code>
 ** are represented by 8 double values 
 ** <code>aRe, aIm, bRe, bIm, cRe, cIm, dRe, dIm</code>.
 ** All methods avoid the creation of temporarily used objects in their
 ** internal computations unless stated in their describtions.
 ** <p>
 ** Creation of many temporarily used instances is expensive and stresses 
 ** the garabage collector.
 ** This can cause your java program to be slow and should therefore be avoided.
 ** Thus, operations which result in an instance of <code>Complex2By2</code> 
 ** can be performed either with or without creating a matrix.
 ** To compute, for example, the product of two matrices you can either use
 ** <pre>
 ** Complex2By2 m = a.times( b ),
 ** </pre>
 ** or
 ** <pre>
 ** m.assignTimes( a, b ).
 ** </pre>
 ** This philosophy is applied to operations which result in other instances,
 ** as well. In such a case an instance of the resulting type
 ** can be prescribed as a parameter which is than filled by
 ** the method. To compute the eigenvalues of a matrix, for example,
 ** <code>m</code> you can use either
 ** <pre>
 ** double[] ev = m.getEigenValues()
 ** </pre>
 ** or
 ** <pre> m.getEigenValues( ev )
 ** </pre>
 ** <p>
 ** @author schmies */

public class Complex2By2 extends AbstractComplex2By2 implements Serializable, Cloneable {
    
    private static final long serialVersionUID = 1L;

    /** creates an identity matrix */
    public Complex2By2() { 	
    }

    /** creates a matrix with the prescribed entries */
    public Complex2By2( final double aRe, final double aIm, final double bRe, final double bIm,
			final double cRe, final double cIm, final double dRe, final double dIm) { 
	super( aRe, aIm, bRe, bIm, cRe, cIm, dRe, dIm );
    }

    /** creates a matrix with the prescribed entries */
    public Complex2By2( final Complex a, final Complex b, final Complex c, final Complex d) { 
	super( a, b, c, d );
    }

    /** creates a matrix equal to the prescribed one */
    public Complex2By2( final AbstractComplex2By2 m ) { 
	super( m );
    }

    /** sets entry <code>a</code> of this with <code>a</code> */
    final public void setA( final Complex a ) {
        super.setA(a);
// 	aRe = a.re;
// 	aIm = a.im;
    }

    /** sets entry <code>b</code> of this with <code>b</code> */
    final public void setB( final Complex b ) {
        super.setB(b);
// 	bRe = b.re;
// 	bIm = b.im;
    }

    /** sets entry <code>c</code> of this with <code>c</code> */
    final public void setC( final Complex c ) {
        super.setC(c);
// 	cRe = c.re;
// 	cIm = c.im;
    }
    
    /** sets entry <code>d</code> of this with <code>d</code> */
    final public void setD( final Complex d ) {
        super.setD(d);
// 	dRe = d.re;
// 	dIm = d.im;
    }

    /** returns copy of this */
    final public Complex2By2 copy() { 
	return new Complex2By2( this );
    }
    
    /** assigns this with the zero matrix */
    final public void assignZero() {
	super.assignZero();
    }

    /** assigns this with the identity matrix */
    final public void assignIdentity() {
	super.assignIdentity();
    }

    /** assigns this with the prescribed entries */
    final public void assign( final double aRe, final double aIm, 
			      final double bRe, final double bIm,
			      final double cRe, final double cIm, 
			      final double dRe, final double dIm) {
	super.assign(  aRe,  aIm,  bRe,  bIm, cRe,  cIm,  dRe,  dIm); 
    }

    /** assigns this with the prescribed entries */
    final public void assign( final Complex a, final Complex b, 
                              final Complex c, final Complex d ) {
        super.assign(a, b, c, d);
//         assign( a.re, a.im, b.re, b.im, c.re, c.im, d.re, d.im );
    }
    
    /** assigns this with the prescribed matrix */
    final public void assign( final AbstractComplex2By2 s ) { 
        super.assign(s);
// 	assign( s.aRe, s.aIm, s.bRe, s.bIm, s.cRe, s.cIm, s.dRe, s.dIm ); 
    }

    /** assigns this with the product of itself and a */
    final public void assignTimes( final AbstractComplex2By2 a ) { 
	super.assignTimes(this, a ); 
    }

    /** assigns this with the product of <code>a</code> and <code>b</code> */
    final public void assignTimes( final AbstractComplex2By2 a, final AbstractComplex2By2 b ) {
	super.assignTimes( a, b );
    }

    /** returns the product of this and <code>a</code> */
    final public Complex2By2 times( final AbstractComplex2By2 a ) { 
	Complex2By2 result = new Complex2By2();
	result.assignTimes( this, a );
	return result;
    }

    /** assigns this with the product of itself and <code>r</code> */
    final public void assignTimes( final double r ) {
	assignTimes( this, r );
    }

    /** assigns this with product of <code>m</code> and <code>r</code> */
    final public void assignTimes( final AbstractComplex2By2 m, final double r ) {
	super.assignTimes( m, r );
    }

    /** returns the product of this and <code>r</code> */
    final public Complex2By2 times( final double r ) { 
	Complex2By2 result = new Complex2By2();
	result.assignTimes( this, r );
	return result;
    }

    /** assigns this with the product of itself and <code>x+iy</code> */
    final public void assignTimes( final double x, final double y ) {
	assignTimes( this, x, y );
    }

    /** assigns this with the product of <code>m</code> and <code>x+iy</code> */
    final public void assignTimes( final AbstractComplex2By2 m, final double x, final double y ) {
	super.assignTimes( m, x, y );
    }

    /** returns product of this and <code>x+iy</code> */
    final public Complex2By2 times( final double x, final double y ) { 
	Complex2By2 result = new Complex2By2();
	result.assignTimes( this, x, y );
	return result;
    }

    /** returns the product of this and <code>z</code> */
    final public void assignTimes( final Complex z ) {
	assignTimes( this, z.re, z.im );
    }

    /** assigns this with product of <code>m</code> and <code>z</code> */
    final public void assignTimes( final AbstractComplex2By2 m, final Complex z ) {
	assignTimes( m, z.re, z.im );
    }

    /** returns the product of this and <code>z</code> */
    final public Complex2By2 times( final Complex z ) { 
	Complex2By2 result = new Complex2By2();
	result.assignTimes( this, z.re, z.im );
	return result;
    }

    /** assigns this with the product of itself and inverse of <code>m</code>*/
    final public void assignDivide( final AbstractComplex2By2 m ) { 
	assignDivide( this, m ); 
    }

    /** assigns this with the product of <code>a</code> and inverse of <code>b</code> */
    final public void assignDivide( final AbstractComplex2By2 a, final AbstractComplex2By2 b ) {
	super.assignDivide( a, b );
    }

    /** returns the product of this and inverse of <code>a</code> */
    final public Complex2By2 divide( final AbstractComplex2By2 a ) { 
	Complex2By2 result = new Complex2By2();
	result.assignDivide( this,  a );
	return result;
    }

    /** assigns this with product of itself and <code>1/r</code> */
    final public void assignDivide( final double r ) {
	assignDivide( this, r );
    }

    /** assigns this with product of <code>m</code> and <code>1/r</code> */
    final public void assignDivide( final AbstractComplex2By2 m, final double r ) {
	super.assignDivide( m, r );
    }

    /** returns the product of this and <code>1/r</code> */
    final public Complex2By2 divide( final double r ) {
	Complex2By2 result = new Complex2By2();
	result.assignDivide( this, r );
	return result;
    }

    /** assigns this with the product of itself and <code>1/(x+iy)</code> */
    final public void assignDivide( final double x, final double y ) {
	assignDivide( this, x, y );
    }

    /** assigns this with the product of <code>m</code> and <code>1/(x+iy)</code> */
    final public void assignDivide( final AbstractComplex2By2 m, final double x, final double y ) {
	super.assignDivide( m, x, y );
    }
   
    /** returns the product of this and <code>1/(x+iy)</code> */
    final public Complex2By2 divide( final double x, final double y ) { 
	Complex2By2 result = new Complex2By2();
	result.assignDivide( this, x, y );
	return result;
    }


    /** assigns this with the product of itself and <code>1/z</code> */
    final public void assignDivide( final Complex z ) {
	assignDivide( this, z.re, z.im );
    }

    /** assigns this with the product of <code>m</code> and <code>1/z</code> */
    final public void assignDivide( final AbstractComplex2By2 m, final Complex z ) {
	super.assignDivide( m, z );
    }
   
    /** returns the product of this and <code>1/z</code> */
    final public Complex2By2 divide( final Complex z ) { 
	Complex2By2 result = new Complex2By2();
	result.assignDivide( this, z );
	return result;
    }
    
    /** returns this with the sum of this and a */
    final public Complex2By2 plus( final AbstractComplex2By2 a ) { 
	Complex2By2 result = new Complex2By2();
	result.assignPlus( this, a );
	return result;
    }

    /** assigns this with the sum of itself and a */
    final public void assignPlus( final AbstractComplex2By2 a ) { 
	assignPlus(this, a ); 
    }

    /** assigns this with the sum of <code>a</code> and <code>b</code> */
    final public void assignPlus( final AbstractComplex2By2 a, final AbstractComplex2By2 b ) {
	super.assignPlus( a, b );
    }

    /** adds <code>x</code> to all entries of this and returns the result */
    final public Complex2By2 plus( final double x ) { 
	Complex2By2 result = new Complex2By2();
	result.assignPlus( this, x );
	return result;
    }

    /** adds <code>x</code> to all entries of this */
    final public void assignPlus( final double x ) {
	super.assignPlus( x );
    }

    /** adds <code>x</code> to all entries of <code>m</code> and 
	and assigns this with the result. */
    final public void assignPlus( final AbstractComplex2By2 m, final double x ) {
	super.assignPlus( m, x );
    }

    /** adds <code>x+iy</code> to all entries of this and returns the result */
    final public Complex2By2 plus( final double x, final double y ) { 
	Complex2By2 result = new Complex2By2();
	result.assignPlus( this, x, y );
	return result;
    }

    /** adds <code>x+iy</code> to all entries of this */
    final public void assignPlus( final double x, final double y ) {
	super.assignPlus( x, y );
    }

    /** adds <code>x+iy</code> to all entries of <code>m</code> and 
	and assigns this with the result. */
    final public void assignPlus( final AbstractComplex2By2 m, final double x, final double y ) {
	super.assignPlus( m, x, y );
    }

    /** adds <code>z</code> to all entries of this and returns the result */
    final public Complex2By2 plus( final Complex z ) { 
	Complex2By2 result = new Complex2By2();
	result.assignPlus( this, z );
	return result;
    }

    /** adds <code>z</code> to all entries of this */
    final public void assignPlus( final Complex z ) {
	super.assignPlus( z );
    }

    /** adds <code>z</code> to all entries of <code>m</code> and 
	and assigns this with the result. */
    final public void assignPlus( final AbstractComplex2By2 m, final Complex z ) {
	super.assignPlus( m, z );
    }

    /** subtracts <code>a</code> from this and returns the result */
    final public Complex2By2 minus( final AbstractComplex2By2 a ) { 
	Complex2By2 result = new Complex2By2();
	result.assignMinus( this, a );
	return result;
    }

    /** subtracts <code>a</code> from this */
    final public void assignMinus( final AbstractComplex2By2 a ) { 
	super.assignMinus( a );
    }

    /** subtracts <code>b</code> from <code>a</code> and stores the result it this */
    final public void assignMinus( final AbstractComplex2By2 a, final AbstractComplex2By2 b ) {
	super.assignMinus( a, b );
    }

    /** subtracts <code>x</code> from all entries of this and returns the result */
    final public Complex2By2 minus( final double x ) { 
	Complex2By2 result = new Complex2By2();
	result.assignMinus( this, x );
	return result;
    }

    /** subtracts <code>x</code> from all entries of this */
    final public void assignMinus( final double x ) {
	super.assignMinus( x );
    }

    /** subtracts <code>x</code> from all entries of <code>m</code> and 
	and assigns this with the result. */
    final public void assignMinus( final AbstractComplex2By2 m, final double x ) {
	super.assignMinus( m, x );
    }

    /** subtracts <code>x+iy</code> from all entries of this and returns the result */
    final public Complex2By2 minus( final double x, final double y ) { 
	Complex2By2 result = new Complex2By2();
	result.assignMinus( this, x, y );
	return result;
    }

    /** subtracts <code>x+iy</code> to all entries of this */
    final public void assignMinus( final double x, final double y ) {
	super.assignMinus( x, y );
    }

    /** subtracts <code>x+iy</code> to all entries of <code>m</code> and 
	and assigns this with the result. */
    final public void assignMinus( final AbstractComplex2By2 m, final double x, final double y ) {
	super.assignMinus( m, x, y );
    }

    /** subtracts <code>z</code> from all entries of this and returns the result */
    final public Complex2By2 minus( final Complex z ) { 
	Complex2By2 result = new Complex2By2();
	result.assignMinus( this, z );
	return result;
    }

    /** subtracts <code>z</code> from all entries of this */
    final public void assignMinus( final Complex z ) {
	super.assignMinus( z );
    }

    /** subtracts <code>z</code> to all entries of <code>m</code> and 
	and assigns this with the result. */
    final public void assignMinus( final AbstractComplex2By2 m, final Complex z ) {
	super.assignMinus( m, z );
    }


    /** returns negative of this */
    final public Complex2By2 neg() { 
	Complex2By2 result = new Complex2By2();
	result.assignNeg( this );
	return result;
    }

    /** assigns this with the negative of <code>a</code> */
    final public void assignNeg( final AbstractComplex2By2 a ) {
	super.assignNeg( a );
    }

    /** assigns this with its own negative */
    final public void assignNeg() {
	super.assignNeg( this );
    }
   
    /** returns conjugate of this */
    final public Complex2By2 conjugate() { 
	Complex2By2 result = new Complex2By2();
	result.assignConjugate( this );
	return result;
    }

    /** assigns this with the conjugateative of <code>a</code> */
    final public void assignConjugate( final AbstractComplex2By2 a ) {
	super.assignConjugate();
    }

    /** assigns this with its own conjugate */
    final public void assignConjugate() {
	super.assignConjugate( this );
    }
         
    /** returns transposed of this */
    final public Complex2By2 transpose() { 
	Complex2By2 result = new Complex2By2();
	result.assignTranspose( this );
	return result;
    }

    /** assigns this with the transposed of <code>a</code> */
    final public void assignTranspose( final AbstractComplex2By2 a ) {
	super.assignTranspose();
    }

    /** assigns this with its own transposed */
    final public void assignTranspose() {
	super.assignTranspose( this );
    }

         
    /** returns transposed and conjugated of this */
    final public Complex2By2 star() { 
	Complex2By2 result = new Complex2By2();
	result.assignStar( this );
	return result;
    }

    /** assigns this with the transposed and conjugated of <code>a</code> */
    final public void assignStar( final AbstractComplex2By2 a ) {
	super.assignStar( a );
    }

    /** assigns this with its own transposed and conjugated */
    final public void assignStar() {
	super.assignStar( this );
    }
      
    /** assigns this with the adjugate of <code>a</code> */
    final public void assignAdjugate( final AbstractComplex2By2 a ) {
	super.assignAdjugate( a );
    }

    /** assigns this with its own adjugate */
    final public void assignAdjugate() {
	assignAdjugate( this );
    }

    /** returns the adjugate of this */
    final public Complex2By2 adjugate() { 
	Complex2By2 result = new Complex2By2();
	result.assignAdjugate( this );
	return result;
    }
    /** assigns this with the adjoined of <code>a</code> */
    final public void assignAdjoined( final AbstractComplex2By2 a ) {
	super.assignAdjoined( a );
    }

    /** assigns this with its own adjoined */
    final public void assignAdjoined() {
	assignAdjoined( this );
    }

    /** returns the adjoined of this */
    final public Complex2By2 adjoined() { 
	Complex2By2 result = new Complex2By2();
	result.assignAdjoined( this );
	return result;
    }

    /**
     * Assign <code>this</code> with <code>t this t<sup>*</sup>.
     * 
     * @param t an {@link AbstractComplex2By2}: the matrix to adjoin
     */
    public void assignAdjoinedWith(final AbstractComplex2By2 t) {
        super.assignAdjoinedWith(t);
    }

    /** assigns this with the inverse of <code>m</code> */
    final public void assignInvert( final AbstractComplex2By2 m ) {
	super.assignInvert( m );
    }

    /** assigns this with its own inverse */
    final public void assignInvert() {
	assignInvert( this );
    }

    /** returns inverse of this */
    final public Complex2By2 invert() {
	Complex2By2 result = new Complex2By2();
	result.assignInvert( this );
	return result;
    }

    /** returns a normalized (determiant equals 1) copy of this */
    final public Complex2By2 normalizeDeterminant() {
	Complex2By2 result = new Complex2By2();
	result.assignNormalizeDeterminant( this );
	return result;
    }

    /** assigns this with <code>m</code> and normalizes (determinant equals 1) it */
    final public void assignNormalizeDeterminant( final AbstractComplex2By2 m ) {
	super.assignNormalizeDeterminant( m );
    }

    /** normalizes (determinant equals 1) this */
    final public void assignNormalizeDeterminant() {
	assignNormalizeDeterminant( this );
    }

    /**
     * <p>Assign <code>this</code> by column.</p>
     *
     * <p>Simply calls {@link AbstractComplex2By2#assignByColumn(Complex2,
     * Complex2)}.</p>
     *
     * @param column1 a {@link Complex2}: the first column. Must not be
     * <code>null</code>
     * @param column2 a {@link Complex2}: the second column. Must not be
     * <code>null</code>
     */
    final public void assignByColumn(final Complex2 column1, 
                                     final Complex2 column2) {
        super.assignByColumn(column1, column2);
    }

    /**
     * <p>Assign <code>this</code> by prescribing eigenvectors and
     * eigenvalues. </p>
     *
     * <p> The given eigenvectors must be linearly independent. </p>
     *
     * <p>Simply calls {@link
     * AbstractComplex2By2#assignByEigenvectors(double, double, Complex2,
     * double, double, Complex2)}.</p>
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
    final public void assignByEigenvectors(final double eigenvalue1Re,
                                           final double eigenvalue1Im,
                                           final Complex2 eigenvector1,
                                           final double eigenvalue2Re,
                                           final double eigenvalue2Im,
                                           final Complex2 eigenvector2) {
        super.assignByEigenvectors(eigenvalue1Re, eigenvalue1Im, eigenvector1,
                                   eigenvalue2Re, eigenvalue2Im, eigenvector2);
    }

    /**
     * <p>Assign <code>this</code> by prescribing eigenvectors and
     * eigenvalues. </p>
     *
     * <p> The given eigenvectors must be linearly independent. </p>
     *
     * <p>Simply calls {@link
     * AbstractComplex2By2#assignByEigenvectors(Field.Complex, Complex2, Field.Complex,
     * Complex2)}.</p>
     *
     * @param eigenvalue1 a {@link Field.Complex}: first eigenvalue
     * @param eigenvector1 a {@link Complex2}: first eigenvector
     * @param eigenvalue2 a {@link Field.Complex}: second eigenvalue
     * @param eigenvector2 a {@link Complex2}: second eigenvector
     */
    final public void assignByEigenvectors(final Complex eigenvalue1,
                                           final Complex2 eigenvector1,
                                           final Complex eigenvalue2,
                                           final Complex2 eigenvector2) {
        super.assignByEigenvectors(eigenvalue1, eigenvector1,
                                   eigenvalue2, eigenvector2);
    }



    public void assignAsEigenvectorMatrixOf( HermitianComplex2By2 m ) {

	final double absOfB   = Math.sqrt( m.bRe * m.bRe + m.bIm * m.bIm );

	final double cosOfPsi =   m.bRe / absOfB;
	final double sinOfPsi = - m.bIm / absOfB;

	final double t = - absOfB / ( m.aRe - m.dRe );

	final double l =( t*t ) / ( 4 * t*t + 1 );

	System.out.println( "t=" + t );
	System.out.println( "l=" + l );

	final double cc = 0.5 + Math.sqrt( 0.25 - l );

	final double cosOfPhi = Math.sqrt(     cc );
	final double sinOfPhi = Math.sqrt( 1 - cc );

	assign( cosOfPhi,               0,
		sinOfPhi * (-cosOfPsi), sinOfPhi * ( sinOfPsi), 
		sinOfPhi * ( cosOfPsi), sinOfPhi * ( sinOfPsi), 
		cosOfPhi,               0 );

    }


}
