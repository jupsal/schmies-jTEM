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

package de.jtem.mfc.vector;

import java.io.Serializable;

import de.jtem.mfc.field.Complex;
import de.jtem.mfc.matrix.AbstractComplex2By2;

/**
 ** This class represents a complex 2-vector
 ** <pre>
 ** [  a  ]     [ aRe + i aIm ] 
 ** [  b  ]  =  [ bRe + i bIm ]
 ** </pre>
 ** <p>
 ** To minimize the object count its 2 complex entries <code>a,b</code>
 ** are represented by 4 double values 
 ** <code>aRe, aIm, bRe, bIm</code>.
 ** All methods avoid the creation of temporarily used objects in their
 ** internal computations unless stated in their describtions.
 ** <p>
 ** Creation of many temporarily used instances is expensive and stresses 
 ** the garabage collector.
 ** This can cause your java program to be slow and should therefore be avoided.
 ** Thus, operations which result in an instance of <code>Field.Complex</code> 
 ** can be performed either with or without creating a vector.
 ** To compute, for example, the product of matrix with a vector  you can either use
 ** <pre>
 ** Complex2 v = m.times( w ),
 ** </pre>
 ** or
 ** <pre>
 ** v.assignTimes( m, w ).
 ** </pre>
 ** This philosophy is applied to operations which result in other instances,
 ** as well. In such a case an instance of the resulting type
 ** can be prescribed as a parameter which is than filled by
 ** the method.
 ** <p>
 ** @author schmies */

public class Complex2 implements Serializable, Cloneable {

    private static final long serialVersionUID = 1L;

     /** threashold for zero tests; default is 1e-14 */
     protected static final double EPS = 1.e-14;

    /** entry of the vector */
    public double aRe, aIm, bRe, bIm;

    /** creates zero vector */
    public Complex2() { 
    }

    /** creates a vector with the prescribed entries */
    public Complex2( final double aRe, final double aIm, 
		     final double bRe, final double bIm ) {

	assign( aRe, aIm, bRe, bIm );
    }

    /** creates a vector with the prescribed entries */
    protected Complex2( final Complex a, final Complex b ) {
	assign( a.re, a.im, b.re, b.im ); 
    }

    /** creates a vector equal to the prescribed one */
    public Complex2( final Complex2 v ) { 
	assign( v );
    }

    /** assigns this with zero */
    public void assignZero() {
	assign( 0, 0, 0, 0 );
    }

    /** assigns this with the prescribed entries */
    public void assign( final double aRe, final double aIm, 
			   final double bRe, final double bIm ) {
	this.aRe = aRe; this.aIm = aIm;
	this.bRe = bRe; this.bIm = bIm;
    }
    
    /** returns entry <code>a</code> of this */
    final public Complex getA() {
	return new Complex( aRe, aIm );
    }
    
    /** assigns <code>a</code> with entry <code>a</code> of this
     ** @param a complex value which is assigned with entry <code>a</code>. */
    final public void getA( final Complex a ) {
	a.assign( aRe, aIm );
    }

    /** sets entry <code>a</code> of this with <code>a</code> */
    final public void setA( final Complex a ) {
	aRe = a.re;
	aIm = a.im;
    }

    /** returns entry <code>b</code> of this */
    final public Complex getB() {
	return new Complex( bRe, bIm );
    }

    /** assigns <code>b</code> with entry <code>b</code> of this
     ** @param b complex value which is assigned with entry <code>b</code>.*/
    final public void getB( final Complex b ) {
	b.assign( bRe, bIm );
    }

    /** sets entry <code>b</code> of this with <code>b</code> */
    final public void setB( final Complex b ) {
	bRe = b.re;
	bIm = b.im;
    }

    /** returns copy of this */
    final public Complex2 copy() { 
	return new Complex2( this );
    }

    /** assigns this with the prescribed entries */
    final public void assign( final Complex a, final Complex b ) {
	assign( a.re, a.im, b.re, b.im );
    }
    
    /** assigns this with the prescribed vector */
    final public void assign( final Complex2 s ) { 
	assign( s.aRe, s.aIm, s.bRe, s.bIm ); 
    }

    /** assign this with product of matrix <code>m</code> and vector
     * <code>v</code>. */
    final public void assignTimes( final AbstractComplex2By2 m, 
                                   final Complex2 v ) {
	AbstractComplex2By2.times( m, v, this );
    }

    /** 
     * Assign <code>this</code> with the product of a complex number and a
     * complex 2-vector.
     */
    final public void assignTimes(final Complex z, final Complex2 v) {

        /* Use dummy variables because z might be this. */
        final double dummyARe = z.re * v.aRe - z.im * v.aIm;
        final double dummyAIm = z.re * v.aIm + z.im * v.aRe;
        final double dummyBRe = z.re * v.bRe - z.im * v.bIm;
        final double dummyBIm = z.re * v.bIm + z.im * v.bRe;
        
        this.aRe = dummyARe;
        this.aIm = dummyAIm;
        this.bRe = dummyBRe;
        this.bIm = dummyBIm;
    }

    public String toString() {
	StringBuffer sb=new StringBuffer().append('(').append('(').append('(').append(aRe);
	if(aIm>=0.0) 
	    sb.append('+');
	sb.append(aIm).append('i').append(')').append(',').append('(').append(bRe);
	if(bIm>=0.0) 
	    sb.append('+');	
	sb.append(bIm).append('i').append(')').append(')').append(')');

	return sb.toString();
    }

}
