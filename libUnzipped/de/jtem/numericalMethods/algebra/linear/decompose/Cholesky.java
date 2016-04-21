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

package de.jtem.numericalMethods.algebra.linear.decompose;

/** Computes the cholesky decomposition of a symmetric positive-definite
 *  Matrix and can therefore test whether a symmetric matrix is positive-definite.
 */

public class Cholesky {
    
    private Cholesky() {};

    private static IllegalArgumentException notPositiveDefinte = 
	new IllegalArgumentException("matrix is not positive definite");

    /** Computes the cholesky decomposition of matrix <code>A</code> and stores
     *	the result in <code>L</code> . 
     *	<code>L</code> may be null and also coincide with <code>A</code> .
     *	The method only uses the upper right part of <code>A</code> and does not
     *	test its symmetrie.
     *	The method considers the length of <code>A</code> as the dimension of the
     *	space and does not check any ranges.
     *  <p>
     *  @param	A   The matrix to decompose.
     *  @param	L   The matrix to store result in. 
     *  @exception IllgealArgumentException if A is not positive-definite
     */

    public static void decompose( double [][] A, double [][] L ) {

	int n = A.length;

	for( int i=0; i<n; i++ ) {
	    double [] ithRowOfL = L[i];
	    double [] ithRowOfA = A[i];

	    for( int j=i; j<n; j++ ) {
		double [] jthRowOfL = L[j];
	    
		double sum = ithRowOfA[j];
		for( int k=i-1; k>=0; k-- )
		    sum -= ithRowOfL[k] * jthRowOfL[k];

		if( i == j ) {
		    if( sum <= 0 )
			throw notPositiveDefinte;
		    
		    if( L != null )
			L[j][i] = Math.sqrt( sum );

		} else {
		    if( L != null )
			L[j][i] = sum / L[i][i]; 
		}		
		
		if( i<j )
		    L[i][j] = 0;  // zero 
	    }
	}
    }

    /** Checks if <code>A</code> has a cholesky decomposition, which is, if
     *	<code>A</code> is symmetric and positive definite.
     *  For the symmetric part it allows a tolerance of  <code>eps</code>.
     *	To check for positive definiteness it uses 
     *	{@link #decompose(double[][],double[][])}.
     *	The method considers the length of <code>A</code> as the dimension of the
     *	space and does not check any ranges.
     *  <p>
     *  @param   A     The matrix to check
     *  @param   eps   The tolerance of check
     */

    public static boolean decomposable( double [][] A, double eps ) {
	int n = A.length;

	for( int i=0; i<n; i++ )
	    for( int j=i+1; j<n; j++ )
		if( Math.abs( A[i][j] - A[j][i] ) > eps )
		    return false;
    
	try {
	    decompose( A, null );
	} catch( IllegalArgumentException e ) {
	    return false;
	}
	
	return true;
    }

    /** Invokes {@link #decomposable(double [][], double ) } with a tolerance
     *	of  <code>eps = 1e-14</code>.
     *  <p>
     *  @param    A   The Matrix to check decomposition
     */
    public static boolean decomposable( double [][] A ) {
	return decomposable( A, 1e-14 );
    }
}
