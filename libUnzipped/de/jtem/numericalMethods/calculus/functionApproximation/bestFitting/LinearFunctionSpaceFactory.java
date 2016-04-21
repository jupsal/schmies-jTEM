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

package de.jtem.numericalMethods.calculus.functionApproximation.bestFitting;

/**
 * @author schmies
 *
 */
public class LinearFunctionSpaceFactory {
	
	LinearFunctionSpaceFactory () {		
	}
	
	/**
	 * Creates sum of two functions spaces. 
	 * Two function spaces can only be added if the parameter and value space of their
	 * elements coinside. 
	 * @param f1
	 * @param f2
	 * @return sum of two spaces.
	 */
	public static LinearFunctionSpace createSum
		(  LinearFunctionSpace f, LinearFunctionSpace g ) {
		if( f.n  != g.n )
			throw new IllegalArgumentException( "function space do not match for this operation");
		if( f.m  != g.m )
			throw new IllegalArgumentException( "function space do not match for this operation");
	
		final int m = f.m;
		
		BasisFunction [][] basis 
			= new BasisFunction[f.dim+g.dim][m];
		
		for( int i=0; i<f.dim; i++ )
			System.arraycopy( f.f[i], 0, basis[i], 0, m);
		for( int i=0, k=f.dim; i<g.dim; i++, k++ )
			System.arraycopy( g.f[i], 0, basis[k], 0, m);
		
		return new LinearFunctionSpace( basis );
	}
	
	/**
	 * Creates power of a function spaces. 
	 * elements coinside. 
	 * @param f
	 * @param power
	 * @return sum of two spaces.
	 */
	public static LinearFunctionSpace createPower
	(  LinearFunctionSpace f, int power ) {
		
		final int m = f.m;
		
		BasisFunction [][] basis
		= new BasisFunction[f.dim][power*m];
		
		for( int i=0; i<f.dim; i++ ) 
			for( int j=0, k=0; j<power; j++, k+=m)
				System.arraycopy( f.f[i], 0, basis[i], k, m);
		
		return new LinearFunctionSpace( basis );
	}
	
	/**
	 * Creates product of two function spaces. 
	 * Two function spaces can only be added if their dimensino and the value space of their
	 * elements coinside. 
	 * @param f
	 * @param g
	 * @return sum of two spaces.
	 */
	public static LinearFunctionSpace createProduct
	(  LinearFunctionSpace f,  LinearFunctionSpace g ) {
		
		if( f.n  != g.n )
			throw new IllegalArgumentException( "function space do not match for this operation");
		if( f.dim  != g.dim )
			throw new IllegalArgumentException( "function space do not match for this operation");
		
		final int n = f.n;
		final int p = f.dim;
		
		BasisFunction [][] basis
		= new BasisFunction[f.dim][f.m+g.m];
		
		for( int i=0; i<f.dim; i++ ) {
				System.arraycopy( f.f[i], 0, basis[i], 0, f.m);
				System.arraycopy( f.f[i], 0, basis[i], f.m, g.m);
		}
		
		return new LinearFunctionSpace( basis );
	}
	

	public static LinearFunctionSpace createPartialDerivative( LinearFunctionSpace f, int k ) {
		
		if( k>=f.n)
			throw new IllegalArgumentException( "wrong dimension");
		
		BasisFunction [][] dBasis = new BasisFunction [f.dim][f.m];
		
		for( int i=0; i<f.f.length; i++ ) {
			for( int j=0; j<f.f[i].length; j++ ) {
				dBasis[i][j] = f.f[i][j].getPartialDerivative(k);
			}
		}
		return new  LinearFunctionSpace( dBasis );
	}
	
	public static LinearFunctionSpace createConstant( int n, int m ) {
		Monomial [][] basis =  new Monomial[1][1];
		basis[0][0] = new Monomial.Constant(n);
		return createPower( new  LinearFunctionSpace( basis ), m );
	}
	
	public static LinearFunctionSpace createLinear( int n, int m ) {
		Monomial [][] basis =  new Monomial[n][1];
		for( int i=0; i<n; i++ ) {
			basis[i][0] = new Monomial.Linear( i, n);
		}
		return createPower( new  LinearFunctionSpace( basis ), m );
	}

	public static LinearFunctionSpace createAffine( int n, int m ) {
		return createPower( createSum( createConstant(n,1), createLinear(n,1) ), m );
	}
	
	public static LinearFunctionSpace createExactOrder2( int n ) {
		Monomial [][] basis =  new Monomial[n*(n+1)/2][1];
		for( int i=0, k=0; i<n; i++ ) {
			for( int j=0; j<=i; j++, k++ ) {
				basis[k][0] = new Monomial.Quadratic( j, i, n);
			}	
		}
		return new LinearFunctionSpace( basis );
	}
	
	public static LinearFunctionSpace createQuadratic( int n, int m ) {
		return createPower( createSum( createAffine(n,1), createExactOrder2(n) ), m );
	}
	
	public static LinearFunctionSpace createExactOrder3( int n ) {
		Monomial [][] basis =  new Monomial[n*(n+1)*(2*n+4)/12][1];
		for( int i=0, l=0; i<n; i++ ) {
			for( int j=0; j<=i; j++) {
				for( int k=0; k<=j; k++, l++) {
					basis[l][0] = new Monomial.Cubic( j, i, k, n);
				}
			}	
		}
		return new LinearFunctionSpace( basis );
	}
		
		public static LinearFunctionSpace createCubic( int n, int m ) {
			return createPower( createSum( createQuadratic(n,1), createExactOrder3(n) ), m );
		}
		
}
