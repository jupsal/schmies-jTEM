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

package de.jtem.numericalMethods.geometry.latticeReduction;


public class LLL {

    private static final double EPS = 1e-13;

    static int dim = 0;

    static double [][] A;
    static double [][] mu;
    static double []   a;

    
    private static void setDimension( int n ) {

	if( n != dim ) {

	    dim = n;

	    A  = new double[dim][dim];
	    mu = new double[dim][dim];
	    a  = new double[dim];
	}
    }

    private static void initLLL( final double [][] B, final int [][] U, final double [][] U_, 
				 final double [][] mu, final double [] a ) {

	int n = B.length;
	
	if( U != null ) // set U to id
	    for( int i=0; i<n; i++ )
		for( int j=0; j<n; j++ )
		    U[i][j] = i == j ? 1 : 0;

	if( U_ != null ) // set U_ to id
	    for( int i=0; i<n; i++ )
		for( int j=0; j<n; j++ )
		    U_[i][j] = i == j ? 1 : 0;

	// copy B into A
	for( int i=0; i<n; i++ )
	    for( int j=0; j<n; j++ )
		A[i][j] = B[i][j]; //
	
	// init mu and a
	for( int i=0; i<n; i++ ) {

	    for( int j=0; j<i; j++ ) {

		double t = mu[i][j] = dot( B[i], A[j] )  / a[j];

		for( int k=0; k<n; k++ )
		    A[i][k] -= t * A[j][k]; 
	    }

	    a[i] = dot( A[i], A[i] );
	}
    }

    public static void reduce( final double [][] B ) {
	reduce( B, null, null, .75 );
    }

    public static void reduce( final double [][] B, final int [][] U ) {
	reduce( B, U, null, .75 );
    }

    public static void reduce( final double [][] B, final double [][] U ) {
	reduce( B, null, U, .75 );
    }

    public static void reduce( final double [][] B, final int [][] U, final double [][] U_, double c ) {
	
	if( c < 0.25 || c > 1 )
	    throw new IllegalArgumentException( "constant is not in [1/4,1]" );

	int n = B.length;

	setDimension( n );

	initLLL( B, U, U_, mu, a );

	if( n < 2 )
	    return;

	int k = 1;
	
	while( true ) {

	    transform( B, U, U_, mu, k, k-1 ); // l=k-1

	    if( a[k] + EPS < ( c - mu[k][k-1] * mu[k][k-1] ) * a[k-1] ) {
		
		swap( B, U, U_, mu, a, k );

		if( k > 1 )
		    k--;

	    } else {

		for( int l=k-2; l>=0; l-- )
		    transform( B, U, U_, mu, k, l );
				
		k++;

		if( k == n )
		    return;
	    }
	}
    }

    private static void swap( final double [][] B,  final int [][] U, final double [][] U_, 
			      final double [][] mu, final double [] a, final int k ) {

	final int n = B.length;

	final double MU = mu[k][k-1] ;

	double tmp = a[k] + MU * MU * a[k-1];

	mu[k][k-1] *= a[k-1] / tmp;
	a [k]      *= a[k-1] / tmp;
	
	a[k-1] = tmp;

	swap( B, k, k-1 );

	if( U != null )
	    swap( U, k, k-1 );

	if( U_ != null )
	    swap( U_, k, k-1 );

	for( int j=0; j<k-1; j++ ) {
	    tmp        = mu[k-1][j];
	    mu[k-1][j] = mu[k  ][j];
	    mu[k  ][j] = tmp;	    
	}

	for( int i=k+1; i<n; i++ ) {

	    double t2 = mu[i][k-1] - MU * mu[i][k  ];
	    double t1 = mu[i][k  ] + t2 * mu[k][k-1];

	    mu[i][k-1] = t1;
	    mu[i][k  ] = t2;
	}

    }

    private static void transform( final double [][] B,  final int [][] U, final double [][] U_,  
				   final double [][] mu, final int k, final int l ) {

	if( Math.abs( mu[k][l] ) > 0.5 ) {

	    final int r = (int)Math.floor( mu[k][l] + 0.5 );
	    
	    final int n = B.length;
	    
	    for( int i=0; i<n; i++ ) {
		B[k][i] -= r * B[l][i];

		if( U != null )
		    U[k][i] -= r * U[l][i];

		if( U_ != null )
		    U_[k][i] -= r * U_[l][i];
	    }

	    for( int j=0; j<l; j++ )
		mu[k][j] -= r * mu[l][j];
	    
	    mu[k][l] -= r;
	}
    }

    private static double dot( final double [] a, final double [] b ) {
		
	double result = 0;

	for( int i=0; i<a.length; i++ )
	    result += a[i] * b[i];

	return result;
    }

    private static void swap( double [][] a, int k, int l ) {
	double [] tmp = a[k];

	a[k] = a[l];
	a[l] = tmp;
    }

    private static void swap( int [][] a, int k, int l ) {
	int [] tmp = a[k];

	a[k] = a[l];
	a[l] = tmp;
    }

    private static void print( double [] a, String name ) {
	
	System.out.print( name );

	    for( int i=0; i<a.length; i++ )
		System.out.print( "  " + a[i] );

	System.out.println( "  " );
    }

}
