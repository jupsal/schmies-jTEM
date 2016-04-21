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


    /** Transforms a symmetric matrix into triadiagonal form.
     */
public class Tridiagonal {

    private Tridiagonal() {};

    /** Transforms symmetric matrix <code>A</code> into triadiagonal form. 
     *	In <code>Q</code> the orthogonal matrix is returned which realizes   
     *	the change of basis. <code>A</code> and <code>Q</code> may be the same.
     * <p>
     *  @param  A   The symmetric matrix which to transform in triadiagonal form. 
     *  @param  Q   The orthogonal matrix which stores the changes of basis.
     *  @param  d   Consists of the diagonal elements of the tridiagonal matrix.
     *	@param  e   The off-diagonal elements starting from the second element (e[0]=0).
     */
    public static void transform( double [][] A, double [][] Q, double [] d, double [] e ) {

	final int n = A.length;
	
	if( A != Q )
	    for( int i=0; i<n; i++ )
		System.arraycopy( A[i], 0, Q[i], 0, n );
	
	for( int i=n-1;i>=1;i--) {
	    final int l=i-1;

	    double h = 0;
		
	    if (l > 1) {

		double scale = 0;

		for( int k=0;k<=l;k++)
		    scale += Math.abs(Q[i][k]);

		if (scale == 0.0) {
		    e[i]=Q[i][l];
		} else {

		    for( int k=0;k<=l;k++) {
			Q[i][k] /= scale;
			h += Q[i][k]*Q[i][k];
		    }
		    
		    final double f=Q[i][l];
		    final double g=(f >= 0.0 ? -Math.sqrt(h) : Math.sqrt(h));
		    
		    e[i]    = scale*g;
		    Q[i][l] = f-g;
		    h -= f*g;


		    double a=0.0;
		    for( int j=0;j<=l;j++) {
			Q[j][i]=Q[i][j]/h;
			double b=0.0;
			for( int k=0;k<=j;k++)
			    b += Q[j][k]*Q[i][k];
			for( int k=j+1;k<=l;k++)
			    b += Q[k][j]*Q[i][k];
			e[j]=b/h;
			a += e[j]*Q[i][j];
		    }
		    
		    final double aBy2h = a/(h+h);

		    for( int j=0;j<=l;j++) {
			final double b = Q[i][j];
			final double c = e[j] -= aBy2h * b;
			for( int k=0;k<=j;k++)
			    Q[j][k] -= (b*e[k]+c*Q[i][k]);
		    }
		}
	    } else {
		e[i]=Q[i][l];
	    }

	    d[i] = h;
	}

	d[1]=0.0;
	e[1]=0.0;

	for( int i=0;i<n;i++) {
	    final int l=i-1;

	    if (d[i] != 0 ) {
		for( int j=0;j<=l;j++) {
		    double g=0.0;
		    for( int k=0;k<=l;k++)
			g += Q[i][k]*Q[k][j];
		    for( int k=0;k<=l;k++)
			Q[k][j] -= g*Q[k][i];
		}
	    }

	    d[i]   = Q[i][i];
	    Q[i][i]= 1;

	    for( int j=0;j<=l;j++) 
		Q[j][i] = Q[i][j] =0;
	}

	    
    }

    /** Transforms symmetric matrix <code>A</code> into triadiagonal form
     *	which is returned in <code>T</code>.
     *	<code>Q</code> consists of the orthogonal matrix realizing
     *	the change of basis. 
     *	<code>A</code> and <code>Q</code> may be the same.
     *  <p>
     *  @param    A    The symmetric matrix which to transform.
     *  @param    Q    The orthogonal matrix which saves the changes.
     *  @param    T    The matrix to store result in.
     */
    public static void transform( double [][] A, double [][] Q, double [][] T ) {

	final double [] d = T[0];
	final double [] e = T[1];

	transform( A, Q, d, e );

	final int n = A.length;

	// set 1st and 2nd row
	T[1][0] = e[1];
	T[1][1] = d[1]; 
	T[0][1] = T[1][0];

	// set all other rows
	for( int i=2; i<n; i++ ) {
	    T[i][i  ] =             d[ i ];
	    T[i][i-1] = T[i-1][i] = e[ i ];
	}

	// zero the off-tri-diagonal elements
	for( int i=0; i<n; i++ ) {
	    for( int k=0; k<i-1; k++ )
		T[i][k] = 0;
	    for( int k=i+2; k<n; k++ )
		T[i][k] = 0;
	}


    }


}
