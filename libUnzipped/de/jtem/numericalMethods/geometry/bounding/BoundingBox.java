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

package de.jtem.numericalMethods.geometry.bounding;

/** Provides static method for computing bounding box for point clouds.
 ** @version 1.0
 ** @author Markus Schmies
 */

public class BoundingBox {

    /** Protect constructor since it is static only class. */
    protected BoundingBox() {
    }

    final static double EPS = 1e-14;
    
    private static void resultIsPoint( double [] point, int dim, 
				       double [][] transform,
				       double [] min, double [] max ) {

	System.arraycopy( point, 0, min, 0, dim );
	System.arraycopy( point, 0, max, 0, dim );

	for( int i=0; i<dim; i++ )
	    for( int j=0; j<dim; j++ )
		transform[i][j] = i == j ? 1 : 0;
    }

    /** Completes orthonoormal basis. */
    private static void gramSchmitt( double [][] basis, int k ) {
	
	final int dim = basis.length;

	final double [] n = basis[k];

	for( int l=0; l<dim; l++ )
	    n[l] = 1 + Math.random();
	for( int j=0; j<k; j++ ) {
	    final double [] m = basis[j];
	    
	    double product = 0;
	    for( int l=0; l<dim; l++ )
		product += m[l] * n[l];
	    for( int l=0; l<dim; l++ )
		n[l] -= product * m[l];		
	}

	double lengthSqr = 0;
	for( int l=0; l<dim; l++ )
	    lengthSqr += n[l] * n[l];
	    
	final double length = Math.sqrt( lengthSqr );

	if( length < EPS )
	    throw new RuntimeException( "gram schmitt failed" );
	    
	for( int l=0; l<dim; l++ )
	    n[l] /= length;
    }

    /** Computes bounding box for point cloud in arbitrary dimension.
     ** The coordinates are given pointwise, therefore must the length
     ** of the <code>point</code> array be a multiple of <code>dim</code>
     ** which is the dimension of the surrounding space.
     ** On output the content of <code>point</code> will be destroyed.
     ** The bounding box is given as an interval 
     ** <code>[min[0],max[0]x...x[min[dim-1],max[dim-1]</code> which is
     ** transformed by the orthonormal <code>dim</code> by <code>dim</code>
     ** matrix which is return in <code>transform</code>.
     ** In other words, the inverse tramsformation, which is
     ** is given by the transposed of <code>transform</code>, maps the points
     ** into the above interval.
     ** The length of the interval are sorted starting with the the largest.
     ** The routine supports degenerated cases, where the points lay in a subspace.
     ** The algorithm does not find the smallest bound and delivers the
     ** worst results in the cases of generalized cubes.
     ** @param point array with length which must be multiple of dim, data is destroyed on output
     ** @param dim dimension of space
     ** @param transform orthonormal matrix (first index rows)
     ** @param lower bound of interval
     ** @param upper bound of interval */
    public static void compute( double [] point, int dim, 
				double [][] transform,
				double [] min, double [] max ) {
	
	if( point.length % dim != 0 )
	    throw new IllegalArgumentException("length of point array is not multiple of dim" );

	final int numOfPoints = point.length / dim ;

	for( int k=0; k<dim; k++ ) {
	  
	    // compute biggest deviation in subspace
	    double maxProduct = 0;

	    int maxI = 0;
	    int maxJ = 1;

	    for( int i=0; i<numOfPoints-1; i++ ) {
		final int I = dim * i;
		for( int j=i+1, J=j*dim; j<numOfPoints; j++, J+=dim ) {	       
		    double product = 0;
		    for( int l=0; l<dim; l++ )
			product += ( point[I+l] - point[J+l] ) 
			    * ( point[I+l] - point[J+l] ); 

		    if( maxProduct < product ) {
			maxProduct = product;
			maxI = I;
			maxJ = J;
		    }
		}
	    }		
	  
	    final double [] n = transform[k];

	    if( maxProduct < EPS ) {
		if( k == 0 ) {
		    resultIsPoint( point, dim, transform, min, max );
		    return;
		}

		gramSchmitt( transform, k );
		
		double maxProjection = -Double.MAX_VALUE;
		double minProjection =  Double.MAX_VALUE;
		for( int i=0, I=0; i<numOfPoints; i++, I+=dim ) {
		    double product = 0;
		    for( int l=0; l<dim; l++ )
			product += point[I+l] * n[l];
		    
		    if( product < minProjection )
			minProjection = product;

		    if( product > maxProjection )
			maxProjection = product;
		}
		
		min[k] = minProjection;
		max[k] = maxProjection;
		
	    } else {
		
		final double length = Math.sqrt( maxProduct );
		for( int l=0; l<dim; l++ )
		    n[l] = ( point[ maxI + l ] -  point[ maxJ + l ] ) / length;

		min[k] = 0;
		max[k] = 0;
		for( int l=0; l<dim; l++ ) {
		    min[k] += point[ maxJ + l ] * n[l];
		    max[k] += point[ maxI + l ] * n[l];
		}
	    }
	    
	    // project points to subspace
	    for( int i=0, I=0; i<numOfPoints; i++ ) {
		double product = 0;
		for( int l=0; l<dim; l++ )
		    product += point[ I + l ] * n[l];
		
		for( int l=0; l<dim; l++ )
		    point[I++] -= product * n[l];
	    }
	}
    }
}





