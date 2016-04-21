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

package de.jtem.numericalMethods.geometry.hyperSurface;

import java.io.Serializable;

public class MarchingCubes implements Serializable {

    private static final long serialVersionUID = 1L;
    
    double [] point;
    int    [] index;

    double [] min;
    double [] max;

    public MarchingCubes() {

	min = new double [] { -0.5, -0.5, -0.5 };
	max = new double [] {  0.5,  0.5,  0.5 };

    }

    public MarchingCubes( double [][][] field, double level, double [] min, double [] max ) {
	this.min = min;
	this.max = max;

	triangulate( field, level );
    }

    public double getMin( int index ) {
	return min[index];
    }

    public void setMin( int index, double min ) {
	if( this.min[index] == min )
	    return;


	if( this.max[index] <= min )
		throw new IllegalArgumentException(" max must be bigger then min");
	
	rescale( index, min, this.min[index], this.max[index], this.max[index] );
	
	this.min[index] = min;
    }

    public double getMax( int index ) {
	return max[index];
    }

    public void setMax( int index, double max ) {
	if( this.max[index] == max )
	    return;
	
	if( this.min[index] >= max )
		throw new IllegalArgumentException(" max must be bigger then min");
	
	rescale( index, this.min[index], this.min[index], max, this.max[index] );

	this.max[index] = max;
    }

    private void rescale( int index, 
			  double newMin, double oldMin, 
			  double newMax, double oldMax ) {

    if( point == null ) return;
    
	final double factor = ( newMax - newMin ) / ( oldMax - oldMin );
	final double offset = newMin - oldMin * factor;
	
	for( int i=index; i<point.length; i+=3 )
	    point[i] = offset + point[i] * factor;
    }

    public double [] getPoints() {
	return point;
    }

    public int [] getIndices() {
	return index;
    }

    static abstract class RealFunctionOnTripleIndex {
    	int xSize, ySize, zSize;

    	abstract double valueAt( int i, int j, int k );
    }
   
    public void triangulate( final double[][][] field, final double level ) {

		RealFunctionOnTripleIndex function = new RealFunctionOnTripleIndex() {

			{
				xSize = field      .length;
				ySize = field[0]   .length;
				zSize = field[0][0].length;
			}

			double valueAt( int i, int j, int k ) {
				return field[i][j][k];
			}
		};
		
		triangulate( function, level );
    }
    
    public static interface  RealFunctionOnReal3 {
    	public double valueAt( double x, double y, double z );
    }
    
    public void triangulate( final RealFunctionOnReal3 f, final int [] sizes, final double level ) {

    	RealFunctionOnTripleIndex function = new RealFunctionOnTripleIndex() {
    		{
    			xSize = sizes[0];
    			ySize = sizes[1];
    			zSize = sizes[2];
    		}

    		final double xs = 1./(xSize-1);
    		final double ys = 1./(ySize-1);
    		final double zs = 1./(zSize-1);
    		
    		final double xMin = min[0];
    		final double yMin = min[1];
    		final double zMin = min[2];
    		
    		final double xMax = max[0];
    		final double yMax = max[1];
    		final double zMax = max[2];
    		
    		double valueAt( int i, int j, int k ) {
    			final double x = xMin * ( 1 - i * xs ) + xMax * i * xs;
    			final double y = yMin * ( 1 - j * ys ) + yMax * j * ys;
    			final double z = zMin * ( 1 - k * zs ) + zMax * k * zs;
    			
    			return f.valueAt(x,y,z);
    		}
    	};
    	
    	triangulate( function, level );
    }
    
    
	void triangulate( final RealFunctionOnTripleIndex field, final double level ) {

		final int xSize = field.xSize;
		final int ySize = field.ySize;
		final int zSize = field.zSize;
		
	final int [][][] xIndex = new int[ xSize-1 ][ ySize   ][ zSize   ];
	final int [][][] yIndex = new int[ xSize   ][ ySize-1 ][ zSize   ];
	final int [][][] zIndex = new int[ xSize   ][ ySize   ][ zSize-1 ];

	int numOfPoints = 0;

	// computing z-lines-intersections
	for( int i=0; i<xSize; i++ ) {
	    for( int j=0; j<ySize; j++ ) {

		double p = field.valueAt(i,j,0); // field[i][j][0];

		for( int k=1; k<zSize; k++ ) {

		    double q = field.valueAt(i,j,k); // field[i][j][k];

		    if( q < level && p >= level || p < level && q>= level ) {

			final double t = ( level - p ) / ( q - p );  
		    
			final double z = k-1 + t;
		    
			setPoint( numOfPoints, i, j, z );

			zIndex[i][j][k-1] = numOfPoints++;
		    } else {
			zIndex[i][j][k-1] = -1;
		    }

		    p = q;
		}
	    }
	}


	// computing x-lines-intersections
	for( int j=0; j<ySize; j++ ) {
	    for( int k=0; k<zSize; k++ ) {

		double p = field.valueAt(0,j,k); //field[0][j][k];

		for( int i=1; i<xSize; i++ ) {

		    double q = field.valueAt(i,j,k); //field[i][j][k];

		    if( q < level && p >= level || p < level && q >= level ) {

			final double t = (   level - p ) / ( q - p );  
		    
			final double x = i-1 + t;
		    
			setPoint( numOfPoints, x, j, k );

			xIndex[i-1][j][k] = numOfPoints++;
		    } else {
			xIndex[i-1][j][k] = -1;
		    }

		    p = q;
		}
	    }
	}


	// computing y-lines-intersections
	for( int i=0; i<xSize; i++ ) {
	    for( int k=0; k<zSize; k++ ) {

		double p = field.valueAt(i,0,k); // field[i][0][k];

		for( int j=1; j<ySize; j++ ) {

		    double q = field.valueAt(i,j,k); //field[i][j][k];

		    if( q < level && p >= level || p < level && q>= level ) {
		    
			final double t = ( level - p ) / ( q - p );  
		    
			final double y = j-1 + t;
		    
			setPoint( numOfPoints, i, y, k );

			yIndex[i][j-1][k] = numOfPoints++;
		    } else {
			yIndex[i][j-1][k] = -1;
		    }

		    p = q;
		}
	    }
	}

	setNumOfPoints( numOfPoints );

	int numOfTriangles = 0;

	int [] vertexOnEdge = new int[12];

	for( int i=0; i<xSize-1; i++ ) {
	    for( int j=0; j<ySize-1; j++ ) {
		for( int k=0; k<zSize-1; k++ ) {
	
		    int cubeindex = 0;

		    if( field.valueAt(i  ,j  ,k  ) < level) cubeindex |= 1;   
		    if( field.valueAt(i+1,j  ,k  ) < level) cubeindex |= 2;
		    if( field.valueAt(i+1,j  ,k+1) < level) cubeindex |= 4;
		    if( field.valueAt(i  ,j  ,k+1) < level) cubeindex |= 8;
		    if( field.valueAt(i  ,j+1,k  ) < level) cubeindex |= 16;
		    if( field.valueAt(i+1,j+1,k  ) < level) cubeindex |= 32;
		    if( field.valueAt(i+1,j+1,k+1) < level) cubeindex |= 64;
		    if( field.valueAt(i  ,j+1,k+1) < level) cubeindex |= 128;
//
//		    if( field[i  ][j  ][k  ] < level) cubeindex |= 1;   
//		    if( field[i+1][j  ][k  ] < level) cubeindex |= 2;
//		    if( field[i+1][j  ][k+1] < level) cubeindex |= 4;
//		    if( field[i  ][j  ][k+1] < level) cubeindex |= 8;
//		    if( field[i  ][j+1][k  ] < level) cubeindex |= 16;
//		    if( field[i+1][j+1][k  ] < level) cubeindex |= 32;
//		    if( field[i+1][j+1][k+1] < level) cubeindex |= 64;
//		    if( field[i  ][j+1][k+1] < level) cubeindex |= 128;
//		    

		    /* Cube is entirely in/out of the surface */
		    if( cubeindex == 0 || cubeindex == 255 )
			continue;
	
		    // edges in x-direction
		    vertexOnEdge[ 0] = xIndex[i  ][j  ][k  ];
		    vertexOnEdge[ 2] = xIndex[i  ][j  ][k+1];
		    vertexOnEdge[ 4] = xIndex[i  ][j+1][k  ];
		    vertexOnEdge[ 6] = xIndex[i  ][j+1][k+1];

		    // edges in y-direction
		    vertexOnEdge[ 8] = yIndex[i  ][j  ][k  ];
		    vertexOnEdge[ 9] = yIndex[i+1][j  ][k  ];
		    vertexOnEdge[10] = yIndex[i+1][j  ][k+1];
		    vertexOnEdge[11] = yIndex[i  ][j  ][k+1];

		    // edges in z-direction
		    vertexOnEdge[ 1] = zIndex[i+1][j  ][k  ];
		    vertexOnEdge[ 3] = zIndex[i  ][j  ][k  ];
		    vertexOnEdge[ 5] = zIndex[i+1][j+1][k  ];
		    vertexOnEdge[ 7] = zIndex[i  ][j+1][k  ];
		    
		    /* Create the triangles */
		    for( int s=0; s<triTable[cubeindex].length; s+=3, numOfTriangles++ ) {
			setTriangle( numOfTriangles, 
				     vertexOnEdge[ triTable[cubeindex][s  ] ],
				     vertexOnEdge[ triTable[cubeindex][s+1] ],
				     vertexOnEdge[ triTable[cubeindex][s+2] ] );		    }
		}
	    }
	}

	setNumOfTriangles( numOfTriangles );

	rescale( 0, min[0], 0, max[0], xSize-1 );
	rescale( 1, min[1], 0, max[1], ySize-1 );
	rescale( 2, min[2], 0, max[2], zSize-1 );
    }
    
    private void setPoint( int pos, double x, double y, double z ) {

	if( point == null || point.length <= 3 * pos ) {

	    double [] newPoint = new double[ 3 * pos + 300 ];

	    if( point != null )
		System.arraycopy( point, 0, newPoint, 0, point.length );

	    point = newPoint;
	}
	
	point[ 3*pos     ] = x;
	point[ 3*pos + 1 ] = y;
	point[ 3*pos + 2 ] = z;
    }

    private void setNumOfPoints( int nop ) {
	
	if( point != null && point.length == 3*nop )
	    return;
	
	double [] newPoint = new double[ 3 * nop ];
	
	if( point != null )
	    System.arraycopy( point, 0, newPoint, 0, Math.min( 3*nop, point.length ) );
	
	point = newPoint;
    }
    
    private void setTriangle( int pos, int p, int q, int r ) {

	if( index == null || index.length <= 3 * pos ) {

	    int [] newIndex = new int[ 3 * pos + 300 ];

	    if( index != null )
		System.arraycopy( index, 0, newIndex, 0, index.length );

	    index = newIndex;
	}
	
	index[ 3*pos     ] = p;
	index[ 3*pos + 1 ] = q;
	index[ 3*pos + 2 ] = r;
    }

    private void setNumOfTriangles( int not ) {

	if( index != null && index.length == 3*not )
	    return;
	
	int [] newIndex = new int[ 3 * not ];

	if( index != null )
	    System.arraycopy( index, 0, newIndex, 0, Math.min( 3*not, index.length ) );

	index = newIndex;
    }

    static int [][] triTable = new int[][] {
	{},
	{0, 8, 3},
	{0, 1, 9},
	{1, 8, 3, 9, 8, 1},
	{1, 2, 10},
	{0, 8, 3, 1, 2, 10},
	{9, 2, 10, 0, 2, 9},
	{2, 8, 3, 2, 10, 8, 10, 9, 8},
	{3, 11, 2},
	{0, 11, 2, 8, 11, 0},
	{1, 9, 0, 2, 3, 11},
	{1, 11, 2, 1, 9, 11, 9, 8, 11},
	{3, 10, 1, 11, 10, 3},
	{0, 10, 1, 0, 8, 10, 8, 11, 10},
	{3, 9, 0, 3, 11, 9, 11, 10, 9},
	{9, 8, 10, 10, 8, 11},
	{4, 7, 8},
	{4, 3, 0, 7, 3, 4},
	{0, 1, 9, 8, 4, 7},
	{4, 1, 9, 4, 7, 1, 7, 3, 1},
	{1, 2, 10, 8, 4, 7},
	{3, 4, 7, 3, 0, 4, 1, 2, 10},
	{9, 2, 10, 9, 0, 2, 8, 4, 7},
	{2, 10, 9, 2, 9, 7, 2, 7, 3, 7, 9, 4},
	{8, 4, 7, 3, 11, 2},
	{11, 4, 7, 11, 2, 4, 2, 0, 4},
	{9, 0, 1, 8, 4, 7, 2, 3, 11},
	{4, 7, 11, 9, 4, 11, 9, 11, 2, 9, 2, 1},
	{3, 10, 1, 3, 11, 10, 7, 8, 4},
	{1, 11, 10, 1, 4, 11, 1, 0, 4, 7, 11, 4},
	{4, 7, 8, 9, 0, 11, 9, 11, 10, 11, 0, 3},
	{4, 7, 11, 4, 11, 9, 9, 11, 10},
	{9, 5, 4},
	{9, 5, 4, 0, 8, 3},
	{0, 5, 4, 1, 5, 0},
	{8, 5, 4, 8, 3, 5, 3, 1, 5},
	{1, 2, 10, 9, 5, 4},
	{3, 0, 8, 1, 2, 10, 4, 9, 5},
	{5, 2, 10, 5, 4, 2, 4, 0, 2},
	{2, 10, 5, 3, 2, 5, 3, 5, 4, 3, 4, 8},
	{9, 5, 4, 2, 3, 11},
	{0, 11, 2, 0, 8, 11, 4, 9, 5},
	{0, 5, 4, 0, 1, 5, 2, 3, 11},
	{2, 1, 5, 2, 5, 8, 2, 8, 11, 4, 8, 5},
	{10, 3, 11, 10, 1, 3, 9, 5, 4},
	{4, 9, 5, 0, 8, 1, 8, 10, 1, 8, 11, 10},
	{5, 4, 0, 5, 0, 11, 5, 11, 10, 11, 0, 3},
	{5, 4, 8, 5, 8, 10, 10, 8, 11},
	{9, 7, 8, 5, 7, 9},
	{9, 3, 0, 9, 5, 3, 5, 7, 3},
	{0, 7, 8, 0, 1, 7, 1, 5, 7},
	{1, 5, 3, 3, 5, 7},
	{9, 7, 8, 9, 5, 7, 10, 1, 2},
	{10, 1, 2, 9, 5, 0, 5, 3, 0, 5, 7, 3},
	{8, 0, 2, 8, 2, 5, 8, 5, 7, 10, 5, 2},
	{2, 10, 5, 2, 5, 3, 3, 5, 7},
	{7, 9, 5, 7, 8, 9, 3, 11, 2},
	{9, 5, 7, 9, 7, 2, 9, 2, 0, 2, 7, 11},
	{2, 3, 11, 0, 1, 8, 1, 7, 8, 1, 5, 7},
	{11, 2, 1, 11, 1, 7, 7, 1, 5},
	{9, 5, 8, 8, 5, 7, 10, 1, 3, 10, 3, 11},
	{5, 7, 0, 5, 0, 9, 7, 11, 0, 1, 0, 10, 11, 10, 0},
	{11, 10, 0, 11, 0, 3, 10, 5, 0, 8, 0, 7, 5, 7, 0},
	{11, 10, 5, 7, 11, 5},
	{10, 6, 5},
	{0, 8, 3, 5, 10, 6},
	{9, 0, 1, 5, 10, 6},
	{1, 8, 3, 1, 9, 8, 5, 10, 6},
	{1, 6, 5, 2, 6, 1},
	{1, 6, 5, 1, 2, 6, 3, 0, 8},
	{9, 6, 5, 9, 0, 6, 0, 2, 6},
	{5, 9, 8, 5, 8, 2, 5, 2, 6, 3, 2, 8},
	{2, 3, 11, 10, 6, 5},
	{11, 0, 8, 11, 2, 0, 10, 6, 5},
	{0, 1, 9, 2, 3, 11, 5, 10, 6},
	{5, 10, 6, 1, 9, 2, 9, 11, 2, 9, 8, 11},
	{6, 3, 11, 6, 5, 3, 5, 1, 3},
	{0, 8, 11, 0, 11, 5, 0, 5, 1, 5, 11, 6},
	{3, 11, 6, 0, 3, 6, 0, 6, 5, 0, 5, 9},
	{6, 5, 9, 6, 9, 11, 11, 9, 8},
	{5, 10, 6, 4, 7, 8},
	{4, 3, 0, 4, 7, 3, 6, 5, 10},
	{1, 9, 0, 5, 10, 6, 8, 4, 7},
	{10, 6, 5, 1, 9, 7, 1, 7, 3, 7, 9, 4},
	{6, 1, 2, 6, 5, 1, 4, 7, 8},
	{1, 2, 5, 5, 2, 6, 3, 0, 4, 3, 4, 7},
	{8, 4, 7, 9, 0, 5, 0, 6, 5, 0, 2, 6},
	{7, 3, 9, 7, 9, 4, 3, 2, 9, 5, 9, 6, 2, 6, 9},
	{3, 11, 2, 7, 8, 4, 10, 6, 5},
	{5, 10, 6, 4, 7, 2, 4, 2, 0, 2, 7, 11},
	{0, 1, 9, 4, 7, 8, 2, 3, 11, 5, 10, 6},
	{9, 2, 1, 9, 11, 2, 9, 4, 11, 7, 11, 4, 5, 10, 6},
	{8, 4, 7, 3, 11, 5, 3, 5, 1, 5, 11, 6},
	{5, 1, 11, 5, 11, 6, 1, 0, 11, 7, 11, 4, 0, 4, 11},
	{0, 5, 9, 0, 6, 5, 0, 3, 6, 11, 6, 3, 8, 4, 7},
	{6, 5, 9, 6, 9, 11, 4, 7, 9, 7, 11, 9},
	{10, 4, 9, 6, 4, 10},
	{4, 10, 6, 4, 9, 10, 0, 8, 3},
	{10, 0, 1, 10, 6, 0, 6, 4, 0},
	{8, 3, 1, 8, 1, 6, 8, 6, 4, 6, 1, 10},
	{1, 4, 9, 1, 2, 4, 2, 6, 4},
	{3, 0, 8, 1, 2, 9, 2, 4, 9, 2, 6, 4},
	{0, 2, 4, 4, 2, 6},
	{8, 3, 2, 8, 2, 4, 4, 2, 6},
	{10, 4, 9, 10, 6, 4, 11, 2, 3},
	{0, 8, 2, 2, 8, 11, 4, 9, 10, 4, 10, 6},
	{3, 11, 2, 0, 1, 6, 0, 6, 4, 6, 1, 10},
	{6, 4, 1, 6, 1, 10, 4, 8, 1, 2, 1, 11, 8, 11, 1},
	{9, 6, 4, 9, 3, 6, 9, 1, 3, 11, 6, 3},
	{8, 11, 1, 8, 1, 0, 11, 6, 1, 9, 1, 4, 6, 4, 1},
	{3, 11, 6, 3, 6, 0, 0, 6, 4},
	{6, 4, 8, 11, 6, 8},
	{7, 10, 6, 7, 8, 10, 8, 9, 10},
	{0, 7, 3, 0, 10, 7, 0, 9, 10, 6, 7, 10},
	{10, 6, 7, 1, 10, 7, 1, 7, 8, 1, 8, 0},
	{10, 6, 7, 10, 7, 1, 1, 7, 3},
	{1, 2, 6, 1, 6, 8, 1, 8, 9, 8, 6, 7},
	{2, 6, 9, 2, 9, 1, 6, 7, 9, 0, 9, 3, 7, 3, 9},
	{7, 8, 0, 7, 0, 6, 6, 0, 2},
	{7, 3, 2, 6, 7, 2},
	{2, 3, 11, 10, 6, 8, 10, 8, 9, 8, 6, 7},
	{2, 0, 7, 2, 7, 11, 0, 9, 7, 6, 7, 10, 9, 10, 7},
	{1, 8, 0, 1, 7, 8, 1, 10, 7, 6, 7, 10, 2, 3, 11},
	{11, 2, 1, 11, 1, 7, 10, 6, 1, 6, 7, 1},
	{8, 9, 6, 8, 6, 7, 9, 1, 6, 11, 6, 3, 1, 3, 6},
	{0, 9, 1, 11, 6, 7},
	{7, 8, 0, 7, 0, 6, 3, 11, 0, 11, 6, 0},
	{7, 11, 6},
	{7, 6, 11},
	{3, 0, 8, 11, 7, 6},
	{0, 1, 9, 11, 7, 6},
	{8, 1, 9, 8, 3, 1, 11, 7, 6},
	{10, 1, 2, 6, 11, 7},
	{1, 2, 10, 3, 0, 8, 6, 11, 7},
	{2, 9, 0, 2, 10, 9, 6, 11, 7},
	{6, 11, 7, 2, 10, 3, 10, 8, 3, 10, 9, 8},
	{7, 2, 3, 6, 2, 7},
	{7, 0, 8, 7, 6, 0, 6, 2, 0},
	{2, 7, 6, 2, 3, 7, 0, 1, 9},
	{1, 6, 2, 1, 8, 6, 1, 9, 8, 8, 7, 6},
	{10, 7, 6, 10, 1, 7, 1, 3, 7},
	{10, 7, 6, 1, 7, 10, 1, 8, 7, 1, 0, 8},
	{0, 3, 7, 0, 7, 10, 0, 10, 9, 6, 10, 7},
	{7, 6, 10, 7, 10, 8, 8, 10, 9},
	{6, 8, 4, 11, 8, 6},
	{3, 6, 11, 3, 0, 6, 0, 4, 6},
	{8, 6, 11, 8, 4, 6, 9, 0, 1},
	{9, 4, 6, 9, 6, 3, 9, 3, 1, 11, 3, 6},
	{6, 8, 4, 6, 11, 8, 2, 10, 1},
	{1, 2, 10, 3, 0, 11, 0, 6, 11, 0, 4, 6},
	{4, 11, 8, 4, 6, 11, 0, 2, 9, 2, 10, 9},
	{10, 9, 3, 10, 3, 2, 9, 4, 3, 11, 3, 6, 4, 6, 3},
	{8, 2, 3, 8, 4, 2, 4, 6, 2},
	{0, 4, 2, 4, 6, 2},
	{1, 9, 0, 2, 3, 4, 2, 4, 6, 4, 3, 8},
	{1, 9, 4, 1, 4, 2, 2, 4, 6},
	{8, 1, 3, 8, 6, 1, 8, 4, 6, 6, 10, 1},
	{10, 1, 0, 10, 0, 6, 6, 0, 4},
	{4, 6, 3, 4, 3, 8, 6, 10, 3, 0, 3, 9, 10, 9, 3},
	{10, 9, 4, 6, 10, 4},
	{4, 9, 5, 7, 6, 11},
	{0, 8, 3, 4, 9, 5, 11, 7, 6},
	{5, 0, 1, 5, 4, 0, 7, 6, 11},
	{11, 7, 6, 8, 3, 4, 3, 5, 4, 3, 1, 5},
	{9, 5, 4, 10, 1, 2, 7, 6, 11},
	{6, 11, 7, 1, 2, 10, 0, 8, 3, 4, 9, 5},
	{7, 6, 11, 5, 4, 10, 4, 2, 10, 4, 0, 2},
	{3, 4, 8, 3, 5, 4, 3, 2, 5, 10, 5, 2, 11, 7, 6},
	{7, 2, 3, 7, 6, 2, 5, 4, 9},
	{9, 5, 4, 0, 8, 6, 0, 6, 2, 6, 8, 7},
	{3, 6, 2, 3, 7, 6, 1, 5, 0, 5, 4, 0},
	{6, 2, 8, 6, 8, 7, 2, 1, 8, 4, 8, 5, 1, 5, 8},
	{9, 5, 4, 10, 1, 6, 1, 7, 6, 1, 3, 7},
	{1, 6, 10, 1, 7, 6, 1, 0, 7, 8, 7, 0, 9, 5, 4},
	{4, 0, 10, 4, 10, 5, 0, 3, 10, 6, 10, 7, 3, 7, 10},
	{7, 6, 10, 7, 10, 8, 5, 4, 10, 4, 8, 10},
	{6, 9, 5, 6, 11, 9, 11, 8, 9},
	{3, 6, 11, 0, 6, 3, 0, 5, 6, 0, 9, 5},
	{0, 11, 8, 0, 5, 11, 0, 1, 5, 5, 6, 11},
	{6, 11, 3, 6, 3, 5, 5, 3, 1},
	{1, 2, 10, 9, 5, 11, 9, 11, 8, 11, 5, 6},
	{0, 11, 3, 0, 6, 11, 0, 9, 6, 5, 6, 9, 1, 2, 10},
	{11, 8, 5, 11, 5, 6, 8, 0, 5, 10, 5, 2, 0, 2, 5},
	{6, 11, 3, 6, 3, 5, 2, 10, 3, 10, 5, 3},
	{5, 8, 9, 5, 2, 8, 5, 6, 2, 3, 8, 2},
	{9, 5, 6, 9, 6, 0, 0, 6, 2},
	{1, 5, 8, 1, 8, 0, 5, 6, 8, 3, 8, 2, 6, 2, 8},
	{1, 5, 6, 2, 1, 6},
	{1, 3, 6, 1, 6, 10, 3, 8, 6, 5, 6, 9, 8, 9, 6},
	{10, 1, 0, 10, 0, 6, 9, 5, 0, 5, 6, 0},
	{0, 3, 8, 5, 6, 10},
	{10, 5, 6},
	{11, 5, 10, 7, 5, 11},
	{11, 5, 10, 11, 7, 5, 8, 3, 0},
	{5, 11, 7, 5, 10, 11, 1, 9, 0},
	{10, 7, 5, 10, 11, 7, 9, 8, 1, 8, 3, 1},
	{11, 1, 2, 11, 7, 1, 7, 5, 1},
	{0, 8, 3, 1, 2, 7, 1, 7, 5, 7, 2, 11},
	{9, 7, 5, 9, 2, 7, 9, 0, 2, 2, 11, 7},
	{7, 5, 2, 7, 2, 11, 5, 9, 2, 3, 2, 8, 9, 8, 2},
	{2, 5, 10, 2, 3, 5, 3, 7, 5},
	{8, 2, 0, 8, 5, 2, 8, 7, 5, 10, 2, 5},
	{9, 0, 1, 5, 10, 3, 5, 3, 7, 3, 10, 2},
	{9, 8, 2, 9, 2, 1, 8, 7, 2, 10, 2, 5, 7, 5, 2},
	{1, 3, 5, 3, 7, 5},
	{0, 8, 7, 0, 7, 1, 1, 7, 5},
	{9, 0, 3, 9, 3, 5, 5, 3, 7},
	{9, 8, 7, 5, 9, 7},
	{5, 8, 4, 5, 10, 8, 10, 11, 8},
	{5, 0, 4, 5, 11, 0, 5, 10, 11, 11, 3, 0},
	{0, 1, 9, 8, 4, 10, 8, 10, 11, 10, 4, 5},
	{10, 11, 4, 10, 4, 5, 11, 3, 4, 9, 4, 1, 3, 1, 4},
	{2, 5, 1, 2, 8, 5, 2, 11, 8, 4, 5, 8},
	{0, 4, 11, 0, 11, 3, 4, 5, 11, 2, 11, 1, 5, 1, 11},
	{0, 2, 5, 0, 5, 9, 2, 11, 5, 4, 5, 8, 11, 8, 5},
	{9, 4, 5, 2, 11, 3},
	{2, 5, 10, 3, 5, 2, 3, 4, 5, 3, 8, 4},
	{5, 10, 2, 5, 2, 4, 4, 2, 0},
	{3, 10, 2, 3, 5, 10, 3, 8, 5, 4, 5, 8, 0, 1, 9},
	{5, 10, 2, 5, 2, 4, 1, 9, 2, 9, 4, 2},
	{8, 4, 5, 8, 5, 3, 3, 5, 1},
	{0, 4, 5, 1, 0, 5},
	{8, 4, 5, 8, 5, 3, 9, 0, 5, 0, 3, 5},
	{9, 4, 5},
	{4, 11, 7, 4, 9, 11, 9, 10, 11},
	{0, 8, 3, 4, 9, 7, 9, 11, 7, 9, 10, 11},
	{1, 10, 11, 1, 11, 4, 1, 4, 0, 7, 4, 11},
	{3, 1, 4, 3, 4, 8, 1, 10, 4, 7, 4, 11, 10, 11, 4},
	{4, 11, 7, 9, 11, 4, 9, 2, 11, 9, 1, 2},
	{9, 7, 4, 9, 11, 7, 9, 1, 11, 2, 11, 1, 0, 8, 3},
	{11, 7, 4, 11, 4, 2, 2, 4, 0},
	{11, 7, 4, 11, 4, 2, 8, 3, 4, 3, 2, 4},
	{2, 9, 10, 2, 7, 9, 2, 3, 7, 7, 4, 9},
	{9, 10, 7, 9, 7, 4, 10, 2, 7, 8, 7, 0, 2, 0, 7},
	{3, 7, 10, 3, 10, 2, 7, 4, 10, 1, 10, 0, 4, 0, 10},
	{1, 10, 2, 8, 7, 4},
	{4, 9, 1, 4, 1, 7, 7, 1, 3},
	{4, 9, 1, 4, 1, 7, 0, 8, 1, 8, 7, 1},
	{4, 0, 3, 7, 4, 3},
	{4, 8, 7},
	{9, 10, 8, 10, 11, 8},
	{3, 0, 9, 3, 9, 11, 11, 9, 10},
	{0, 1, 10, 0, 10, 8, 8, 10, 11},
	{3, 1, 10, 11, 3, 10},
	{1, 2, 11, 1, 11, 9, 9, 11, 8},
	{3, 0, 9, 3, 9, 11, 1, 2, 9, 2, 11, 9},
	{0, 2, 11, 8, 0, 11},
	{3, 2, 11},
	{2, 3, 8, 2, 8, 10, 10, 8, 9},
	{9, 10, 2, 0, 9, 2},
	{2, 3, 8, 2, 8, 10, 0, 1, 8, 1, 10, 8},
	{1, 10, 2},
	{1, 3, 8, 9, 1, 8},
	{0, 9, 1},
	{0, 3, 8},
	{} };

}

