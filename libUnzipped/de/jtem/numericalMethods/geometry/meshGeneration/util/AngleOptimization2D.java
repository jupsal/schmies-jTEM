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

package de.jtem.numericalMethods.geometry.meshGeneration.util;

import java.io.Serializable;

import de.jtem.numericalMethods.calculus.function.RealFunctionOfSeveralVariablesWithGradient;
import de.jtem.numericalMethods.calculus.minimizing.ConjugateGradient;

public class AngleOptimization2D
    implements RealFunctionOfSeveralVariablesWithGradient, Serializable {

    double [] length;
    double [] cosAngle;
    double [] edge;

    boolean [] pointIsFix;

    int [] index;

    int numOfElements, numOfPoints;

    double cosMinAngle = .5;

    int numOfEvaluations = 0;

    public AngleOptimization2D( int [] index, boolean [] pointIsFix ) {

	this.index = index;

	this.pointIsFix = pointIsFix;

	numOfPoints   = pointIsFix.length / 2;
	numOfElements = index.length / 3;
    }

    public static boolean [] pointIsFix( int [] index, int [] neighbour, int numOfPoints ) {

	boolean [] pointIsFix = new boolean[ numOfPoints ];

	int numOfElements = index.length / 3;

	for( int i=0, k=0; i < numOfElements; i++, k+=3 )
	    for( int j=0; j<3; j++ )
		if( neighbour[k+j] == -1 ) {
		    pointIsFix[ index[ k+(j+1)%3 ] ] = true;
		    pointIsFix[ index[ k+(j+2)%3 ] ] = true;
		}

	return pointIsFix;
    }

    public AngleOptimization2D( int [] index, int [] neighbour, int numOfPoints ) {

	this( index, pointIsFix( index, neighbour, numOfPoints ) );
    }

    public final double edge( int element, int local, int coord ) {
	return edge[ ( element*3 + local )*2 + coord ];
    }

    public void updateEdge( double [] point ) {

	if( edge == null || edge.length != 6 * numOfElements )
	    edge = new double[ 6 * numOfElements ];

	for( int i=0, n=0; i<numOfElements; i++ )
	    for( int j=0; j<3; j++ )
		for( int k=0; k<2; k++)
		    edge[ n++ ] =
			point[ index[ i*3+(j+2)%3 ] * 2 + k ] -
			point[ index[ i*3+(j+1)%3 ] * 2 + k ];
    }

    public final double length( int element, int local ) {
	return length[ element * 3 + local ];
    }

    public void updateLength() {

	if( length == null || length.length != 3 * numOfElements )
	    length = new double[ 3 * numOfElements ];

	for( int i=0, n=0; i<numOfElements; i++ )
	    for( int j=0; j<3; j++ )
		length[n++] = Math.sqrt( edge(i,j,0) * edge(i,j,0) + edge(i,j,1) * edge(i,j,1) );
    }

    public double cosAngle( int element, int local ) {
	return cosAngle[ element * 3 + local ];
    }

    public void updateCosAngle() {

	if( cosAngle == null || cosAngle.length != 3 * numOfElements )
	    cosAngle = new double[ 3 * numOfElements ];

	for( int i=0, n=0; i<numOfElements; i++ )
	    for( int j=0; j<3; j++ ) {
		final double a = length( i, (j+2)%3 );
		final double b = length( i, (j+1)%3 );
		final double c = length( i,  j      );

		cosAngle[ i * 3 + j ] = ( a*a + b*b - c*c ) / ( 2*a*b );
	    }
    }


    public int getNumberOfVariables() {
      return numOfPoints * 2;
    }

    public double eval(double[] p, double[] grad) {
      setDoubleArrayParameter(p, 0);
      getDoubleArrayValue(grad, 0);
      return getDoubleValue();
    }

    public double eval(double[] p) {
      setDoubleArrayParameter(p, 0);
      return getDoubleValue();
    }

    public double getDoubleValue() {

	numOfEvaluations++;

	double energie = 0;

	for( int i=0; i<numOfElements; i++ )
	    for( int j=0; j<3; j++ ) {
		final double localForce = cosAngle(i,j) - cosMinAngle;

		if( localForce > 0 )
		    energie += localForce * localForce;
	    }

	return energie / 2;
    }


    public int getDoubleArrayParameterLength() {
	return numOfPoints * 2;
    }

    public void setDoubleArrayParameter( final double [] p, final int offset ) {

	if( offset != 0 )
	    throw new RuntimeException( "offset is not supported" );

	updateEdge( p );
	updateLength();
	updateCosAngle();
    }

    public int getDoubleArrayValueLength() {
	return numOfPoints * 2;
    }

    public void getDoubleArrayValue( final double [] grad, final int offset ) {

	final int maxIndex = offset+numOfPoints*2;
	for( int i=offset; i<maxIndex; i++ )
	    grad[i] = 0;

	final int numOfElements = index.length / 3;

	for( int i=0; i<numOfElements; i++ ) {

	    for( int j=0; j<3; j++ ) {

		final double localForce = cosAngle(i,j) - cosMinAngle;

		if( localForce > 0 ) {

		    final int thisVertex = index[ 3*i +  j      ];
		    final int nextVertex = index[ 3*i + (j+1)%3 ];
		    final int prevVertex = index[ 3*i + (j+2)%3 ];

		    final double a = length( i, (j+2)%3 );
		    final double b = length( i, (j+1)%3 );

		    final double cFactor = localForce / ( a*b );
		    final double aFactor = ( 1 - cosAngle(i,j) * b / a ) * cFactor;
		    final double bFactor = ( 1 - cosAngle(i,j) * a / b ) * cFactor;

		    for( int k=0; k<2; k++ ) {
			grad[ offset + nextVertex * 2 + k ] += edge( i, (j+2)%3, k ) * aFactor;
			grad[ offset + thisVertex * 2 + k ] -= edge( i, (j+2)%3, k ) * aFactor;

			grad[ offset + thisVertex * 2 + k ] += edge( i, (j+1)%3, k ) * bFactor;
			grad[ offset + prevVertex * 2 + k ] -= edge( i, (j+1)%3, k ) * bFactor;

			grad[ offset + prevVertex * 2 + k ] -= edge( i,  j,      k ) * cFactor;
			grad[ offset + nextVertex * 2 + k ] += edge( i,  j,      k ) * cFactor;
		    }
		}
	    }
	}

	nilFixedPoints( grad );
    }

    void nilFixedPoints( double [] grad ) {
	// nil grad at fixed points
	for( int i=0, k=0; i<numOfPoints; i++, k+=2 )
	    if( pointIsFix[ i ] )
		grad[k] = grad[k+1] = 0;
    }

    public double getMinAngle( double [] p ) {

	setDoubleArrayParameter( p, 0 );

	double maxCosAngle = -1;

	for( int i=0, k=0; i<numOfElements; i++ )
	    for( int j=0; j<3; j++, k++ )
		if( cosAngle[k] > maxCosAngle )
		    maxCosAngle = cosAngle[k];

	return Math.acos( maxCosAngle );
    }

    static public void minimize( final double [] p,
				 final int [] index,
				 final int [] neighbour,
				 final double minAngle,
				 final int n, final double ftol ) {

	AngleOptimization2D f = new AngleOptimization2D( index, neighbour, p.length );

	f.cosMinAngle = Math.cos( minAngle / 180 * Math.PI );

	System.out.println( "minimal angle before minimizing: "
			    + ( f.getMinAngle( p ) / Math.PI * 180) );

	ConjugateGradient.search( p, ftol, f, n, false, null);

	System.out.println( "number of evaluations : " + f.numOfEvaluations );

	System.out.println( "minimal angle after minimizing:  "
			    +  ( f.getMinAngle( p ) / Math.PI * 180) );
    }

    static public boolean minimize( final double [] p,
				 final int [] index,
				 final int [] neighbour,
				 final int n, final double ftol ) {

	AngleOptimization2D f = new AngleOptimization2D( index, neighbour, p.length );

	f.cosMinAngle = Math.cos( f.getMinAngle( p ) + Math.PI / 180 );

	double initMinAngle = f.getMinAngle( p ) / Math.PI * 180;

	System.out.println( "minimal angle before minimizing: " + initMinAngle );

	ConjugateGradient.search( p, ftol, f, n, false, null );

	double finalMinAngle = f.getMinAngle( p ) / Math.PI * 180;

	System.out.println( "minimal angle after minimizing:  " + finalMinAngle );

	return finalMinAngle - 0.1 > initMinAngle;

    }

    static public void minimize( final double [] p,
				 final int [] index,
				 final int [] neighbour,
				 final int n, final int m, final double ftol ) {

	for( int i=0; i<m; i++ )
	    if( !minimize( p, index, neighbour, n, ftol ) )
		break;
    }
}
