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

/**
 * This class computes iso lines on triangulations of different
 * 
 * @author schmies
 */
public class MarchingTriangles {

	final int START_CAPACITY = 100;
	
	double [] point;
	int [] index = new int[ 2 *START_CAPACITY];;
	
	int dimension;
	
	int numOfPoints = 0;
	int numOfEdges = 0;
	
	int maxNumOfPoints = START_CAPACITY;
	int maxNumOfEdges = START_CAPACITY;
	
	public MarchingTriangles( int dimension ) {
		
		this.dimension = dimension;
		
	    point = new double[ dimension * START_CAPACITY ];
	}
	
	public double [] getPoints() {
		return point;
	}
	
	public int [] getIndices() {
		return index;
	}
	
	/**
	 * @deprecated
	 * @return
	 */
	public int getNumOfPonits() {
		return numOfPoints;
	}
	public int getNumOfPoints() {
		return numOfPoints;
	}
	public int getNumOfEdges() {
		return numOfEdges;
	}
	
	public void clear() {
		numOfEdges=0;
		numOfPoints=0;
	}
	
	void setMaxNumOfPoints( int maxNumOfPoints ) {
		if( maxNumOfPoints == this.maxNumOfPoints )
			return;
		this.maxNumOfPoints = maxNumOfPoints;
		numOfPoints = Math.min( numOfPoints, maxNumOfPoints );
		double [] newPoint = new double[ maxNumOfPoints * dimension ];
		System.arraycopy( point, 0, newPoint, 0, numOfPoints*dimension );
		point = newPoint;
	}
	
	void setMaxNumOfEdges( int maxNumOfEdges ) {
		if( maxNumOfEdges == this.maxNumOfEdges )
			return;
		this.maxNumOfEdges = maxNumOfEdges;
		numOfEdges = Math.min( numOfEdges, maxNumOfEdges );
		int [] newIndex = new int[ maxNumOfEdges * 2 ];
		System.arraycopy( index, 0, newIndex, 0, numOfEdges*2 );
		index = newIndex;
	}
	
	
	/**
	 * add point ( 1 - tPQ) * P + tPQ * Q to point list.
	 */
	protected void addPointOnEdge( double [] point, int indexP, int indexQ, double tPQ ) {
	
		if( numOfPoints >= maxNumOfPoints  )
			setMaxNumOfPoints( 3 * maxNumOfPoints / 2 + 1);
			
		int IP = dimension * indexP;
		int IQ = dimension * indexQ;
		
		for( int k=0, K=numOfPoints * dimension; k<dimension; k++, K++, IP++, IQ++) {
			this.point[K] =  ( 1 - tPQ) * point[IP] + tPQ * point[IQ];
		}
		
		numOfPoints++;
	}
	
	void addEdge( int indexP, int indexQ ) {
		
		if( numOfEdges >= maxNumOfEdges  )
			setMaxNumOfEdges( 3 * maxNumOfEdges / 2 + 1);
		
		index[ 2*numOfEdges ]    = indexP;
		index[ 2*numOfEdges+1] = indexQ;
		
		numOfEdges++;
	}
	
	public static MarchingTriangles compute( double [] point, int dimension, int [] index, int [] neighbour,
			double [] field, double level) {
		MarchingTriangles il = new MarchingTriangles( dimension);
		il.compute( point, index, neighbour, field, level );
		return il;
	}
	
	public void compute( double [] point, int [] index, int [] neighbour,
			double [] field, double level) {

		final int noe = index.length / 3;

		int [] intersectionIndex = new int[index.length];

		// first loop over all edges, detects all intersections
		for( int i=0, I=0, J=0; i<noe; i++, I+=3 ) {
			for (int j = 0; j < 3; j++, J++ ) {

				// if index of neighour is smaller we take its intersection
				if (neighbour != null && neighbour[J] != -1 && neighbour[J] < i) {
					int iN =3*neighbour[J];
					for( int k=0; k<3; k++ ) {
						if( neighbour[iN+k]==i) {
							intersectionIndex[J] = intersectionIndex[iN+k];
							break;
						}
						if( k==2 ) throw new RuntimeException("wrong connectivity"); 
					}
					continue;
				}	
				
				final int indexP = index[ I + (j + 1) % 3];
				final int indexQ = index[ I + (j + 2) % 3];

				final double fP = field[indexP] - level;
				final double fQ = field[indexQ] - level;

				if (fP * fQ < 0 && fP != fQ) {

					final double tPQ = fP / ( fP-fQ );

					intersectionIndex[J] = this.numOfPoints;

					this.addPointOnEdge( point, indexP, indexQ, tPQ );
					
				} else intersectionIndex[J] = -1;
			}
		}

		this.setMaxNumOfPoints(this.numOfPoints);
	
		// second lopp over all edges generates edges of indexed line set
		for( int i=0, I=0; i<noe; i++, I+=3 ) {
			for (int j = 0; j < 3; j++ ) {

				final int ii = intersectionIndex[ I+(j+1)%3 ];
				final int jj = intersectionIndex[ I+(j+2)%3 ];

				if(  ii != -1 && jj != -1 )
					this.addEdge( ii, jj );
			}
		}

		this.setMaxNumOfEdges(this.numOfEdges);
		
	}
}
