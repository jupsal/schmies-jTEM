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

import java.util.ArrayList;

public class LevelLineUtility {

	
	public static int [][] generateLineSets( int [] edgeIndices ) {
		return generateLineSets( edgeIndices, edgeIndices.length/2);	
	}
	
	public static int [][] generateLineSets( int [] edgeIndices, int numEdges ) {
		
		ArrayList lineList = new ArrayList();
		
		IndexList edges = new IndexList( edgeIndices );
		IndexList line  = new IndexList( edgeIndices.length );
		
		
		while( edges.size() > 0 ) {
			line.removeAll();
			// take last two points from edge list and remove it
			line.append( edges.getLast() ); edges.removeLast();
			line.append( edges.getLast() ); edges.removeLast(); 
			
			stripLine( line, edges );
			
			if( line.get(0) != line.getLast() ) {
				line.reverse();
				stripLine(line, edges );
			}
			lineList.add( line.toArray() );
		}
		
		return toIntArray( lineList );
	}

	private static void stripLine( IndexList line, IndexList edges ) {
		
		while( edges.size() > 0 ) {
			
			final int I = edges.searchForIndexValue( line.getLast() );
			if( I == -1 ) 
				break;
			final int edge = I/2;
			final int other = edges.get(2*edge + ((I%2)+1)%2);
			line.append(other);
			edges.set(2*edge+1, edges.getLast() ); edges.removeLast();
			edges.set(2*edge+0, edges.getLast() ); edges.removeLast();
		}
		
	}
	
	static class IndexList {
		private final int [] index;
		private int numIndices;
		
		IndexList( int maxNumber ) {
			index = new int[maxNumber];
			numIndices=0;
		}
		IndexList( int [] indices ) {
			index = (int[])(indices.clone());
			numIndices=indices.length;
		}
		
		int size() {
			return numIndices;
		}
		
		int getLast() {
			return index[numIndices-1];
		}
		
		int get( int pos ) {
			return index[pos];
		}
		
		void set( int pos, int value ) {
			index[pos]=value;
		}

		void append( int value ) {
			index[numIndices]=value;
			numIndices++;
		}
		
		void removeLast() {
			numIndices--;
		}
		
		void removeAll() {
			numIndices=0;
		}
		
		int searchForIndexValue( int value ) {
			for( int i=0; i<numIndices; i++ ) {
				if( index[i] == value )
					return i;
			}
			
			return -1;
		}
		
		int [] toArray() {
			int [] dst = new int[numIndices];
			System.arraycopy(index,0,dst,0,numIndices);
			return dst;
		}
		
		void reverse() {
			final int n = numIndices/2;
			for( int i=0, j=numIndices-1; i<n; i++, j-- ) {
				final int tmp = index[i];
				index[i]=index[j];
				index[j]=tmp;
			}
		}
	}
	
	static int [][] toIntArray( ArrayList list ) {
		int [][] lineArray = new int[ list.size()][];
		
		for( int i=0; i<lineArray.length; i++ ) {
			lineArray[i] = (int[])list.get(i);
		}
		return lineArray;
	}
	
	public static double lengthOfCurve( double [] points, int dimension ) {
		if( points.length%dimension != 0 )
			throw new IllegalArgumentException( "dimension and size does not match!");
		
		double length=0;
		for( int i=dimension; i<points.length; i+=dimension ) {
			length += dist( points, i-dimension, points, i, dimension );
		}
		return length;
	}
	
	public static double [] resampleCurve( double [] points, int dimension, int N ) {
		if( points.length%dimension != 0 )
			throw new IllegalArgumentException( "dimension and size does not match!");
		
		double [] sample = new double[N*dimension];
		double totalLength = lengthOfCurve( points, dimension );
		double stepSize = totalLength / (N-1);
		
		double length = 0;
		
		System.arraycopy( points,0,sample,0,dimension);
		
		for( int i=dimension, j=1; i<points.length; i+=dimension ) {
			double edgeLength = dist( points, i-dimension, points, i, dimension );
			
			while( j*stepSize >= length && j*stepSize < length+edgeLength ) {
				double weight = (j*stepSize - length) / edgeLength;
				for( int k=0, J=dimension*j; k<dimension; k++, J++ ) {
					sample[J] = weight * points[i] + (1-weight) * points[i-dimension];
				}
				j++;
			}
			
			length += edgeLength;
		}
		

		System.arraycopy( points,points.length-dimension,
				          sample,sample.length-dimension,dimension);
		
		return sample;
	}
	
	static double dist( double [] s1, int os1, double [] s2, int os2, int dimension ) {
		double sum=0;
		for( int k=0; k<dimension; k++, os1++, os2++ ) {
			final double d = s1[os1] - s2[os2];
			sum += d*d;
		}
		return Math.sqrt(sum);
	}
	
	static double dot( double [] s1, int os1, double [] s2, int os2, int dimension ) {
		double sum=0;
		for( int k=0; k<dimension; k++, os1++, os2++ ) {
			final double d = s1[os1]*s2[os2];
			sum += d;
		}
		return sum;
	}
	
	
	
//	
//	static class PointList {
//		double []
//		PointList( double [] points, int dimension ) {
//			
//		}
//	}
	
}
