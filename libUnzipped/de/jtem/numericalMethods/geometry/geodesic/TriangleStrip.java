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

package de.jtem.numericalMethods.geometry.geodesic;

/**
 * @author schmies
 *
 * To change the template for this generated type comment go to
 * Window - Preferences - Java - Code Generation - Code and Comments
 */
public class TriangleStrip {

	TriangleStrip() {
		super();
	}


	/**
	 * Connects two given triangles with triangle strip.
	 * @param elementData array of triangle indices given as triple integers.
	 * @param neighbourData array of neighbour indices given as triple integers.
	 * @param e1 index of first triangle.
	 * @param e2 index of secound triangle.
	 * @return triangle strip between two given element indices, which can be null if no 
	 */
	public static int [] getStripBetweenElements( int [][] elementData, int [][] neighbourData, int e1, int e2  ) {
		
		int noe = elementData.length;
	
		int prevElement, nextElement = -1, i, j, e, frontCounter = 0, nFsize, pFsize = 0;

		if(e1 >= 0 && e1 < noe && e2 >= 0 && e2 < noe )
			throw new IllegalArgumentException(
				"elements are not part of triangulations");

		int [] pFront = new int[noe];
		int [] nFront = new int[noe];
		int [] prev = new int[noe];
		
		for( i=0; i<noe; prev[i++] = -1 );
		
		pFront[ pFsize++ ] = e1;

		while( nextElement != e2 && pFsize > 0 ) {

			frontCounter++;

			nFsize = 0;

			for( j=0; j<pFsize; j++ ) {

				prevElement =  pFront[j];
				
				for( i=0; i<3; i++ ) {
					if( ( nextElement = neighbourData[ prevElement ][i] ) != -1 && prev[nextElement] == -1 ) {
						nFront[ nFsize++ ] = nextElement;
						prev[nextElement] = prevElement;
					}
					
					if( nextElement == e2 )
						break;
				}
				
				if( nextElement == e2 )
					break;
			}

			int [] dummy  = pFront;
			pFront = nFront;    pFsize = nFsize;
			nFront = dummy;

		}

		if( nextElement != e2 )
			return null;
		
		int [] strip = new int[frontCounter+1];

		e = strip[ i = frontCounter ] = e2;

		while( i > 0 )
			e = strip[ --i ] = prev[ e ];
		
		return strip;
	}



//	- (OoArrayi *)getElementWaveFromElement: (int)anElement {
//		
//		int noe = [self getNumElements];
//		
//		int *pFront, *nFront, *dummy, *prev, *waveData; 
//		int prevElement, nextElement = -1, i, j, nFsize, pFsize = 0, waveSize = 1;
//
//		OsVec3i *neighbourData = [neighbour getData];
//
//		OoArrayi *aPFront, *aNFront, *aPrev, *aWave, *aDummy;
//
//		ASSERT( neighbourData , OE_MISC, OE_DEFAULT, 
//				"no neighbour information", return nil );
//		ASSERT( anElement >= 0 && anElement < noe , OE_MISC, OE_DEFAULT,
//				"element is not part of triangulations", return nil);
//
//		pFront = [ [ aPFront = [OoArrayi new] setUsedSize: noe] getData ];
//		nFront = [ [ aNFront = [OoArrayi new] setUsedSize: noe] getData ];
//		prev   = [ [ aPrev   = [OoArrayi new] setUsedSize: noe] getData ];
//		
//		for( i=0; i<noe; prev[i++] = -1 );
//
//		prev[ anElement ] = pFront[ pFsize++ ] = anElement;
//		
//		[ aWave = [ aPFront deepCopy ] setUsedSize: 1 ];
//
//		while( pFsize > 0 ) {
//			
//			nFsize = 0;
//
//			for( j=0; j<pFsize; j++ ) {
//				
//				prevElement =  pFront[j];
//				
//				for( i=0; i<3; i++ ) {
//					if( ( nextElement = neighbourData[ prevElement ][i] ) != -1 && prev[nextElement] == -1 ) {
//						nFront[ nFsize++ ] = nextElement;
//						prev[nextElement] = prevElement;
//					}
//				}
//			}
//
//			/* append the new front to the wave */
//			waveData = [ [aWave setUsedSize: waveSize + nFsize ] getData ];
//
//			ooMemCopy( waveData + waveSize, nFront, nFsize * sizeof(int) ); 
//
//			waveSize += nFsize;
//
//			/* the new front become the previous and vice versa */
//			dummy  = pFront;
//			pFront = nFront;    pFsize = nFsize;
//			nFront = dummy;
//
//			aDummy  = aPFront;
//			aPFront = aNFront;
//			aNFront = aDummy;
//
//		}
//
//		[aPFront free];
//		[aNFront free];
//		[aPrev   free];
//
//		return [ aWave decrRefCount ];
//	}
//
//
//
//
//
//	- (OoArrayVec3 *)getMetricOfStrip: (OoArrayi *)strip  {
//
//		OoArrayVec3 *metric = [OoArrayVec3 new];
//
//		if( ![self computeMetricOfStrip: strip into: metric ] ) {
//			[metric free];
//			return nil;
//		}
//
//		return [metric decrRefCount];
//	}
//
//	- computeMetricOfStrip: (OoArrayi *)strip into: (OoArrayVec3 *)metric {
//		
//		OsVec3i *elementData;
//		OsVec3  *metricData;
//		int *stripData, noe, i, NOE = [self getNumElements];
//
//		char    *vertexData = [vertex getData];
//		char *A, *B, *C;
//
//		OoLinAlgFunctions *f = [vertex getLinAlgFunctions];
//
//		int vsize = [vertex getChunkSize];
//		
//		ASSERT( strip, OE_MISC, OE_DEFAULT, "no strip", return nil );
//		ASSERT( metric, OE_MISC, OE_DEFAULT, "no array for metric", return nil );
//		ASSERT( element, OE_MISC, OE_DEFAULT, "no list with elements", return nil );
//		ASSERT( f, OE_MISC, OE_DEFAULT, "vertices have wrong fibers", return nil );
//
//		[strip incrRefCount];
//
//		[metric setUsedSize: noe = [strip getUsedSize] ]; 
//
//		metricData  = [metric getData];
//		stripData   = [strip getData];
//		elementData = [element getData];
//
//		for( i=0; i < noe; i++ )
//			if ( stripData[i] >= 0 && stripData[i]<NOE ) {
//				A = vertexData + elementData[ stripData[i] ][0] * vsize;
//				B = vertexData + elementData[ stripData[i] ][1] * vsize;
//				C = vertexData + elementData[ stripData[i] ][2] * vsize;
//				
//				metricData[i][0] = f->dist( B, C );
//				metricData[i][1] = f->dist( A, C );
//				metricData[i][2] = f->dist( B, A );
//			} else {
//				metricData[i][0] =  1;
//				metricData[i][1] =  1;
//				metricData[i][2] =  1;
//			}
//
//		[strip free];
//
//		return self;
//	}
//
//	/*@DOC**************************************************************************
//	 * 		
//	 * Name:	- checkStrip: (OoArrayi *)strip
//	 * 		
//	 * Description:	optimize the strip
//	 * 		
//	 * Arguments:	OoArrayi *strip
//	 * 		
//	 * Keywords:	
//	 * 		
//	 * Author:	Markus Schmies
//	 * 		
//	 * History:	created: 
//	 *         	revised: 
//	 * 		
//	 *@END**************************************************************************
//	 */
//
//	- checkStrip: (OoArrayi *)strip {
//
//		OoTriangShape *out = self;
//
//		int i, noe, type;
//
//		OsVec3i *elementData, map;
//
//		int *stripData;
//
//		ASSERT( strip, OE_MISC, OE_DEFAULT, "no strip", return nil );
//
//		if( (noe=[strip getUsedSize]) < 2 )
//			return self;
//
//		stripData     = [ strip getData ];
//
//		elementData   = [ element getData ];
//
//		for( i=0; i<noe-1; i++ ) {
//
//			type = guVec3iGetMap( map, elementData[stripData[i]], elementData[stripData[i+1]] );
//
//			if( !type ) {
//				
//				out = nil;
//
//				printf("polygon on triang is not connected between node %d and %d \n", i, i+1 );
//			}
//		}
//
//		return out;
//	}
//
//	- (BOOL)isClosedStrip: (OoArrayi *)strip ofType: (int)type {
//
//		int noe, *stripData;
//		OsVec3i *elementData, map;
//
//		ASSERT( strip, OE_MISC, OE_DEFAULT, "no strip", return NO );
//
//		if( (noe=[strip getUsedSize]) < 2 )
//			return YES;
//
//		stripData     = [ strip getData ];
//		elementData   = [ element getData ];
//
//		if( guVec3iGetMap( map, elementData[stripData[0]], elementData[stripData[noe-1]] ) < type )
//			return NO;
//
//		return YES;
//	}
//
//	- (BOOL)isClosedStrip: (OoArrayi *)strip {
//		return [self isClosedStrip: strip ofType: 1];
//	}
//
//	- (BOOL)isSimpleClosedStrip: (OoArrayi *)strip {
//		return [self isClosedStrip: strip ofType: 2];
//	}
//
//	- (BOOL)isLoopStrip: (OoArrayi *)strip {
//		int stripLen = [strip getUsedSize];
//
//		if(   stripLen > 1
//				&& [strip getIntegerAt: 0] == [strip getIntegerAt: stripLen - 1 ]
//											   && [self isClosedStrip: strip ofType: 2]  )
//			return YES;
//
//		return NO;
//	}
	/**
	 * Returns smallest connected triangle strip within the given one, but without
	 * changing the positions where the strip is only connected via an isolate vertex.
	 * The array of the input strip is used optimized strip has the length of the input array,
	 * neither way the input data is lost. 
	 * @param elementData
	 * @param stripData
	 * @return optimized triangle strip.
	 */
	public static int [] getOptimizedStrip( int [][] elementData, int [] stripData ) {
	
		int noe = optimizeStrip( elementData, stripData, 2  );
		
		if( noe == stripData.length )
			return stripData;
		
		int [] optimizedStrip = new int[ noe ];
		
		System.arraycopy(stripData,0,optimizedStrip,0,noe);
		
		return optimizedStrip;
	}
	
	/**
	 * Returns smallest connected triangle strip within the given one, but with
	 * possibly adding positions where the strip is only connected via an isolate vertex.
	 * The array of the input strip is used optimized strip has the length of the input array,
	 * neither way the input data is lost. 
	 * @param elementData
	 * @param stripData
	 * @return optimized triangle strip.
	 */
	public static int [] getComplexOptimizedStrip( int [][] elementData, int [] stripData ) {
		
		int noe = optimizeStrip( elementData, stripData, 1  );
		
		if( noe == stripData.length )
			return stripData;
		
		int [] optimizedStrip = new int[ noe ];
		
		System.arraycopy(stripData,0,optimizedStrip,0,noe);
		
		return optimizedStrip;
	}
	
	/**
	 * Optimizes strip using "the" sofisticated algo.
	 * @param elementData
	 * @param stripData result on output
	 * @param threshType (complex==1,normal==2)
	 * @return number of triangles used after optimization used after optimzation.
	 */
	static int  optimizeStrip( int [][] elementData, int [] stripData, int threshType ) {
	
		int  noe = stripData.length;

		int [] map = new int[3];

		while( noe > 1 ) {
			int This = 0;
			int thisElement = stripData[ This ];
			
			for( int next=1; next < noe; next++ ) {

				int nextElement = stripData[ next ];

				int type = IntTriple.getMap( map, elementData[thisElement], elementData[nextElement] );

				if( type == 0 )
					throw new RuntimeException( "strip is not connected" );

				if( next+1 < noe ) {

					int nextType = IntTriple.getMap( map, elementData[nextElement], elementData[ stripData[ next+1 ] ] );       


					for( int test=next+1; test < noe; test++ ) {
						
						int testElement = stripData[ test ];

						int testType = IntTriple.getMap( map, elementData[thisElement], elementData[testElement] );

						if( testType==0 )
							break;

						if( testType > threshType || testType >= type || ( test==next+1 && testType >= nextType ) ) {
							next        = test;
							nextElement = testElement;
							type        = testType;
						}
					}
				}

				stripData[ ++This ] = thisElement = nextElement;

			}

			if( noe == This + 1 )
				break;
			
			noe = This + 1;
		}

		return noe;
	}
	
	
	/**
	 * Removes sequentially doubly occuring triangles from strip.
	 * The array of the input strip is used for ouptut optimized strip has the length of the input array,
	 * neither way the input data is lost. 
	 * @param elementData
	 * @param stripData
	 * @return optimized triangle strip.
	 */
	public static int [] getSimpleOptimizedStrip( int [][] elementData, int [] stripData ) {
		
		int  noe = stripData.length;
		
		while( noe > 1 ) {

			int This = 0;
			int thisElement = stripData[ This ];
			
			int next=1;

			while ( next < noe ) {

				int nextElement = stripData[ next ];
				
				if( nextElement == thisElement ) {
					next++;
					continue;
				}

				int test = next+1;

				if( test < noe ) {

					int testElement = stripData[ test ];
					
					if( testElement != thisElement ) {
						
						stripData[ ++This ] = thisElement = nextElement;
						
						next = test;
						
					} else {

						if( This > 0 ) {

							thisElement = stripData[ --This ];

							next = test;

						} else {

							next = test+1;
						}
					}
				} else {

					stripData[ ++This ] = nextElement;

					break;
				}
			}

			if( noe == This + 1 )
				break;
			
			noe = This + 1;
		}
		
		if( noe == stripData.length )
			return stripData;
		
		int [] optimizedStrip = new int[ noe ];
		
		System.arraycopy(stripData,0,optimizedStrip,0,noe);
		
		return optimizedStrip;
	}

//
//	/*@DOC**************************************************************************
//	 * 		
//	 * Name:	developStrip: (OoArrayi *)strip into: (OoTriangShape *)flat
//	 * 		
//	 * Description:	
//	 * 		
//	 * Arguments:	OoArrayi *strip
//	 *              OoTriangShape *flat
//	 * 		
//	 * Keywords:	
//	 * 		
//	 * Author:	Markus Schmies
//	 * 		
//	 * History:	created: 
//	 *         	revised: 
//	 * 		
//	 *@END**************************************************************************
//	 */
//
//	- developStrip: (OoArrayi *)strip into: (OoTriangShape *)flat {
//
//		int i, thisElement, nextElement, noe, type, Vertex = 0, index, lIndex, rIndex, lVertex, rVertex;
//
//		double phi, *lengths;
//
//		OsVec3i map;
//		OsVec3  angles;
//
//		OsVec3  *flatVertex, *metricData;
//		OsVec3i *elementData, *flatElement, *flatMap;
//
//		int *stripData, *this, *next;
//
//		OoArrayVec3 *metric;
//
//		ASSERT( element, OE_MISC, OE_DEFAULT, "no list with elements", return nil );
//		ASSERT( strip, OE_MISC, OE_DEFAULT, "no strip", return nil );
//		ASSERT( flat, OE_MISC, OE_DEFAULT, "no TriangShape for flattened strip", return nil ); 
//
//		if( (noe=[strip getUsedSize]) < 2 )
//			return self;
//
//		[strip incrRefCount];
//
//		metric = [self getMetricOfStrip: strip];
//
//		ASSERT( metric,  OE_MISC, OE_DEFAULT, "could not compute metric of strip", { [strip free]; return nil; } );
//		
//		if( ![flat getVertices] || strcmp( [ [flat getVertices] getCode ], "ddd" ) )
//			[ flat setVertices: [OoArrayVec3 newOrphan] ];
//		if( ![flat getElements] )
//			[ flat setElements: [OoArrayVec3i newOrphan ] ];
//		if( ![flat getNeighbours] )
//			[ flat setNeighbours: [OoArrayVec3i newOrphan ] ];
//		
//
//
//		metricData = [metric getData];
//		stripData  = [strip getData];
//		
//		elementData = [element getData];
//
//		[flat setNumVertices: noe * 2 + 1];
//		[flat setNumElements: noe ];
//
//		flatVertex  = [ [flat getVertices] getData ];
//		flatElement = [ [flat getElements] getData ];
//		flatMap     = [ [flat getNeighbours] getData ];
//
//		flatVertex[0][0] = flatVertex[0][1] = 0.;
//		Vertex = 0;
//
//		nextElement = 0;
//		next = elementData[ stripData[ 0 ] ];
//
//		i=-1;
//
//		while( nextElement < noe ) {
//
//			thisElement = nextElement; 
//			this        = next;
//
//			nextElement = thisElement + 1;
//			
//			if( nextElement < noe ) {
//
//				next        = elementData[ stripData[ nextElement ] ];
//
//				type = guVec3iGetMap( map, next, this );
//
//				if( type==3 )
//					/* this and next element are the same */
//					continue;
//			}
//
//			gmAnglesFromLengths( angles, lengths = metricData[ thisElement ] );
//
//			if( thisElement == noe-1 || type == 1 ) {
//				/* this element is last element of strip or 
//				 this element of strip joins just one vertex with the next */
//				
//				if( thisElement < noe-1 ) {
//					/* if it is not the last triangle  */
//					for( index=0; map[index] < 0; index++ );
//					
//					if( i!=-1 ) {
//						/* if it is not the first triangle */
//						lIndex = guVec3iGetLocInd( this, elementData[ stripData[i] ][ flatMap[i][0] ] );
//						
//						if( this[ lIndex ] == next[ map[index] ] )
//							/* the prev element joins also just one vertex with the next */
//							continue;
//						
//					} else lIndex = (index+1)%3;
//				} else {
//					/* if it is the last triangle */
//					if( i!=-1 )
//						/* and not also the first */
//						lIndex = guVec3iGetLocInd( this, elementData[ stripData[i] ][ flatMap[i][0] ] );
//					else 
//						lIndex = 1;
//					index = (lIndex+2)%3;
//				}
//				
//				rIndex = ( (index+lIndex) * 2 ) % 3;
//				
//				flatElement[++i][ 1 ] = Vertex;     flatMap[i][1] = lIndex;
//				flatElement[  i][ 2 ] = Vertex+1;   flatMap[i][2] = rIndex; 
//				flatElement[  i][ 0 ] = Vertex+2;   flatMap[i][0] = index;
//
//				stripData[i] = stripData[thisElement];
//
//				flatVertex[++Vertex][0] = 
//					flatVertex[Vertex-1][0] + cos( angles[ lIndex ] ) * lengths[ index ];
//				flatVertex[  Vertex][1] = 
//					flatVertex[Vertex-1][1] + sin( angles[ lIndex ] ) * lengths[ index ];
//				
//				flatVertex[++Vertex][0] = flatVertex[Vertex-2][0] + lengths[ rIndex ];
//				flatVertex[  Vertex][1] = flatVertex[Vertex-2][1];
//				
//				continue;
//			}
//			
//			ASSERT( type == 2, OE_MISC, OE_DEFAULT, "unregular strip", { [strip free]; return nil; } );
//
//			flatMap[++i][0] = index  = guVec3iGetLocInd( map, -1 );   
//			flatMap[  i][1] = lIndex = (index+1)%3;
//			flatMap[  i][2] = rIndex = (index+2)%3;
//			
//			flatElement[i][0] = Vertex;                              
//			flatElement[i][1] = Vertex + 1; lVertex = this[ lIndex ];
//			flatElement[i][2] = Vertex + 2; rVertex = this[ rIndex ]; 
//
//			stripData[i] = stripData[thisElement];
//			
//			flatVertex[++Vertex][0] = flatVertex[Vertex-1][0] + lengths[ rIndex ];
//			flatVertex[  Vertex][1] = flatVertex[Vertex-1][1];
//			
//			flatVertex[++Vertex][0] = 
//				flatVertex[Vertex-2][0] + cos( angles[ index ] ) * lengths[ lIndex ];
//			flatVertex[  Vertex][1] = 
//				flatVertex[Vertex-2][1] + sin( angles[ index ] ) * lengths[ lIndex ];
//
//			phi = M_PI - angles[ lIndex ];
//			
//			while( 1 ) {
//
//				thisElement = nextElement; 
//				this        = next;
//				
//				nextElement = thisElement + 1;
//				next        = elementData[ stripData[ nextElement ] ];
//
//				if( nextElement >= noe )
//					break;
//				type = guVec3iGetMap( map, next, this );
//
//				if( type==3 )
//					/* this and next element are the same */
//					continue;
//
//				if( type < 2 )
//					break;
//
//				gmAnglesFromLengths( angles, lengths = metricData[thisElement] );
//
//				index = guVec3iGetLocInd( map, -1 );  
//				
//				if( lVertex == this[ index ] ) {
//					/* the right vertex did not change */
//
//					flatElement[++i][ 0 ] = flatElement[i-1][1];
//					flatElement[  i][ 1 ] = Vertex + 1;;
//					flatElement[  i][ 2 ] = flatElement[i-1][2];
//
//					rIndex  = guVec3iGetLocInd( this, rVertex );
//					lIndex  = guVec3iGetOppLocInd( this, rVertex, this[index] );
//					lVertex = this[ lIndex ];
//
//					phi += angles[ rIndex ];
//
//					flatVertex[++Vertex][0] = flatVertex[ flatElement[i][2] ][0] - cos( phi ) * lengths[ index ];
//					flatVertex[  Vertex][1] = flatVertex[ flatElement[i][2] ][1] - sin( phi ) * lengths[ index ];
//
//				} else {
//					/* the left vertex did not change */
//
//					flatElement[++i][ 0 ] = flatElement[i-1][2];
//					flatElement[  i][ 1 ] = flatElement[i-1][1];
//					flatElement[  i][ 2 ] = Vertex + 1;;
//
//					lIndex  = guVec3iGetLocInd( this, lVertex );
//					rIndex  = guVec3iGetOppLocInd( this, lVertex, this[index] );
//					rVertex = this[ rIndex ];
//					
//					phi -= angles[ lIndex ];
//
//					flatVertex[++Vertex][0] = flatVertex[ flatElement[i][1] ][0] + cos( phi ) * lengths[ index ];
//					flatVertex[  Vertex][1] = flatVertex[ flatElement[i][1] ][1] + sin( phi ) * lengths[ index ];
//				}
//
//				stripData[i] = stripData[thisElement];
//
//				flatMap[i][0] = index;
//				flatMap[i][1] = lIndex;
//				flatMap[i][2] = rIndex;
//			}
//
//			gmAnglesFromLengths( angles, lengths = metricData[thisElement] );
//
//			flatElement[++i][ 0 ] = Vertex + 1;
//			flatElement[  i][ 1 ] = flatElement[i-1][2];
//			flatElement[  i][ 2 ] = flatElement[i-1][1];
//
//			flatMap[i][2] = lIndex = guVec3iGetLocInd( this, lVertex );
//			flatMap[i][1] = rIndex = guVec3iGetLocInd( this, rVertex );
//			flatMap[i][0] = ((rIndex+lIndex)*2)%3;
//
//			stripData[i] = stripData[thisElement];
//
//			flatVertex[++Vertex][0] = 
//				flatVertex[ flatElement[i][2] ][0] + cos( phi - angles[ lIndex ] ) * lengths[ rIndex ];
//			flatVertex[  Vertex][1] =  
//				flatVertex[ flatElement[i][2] ][1] + sin( phi - angles[ lIndex ] ) * lengths[ rIndex ];
//		}
//
//		[flat setNumVertices: Vertex+1];
//		[flat setNumElements: i+1 ];
//
//		[strip setUsedSize: i+1 ];
//
//		[metric delete];
//		[strip free];
//
//		return self;
//	}
//
//	#include <polygonOnTriang/OoPolygonOnTriang.h>
//
//	/*@DOC**************************************************************************
//	 * 		
//	 * Name:	-  computeShortestInRegularFlatInto: (OoPolygonOnTriang *)short;
//	 * 		
//	 * Description:	
//	 * 		
//	 * 		
//	 * Keywords:	
//	 * 		
//	 * Author:	Markus Schmies
//	 * 		
//	 * History:	created: 
//	 *         	revised: 
//	 * 		
//	 *@END**************************************************************************
//	 */
//
//	#define gmDotVec2( V, W ) ((V)[0]*(W)[0] + (V)[1]*(W)[1] )
////#define gmAbsVec2( V )    sqrt( (V)[0]*(V)[0] +  (V)[1]*(V)[1] )
//
//			- computeShortestInRegularFlatInto: (OoPolygonOnTriang *)shortest {
//
//				int nov, noe;
//				int thisElement, startElement, endElement, lElement, rElement;
//				int thisVertex,  startVertex, endVertex, lVertex, rVertex;
//
//				OsVec3i *flatElement;
//				OsVec3 *flatVertex;
//				OsVec2 L, R, T, pL, pR;
//				double lenL, lenR, lenT, *start, eps = 1e-5;
//
//				ASSERT( shortest, OE_MISC, OE_DEFAULT, "shortest is nil", return nil );
//
//				[shortest setNumVertices: 0 ];
//
//				flatVertex  = [ vertex getData ];
//				flatElement = [ element getData ];
//
//				nov = [ self getNumVertices ];
//				noe = [ self getNumElements ];
//
//				startVertex = startElement = 0;
//
//				while( startElement < noe ) {
//
//					/*    endVertex = startVertex + 1;  to avoid errors in strange situations */
//
//					if( flatElement[ startElement ][0] != startVertex ) {
//						/* startElement joins only a vertex with next element */
//						
//						thisVertex = endVertex  = flatElement[ startElement ][ 0 ];
//						endElement = startElement;
//
//					} else {
//
//						start = flatVertex[ flatElement[ startElement ][0] ];
//						
//						lVertex = flatElement[ startElement ][1]; lElement = startElement;
//						rVertex = flatElement[ startElement ][2]; rElement = startElement;
//						
//						thisVertex  = guVec3iGetMax( flatElement[startElement] ) + 1;
//						thisElement = startElement + 1;
//						
//						gmSubVec2( L, flatVertex[lVertex], start ); lenL = gmAbsVec2( L );
//						gmSubVec2( R, flatVertex[rVertex], start ); lenR = gmAbsVec2( R );
//						
//						gmMakeVec2( pL, L[1], -L[0] );
//						gmMakeVec2( pR, R[1], -R[0] );
//
//						if( gmDotVec2( pL, R ) < 0 ) {
//							/* multiply pL with -1 */
//							gmMakeVec2( pL, -L[1], L[0] );
//						}	
//						if( gmDotVec2( pR, L ) < 0 ) {
//							/* multiply pR with -1 */
//							gmMakeVec2( pR, -R[1], R[0] );
//						}
//						
//						while( thisVertex < nov ) {
//							gmSubVec2( T, flatVertex[thisVertex], start ); lenT = gmAbsVec2( T );
//
//							
//							if(   flatElement[ thisElement ][ 1 ] == thisVertex
//									|| flatElement[ thisElement ][ 0 ] == thisVertex ) {
//								/* thisVertex is on the left side */
//								
//								if( gmDotVec2( pL, T ) > eps * lenL * lenT ) { 
//									/* this Vertex is a new horizon on the left side */
//									
//									if( gmDotVec2( pR, T ) < eps * lenR * lenT ) {
//										/* this Vertex on the left side is hidden by R */
//										endVertex  = rVertex;
//										endElement = rElement;
//										break;
//									}	  
//									lVertex  = thisVertex;
//									lElement = thisElement;
//									gmCopyVec2( L, T ); lenL = lenT;  
//
//									gmMakeVec2( pL, L[1], -L[0] );
//
//									if( gmDotVec2( pL, R ) < 0 ) {
//										/* multiply pL with -1 */
//										gmMakeVec2( pL, -L[1], L[0] );
//									}	
//									
//								}
//							}
//							if(   flatElement[ thisElement ][ 2 ] == thisVertex
//									|| flatElement[ thisElement ][ 0 ] == thisVertex ) {
//								/* thisVertex is on the right side */
//								
//								if( gmDotVec2( pR, T ) > eps * lenR * lenT ) {
//									/* this Vertex is a new horizon on the right side */
//									
//									if( gmDotVec2( pL, T ) < eps * lenL * lenT ) { 
//										/* this Vertex on the right side is hidden by L */
//										endVertex  = lVertex;
//										endElement = lElement;
//										break;
//									}
//									rVertex  = thisVertex;
//									rElement = thisElement;
//									gmCopyVec2( R, T ); lenR = lenT;
//
//									gmMakeVec2( pR, R[1], -R[0] );
//
//									if( gmDotVec2( pR, L ) < 0 ) {
//										/* multiply pR with -1 */
//										gmMakeVec2( pR, -R[1], R[0] );
//									}		    
//								}
//							}
//							
//							if( flatElement[ thisElement ][ 0 ] == thisVertex ) {
//								/* thisVertex is at the end of a regular strip */
//								endVertex  = thisVertex;
//								endElement = thisElement;
//								break;
//							}
//							
//							thisVertex++;
//							thisElement++;
//						}
//					}
//
//					[self computeLineFromVertex: startVertex inElement: startElement toVertex: endVertex appendTo: shortest ];
//					
//					startElement = endElement + 1;
//					startVertex  = endVertex;
//					
//					if( endVertex != thisVertex )
//						/* did not reach the end of a regular strip */
//						while( startElement + 1 < noe && guVec3iGetLocInd( flatElement[ startElement+1 ], startVertex ) != -1 )
//							startElement++;
//				}
//				
//				return self;
//			}
//
//
//			- computeLineFromVertex: (int)s inElement: (int)elem toVertex: (int)e appendTo: (OoPolygonOnTriang *)shortest  {
//
//				int index, nop = [shortest getNumVertices];
//				
//				OsVec3i *flatElement;
//				OsVec3 *flatVertex, *baryData;
//				OsVec2 v, w;
//
//				int *elementData;
//
//				double det, *S, *E, *L, *R;
//
//				flatVertex  = [ vertex getData ];
//				flatElement = [ element getData ];
//
//				gmSubVec2( v, S = flatVertex[ s ], E = flatVertex[ e ] );
//
//				ASSERT( (index = guVec3iGetLocInd( flatElement[ elem ], s ))!=-1, OE_MISC, OE_DEFAULT,
//						"start Vertex is not in element", return nil );
//				
//				[ shortest setNumVertices: nop + e - s + 1];
//				
//				baryData    = [ [ shortest getBaryCoords ] getData ];
//				elementData = [ [ shortest getElements ] getData ];
//
//				if( !nop ) {
//					/* init the start point only if it is the first*/
//					elementData[nop] = elem;
//					gmZeroVec3( baryData[nop] ); baryData[nop++][ index ] = 1.; 
//				}
//
//				/* compute the intersections with the edges */
//
//				while( (index=guVec3iGetLocInd( flatElement[ elem ], e )) == -1 ) {
//
//					gmSubVec2( w, R = flatVertex[ flatElement[ elem ][2] ], L = flatVertex[  flatElement[ elem ][1] ] ); 
//
//					det = v[0]*w[1] - v[1]*w[0];
//
//					baryData[nop][0] = 0;
//					if( det == 0 ) {
//						/* some lines do not intersect, set standart coordinates */
//						baryData[nop][1] = baryData[nop][2] = .5;
//					} else {
//						/* compute barycentric coordinates for the intersection point */
//						baryData[nop][1] = ( ( R[1]-E[1] ) * v[0] - (R[0]-E[0]) * v[1] ) / det;
//						baryData[nop][2] = 1 - baryData[nop][1];
//
//						if( !guBaryIsInsideElement( baryData[nop] ) ) {
//							ASSERT( 0, OE_MISC, OE_DEFAULT, "intersection is not on edge", printf(" %f at Node %d \n", baryData[nop][1], nop  ) );
//
//							baryData[nop][1] = OO_MAX( 0., OO_MIN( 1. ,baryData[nop][1] ) ); 
//							baryData[nop][2] = 1 - baryData[nop][1];
//						}
//					}
//					elementData[nop++] = elem++;
//				}
//
//				/* init the end point */
//				elementData[nop] = elem;
//				gmZeroVec3( baryData[nop] ); baryData[nop++][ index ] = 1.; 
//
//				[shortest setNumVertices: nop];
//				[shortest setTriang: self];
//
//				return self;
//				
//			}
//			
	
}
