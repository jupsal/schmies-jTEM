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

package de.jtem.numericalMethods.geometry.mesh;

import java.util.Collection;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Iterator;
import java.util.LinkedList;
import java.util.List;
import java.util.Set;
import java.util.Vector;

/**
 * This class represents an unoriented edge.
 * 
 * It also provides utility methods for indexed face sets, e.g.
 * computing the connectivity given by a winged edge structure or
 * computing the bounary curves as indexed line set.
 * 
 * @author schmies
 */
public class Edge {

	final static int MAX_VERTEX_INDEX = (int)Math.sqrt( Integer.MAX_VALUE - 2 );
	
	/**
	 * Vertex indices forming the edge.
	 */
	public final int v1, v2;
	
	/**
	 * A winged edge represents an edge in a 2-manifold.
	 * It either can  have one or two ajacent faces.
	 * 
	 * @author schmies
	 */
	public static class Winged extends Edge {
		
		/**
		 * Constant indicating that no face is defined.
		 */
		public final static int NO_FACE = -1;
		
		/**
		 * Indices of adjacent faces of edge.
		 */
		public final int f1, f2;
			
		public Winged( int v1, int v2, int f1 ) {
			this( v1, v2, f1, NO_FACE );
		}
		
		public Winged( int v1, int v2, int f1, int f2 ) {
			super( v1, v2 );
			this.f1 = f1;
			this.f2 = f2;
		}
		
		public String toString() {
			return "(v1="+v1+",v2="+v2 +",f1="+f1 +",f2="+f2 + ")";
		}
	}
	
	/**
	 * Creates edge with given vertex indedices v1 and v2.
	 * @param v1 vertex index
	 * @param v2 vertex index
	 */
	public Edge( int v1, int v2 ) {
		this.v1 = v1;
		this.v2 = v2;
	}
	
	static int hashCode( int v1, int v2 ) {
		return v1 > v2 ? v1*MAX_VERTEX_INDEX + v2 : v2*MAX_VERTEX_INDEX + v1;
	}
	
	final public int hashCode() {
		return hashCode( v1, v2 );
	}
	
	public boolean equals( Object o ) {
		if( !( o instanceof Edge ) )
			return false;
		
		Edge we = (Edge)o;
		
		return v1 == we.v1 && v2 == we.v2 || v1 == we.v2 && v2 == we.v1;
	}
	
	/**
	 * Checks whether vertex with given index is conatained in edge.
	 * @param index of vertex
	 * @return true if vertex is contained in edge, false otherwise.
	 */
	public final boolean containsVertex( int index ) {
		return v1==index || v2==index;
	}
	
	/**
	 * Returns common vertex index of two edges.
	 * If both indices coinside v1 is returned.
	 * @throws IllegalArguementException if edges do not share common vertex
	 * @param we other edege
	 * @return common vertex
	 */
	public final int getCommonVertex( Edge we ) {
		if( we.containsVertex(v1) ) {
			return v1;
		} else if( we.containsVertex(v2)) {
			return v2;
		} else
			throw new IllegalArgumentException( "edges do not share common vertex" );
	}
	
	/**
	 * Return other vertex then the one given by index. 
	 * @throws IllegalArguementException if the given index is not a vertex of the edge
	 * @param index of vertex
	 * @return other vertex
	 */
	int getOtherVertex( int index ) {
		if( v1 == index ) {
			return v2;
		} else if( v2 == index ) {
			return v1;
		} else
			throw new IllegalArgumentException( "given index is not a vertex of the edge" );
	}
	
	/**
	 * Computes winged edge representation of presribed indexed face set.
	 * Inner  and boundary edges are given as seperate hash maps.
	 * @param ifs indexed face set
	 * @param innerEdges Set of inner edges (winged) on ouput
	 * @param boundaryEdges Set of boundary edges (winged)  on output
	 */
	public static void computeEdgesOfIndexedFaceSet(  int [][] ifs, Set innerEdges, Set boundaryEdges ) {
	
		HashMap hm = new HashMap();
		
		if( innerEdges != null ) innerEdges.clear();
		else innerEdges = new HashSet();
		
		if( boundaryEdges != null ) boundaryEdges.clear();
		else boundaryEdges = new HashSet();
		
		
		final int nof = ifs.length;
		
		for( int i=0; i<nof; i++ ) {
			final int [] face = ifs[i];
			
			for( int j=1; j<face.length; j++ ) {
				addEdge(innerEdges, hm, i, face[j-1], face[j]);			
			}
			addEdge(innerEdges, hm, i, face[face.length-1], face[0] );
		}
		
		boundaryEdges.addAll( hm.values() );
	}

	/**
	 * @param innerEdges
	 * @param hm
	 * @param face
	 * @param k vertex index
	 * @param l vertex index
	 */
	private static void addEdge( Set innerEdges, HashMap hm, int face, int k, int l) {
		Edge we1 = new Edge( k, l );
		
		Edge.Winged we2 = (Edge.Winged)hm.get(we1);
		if( we2 != null ) {
			hm.remove(we1);
			innerEdges.add( new Edge.Winged( we2.v1, we2.v2, we2.f1, face));
		} else {
			hm.put( we1, new Edge.Winged( k,l, face) );
		}
	}

	/**
	 * Computes winged edge representation of presribed indexed face set.
	 * @param ifs indexed line set
	 * @return Set of edges (winged)
	 */
	static public HashSet getWingedEdgesOfIndexedFaceSet( int [][] ifs ) {
		HashSet innerEdges  = new HashSet();
		HashSet boundaryEdges = new HashSet(); 
		
		computeEdgesOfIndexedFaceSet(  ifs, innerEdges, boundaryEdges );
		
		innerEdges.addAll( boundaryEdges );
		
		return innerEdges;
	}
	
	/**
	 * Computes inner edges of presribed indexed face set.
	 * @param ifs indexed line set
	 * @return Set of inner edges (winged)
	 */
	static public HashSet getInnerEdgesOfIndexedFaceSet( int [][] ifs ) {
		HashSet innerEdges  = new HashSet();
		HashSet boundaryEdges = new HashSet(); 
		
		computeEdgesOfIndexedFaceSet(  ifs, innerEdges, boundaryEdges );
		
		return innerEdges;
	}
	
	/**
	 * Computes boundary edges of presribed indexed face set.
	 * @param ifs indexed line set
	 * @return Set of boundary edges (winged)
	 */
	static public HashSet getBoundaryEdgesOfIndexedFaceSet( int [][] ifs ) {
		HashSet innerEdges  = new HashSet();
		HashSet boundaryEdges = new HashSet(); 
		
		computeEdgesOfIndexedFaceSet(  ifs, innerEdges, boundaryEdges );
		
		return boundaryEdges;
	}
	
	/**
	 * Returns boundary of face set as line set.
	 * @param ifs indexed face set
	 * @return boundary of face set as line set
	 */
	static public int [][] getBoundaryOfIndexedFaceSet( int [][] ifs ) {
		return stripEdges( getBoundaryEdgesOfIndexedFaceSet( ifs ));
	}
	
	/**
	 * Returns indexed line set of 
	 * @param hs set of edges, which will be destroyed on ouput
	 */
	public static int [][] stripEdges( Set hs ) {
		Vector lines = new Vector();
		stripEdges( hs, lines );
		
		
		int [][] ils = new int[lines.size()][];
		for(int i=0; i<ils.length; i++) {
			List line = (List)lines.get(i);
			ils[i] = new int[line.size()+1];
			Edge firstEdge = (Edge)line.get(0);
			if( ils[i].length  == 2 ) {
				ils[i][0] = firstEdge.v1;
				ils[i][1] = firstEdge.v2;
			}
			
			ils[i][1] = firstEdge.getCommonVertex( (Edge)line.get(1));
			ils[i][0] = firstEdge.getOtherVertex(ils[i][1]);
			for( int j=2; j<ils[i].length; j++ ) {
				ils[i][j] = ((Edge)line.get(j-1)).getOtherVertex(ils[i][j-1]);
			}		
		}
	
		return ils;
	}
	
	/**
	 * Returns indexed line set giving the boundary of the presribed 
	 * indexed face set.
	 * @param ifs indexed face set.
	 * @return indexed line set of boundary curves.
	 */
	static public int [][] getBoundaryLinesOfIndexedFaceSet( int [][] ifs ) {
		return stripEdges(getBoundaryEdgesOfIndexedFaceSet( ifs ));
	}
	
	/**
	 * Strips set of edges to collections of edges giving line segments.
	 * @param edge set, which will be destroyed on ouput
	 * @param lines collection of line segements given by lists of edges (on output)
	 */
	static public void stripEdges( Set edges, Collection lines ) {
		
		 lines.clear(); 
		 
		Vector usedEdges = new Vector();
		
		while( !edges.isEmpty() ) {
	
			Iterator iter = edges.iterator();
			Edge we1 = (Edge)iter.next();
			
			int firstIndex = we1.v1;
			int lastIndex = we1.v2;
			
			LinkedList line   = new LinkedList();
			
			usedEdges.clear();
			usedEdges.add(we1);
			
			line.add( we1);
			
			while( true ) {
				while( iter.hasNext() ) {
					Edge we2 = (Edge)iter.next();
					
					if( we2.containsVertex(firstIndex)) {
						firstIndex = we2.getOtherVertex(firstIndex);
						line.addFirst(we2);
						usedEdges.add(we2);
					} else if( we2.containsVertex(lastIndex)) {
						lastIndex = we2.getOtherVertex(lastIndex);
						line.addLast(we2);
						usedEdges.add(we2);
					} 
				}
				
				edges.removeAll(usedEdges);
				
				if( usedEdges.isEmpty() )
					break;
				
				usedEdges.clear();
				iter = edges.iterator();
			}
			
			lines.add(line);	
		}
	}
	
	public String toString() {
		return "(v1="+v1+",v2="+v2  + ")";
	}
	
	public static void main( String [] arg ) {
	
		int [][] ifs = new int[][] {
				{ 0, 5, 4, 1},
				{ 1, 4, 3 },
				{ 1, 3, 2 },
		};
		
		System.out.println( getBoundaryEdgesOfIndexedFaceSet(ifs));
		
		int [][] ils = stripEdges( getBoundaryEdgesOfIndexedFaceSet( ifs ) );
		
		for( int i=0; i<ils.length; i++ ) {
			for( int j=0; j<ils[i].length; j++ ) {
				System.out.print( "->" + ils[i][j] );
			}
			System.out.println(  );
		}
	}
	
}
