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
public class IntTriple {

	IntTriple() {
		super();
	}
	

	public static int getOppEntry( int [] triple, int entryA, int entryB ) {

		if( triple[0] != entryA && triple[0] != entryB ) return( triple[0] );
		if( triple[1] != entryA && triple[1] != entryB ) return( triple[1] );
		return( triple[2] );
	}

	public static int getOppLocInd( int [] triple, int entryA, int entryB ) {

		if( triple[0] != entryA && triple[0] != entryB ) return( 0 );
		if( triple[1] != entryA && triple[1] != entryB ) return( 1 );
		if( triple[2] != entryA && triple[2] != entryB ) return( 2 );

		return( -1 );
	}

	public static int getLocInd( int [] triple, int entry ) {
		
		if( triple[0] == entry ) return( 0 );
		if( triple[1] == entry ) return( 1 );
		if( triple[2] == entry ) return( 2 );

		return( -1 );
	}

	public static void replaceEntry( int [] triple, int entry, int newEntry ) {

		if( triple[0] == entry ) triple[0] = newEntry;
		if( triple[1] == entry ) triple[1] = newEntry;
		if( triple[2] == entry ) triple[2] = newEntry;
	}

	public static int getMax( int [] triple ) {  
		int max = triple[0];
		if( triple[1] > max )
			max = triple[1];
		if( triple[2] > max )
			max = triple[2];
		return max;
	}

	public static int getMin( int [] triple ) {  
		int min = triple[0];
		if( triple[1] < min )
			min = triple[1];
		if( triple[2] < min )
			min = triple[2];
		return min;
	}

	public static int getMap( int [] map, int [] triple1, int [] triple2 ) {
		int i, j, k = 0;
		for( i=0; i<3; i++ ) {
			map[i] = -1;
			for( j=0;j<3; j++ ) 
				if(  triple1[j] == triple2[i] ) { 
					map[i]=j; /* result: triple1[ map[i] ] == triple2[ i ] */
					k++;
					break;
				}
		}
		return k;
	}

	public static void applyToVec3( int [] map, double [] vec ) {
		double [] h =new double[3];
		h[0]=vec[0]; h[1]=vec[1]; h[2]=vec[2];
		if( map[0]<0 ) 
			vec[0] = 0;
		else vec[0] = h[map[0]];
		if( map[1]<0 )
			vec[1] = 0;
		else vec[1] = h[map[1]];
		if( map[2]<0 )
			vec[2] = 0;
		else vec[2] = h[map[2]];
	}

	public static void applyInvToVec3( int [] map, double[] vec ) {
		double []  h = new double[3];
		h[0]=vec[0]; h[1]=vec[1]; h[2]=vec[2];
		if( map[0]<0 ) 
			vec[0] = 0;
		else  vec[map[0]] = h[0];
		if( map[1]<0 )
			vec[1] = 0;
		else vec[map[1]] = h[1];
		if( map[2]<0 )
			vec[2] = 0;
		else vec[map[2]] = h[2];
	}



	public static double getMapForVec3( int [] map, double [] vec1, double [] vec2 ) {
		int i, j;
		double min, help, out = 0.;
		boolean [] is = new boolean[3];

		for( i=0; i<3; i++ ) {
			min = 1e99;
			for( j=0;j<3; j++ ) 
				if( !is[j] && (help=Math.abs( vec1[i]-vec2[j] )) < min ) { 	
					map[i]=j;
					min = help;
				}
			is[ map[ i ] ] = true;
			out += min;    
		}
		return out / (Math.abs(vec1[0]) + Math.abs(vec1[1]) + Math.abs(vec1[2] ) );
	}

	

}
