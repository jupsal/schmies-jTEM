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

package de.jtem.mfc.set;


public class MultiIndexedDoubleSet extends MultiIndexedSet {

    public static final long serialVersionUID = 1L;

    final double [] set;

    public MultiIndexedDoubleSet( int [] max ) {
	this( new int[ max.length ], max );
    }

    public MultiIndexedDoubleSet( MultiIndexedSet multiIndexedSet ) {
	this( multiIndexedSet.min, multiIndexedSet.max );
    }

    public MultiIndexedDoubleSet( MultiIndexedDoubleSet multiIndexedDoubleSet ) {
	super( multiIndexedDoubleSet.min, multiIndexedDoubleSet.max );

	set = (double [])multiIndexedDoubleSet.set.clone();
    }

    public MultiIndexedDoubleSet( int [] min, int [] max ) {

	super( min, max );

	set = new double[ numOfVertices ];
    }

    public double getValue( int [] index ) {
	return getValue( index, true );
    }

    public double getValue( int [] index, boolean checkBounds ) {
	return set[ getLinearIndexOf( index, checkBounds ) ];
    }

    public double getValue( int linearIndex ) {
	return set[ linearIndex ];
    }


    public void setValue( int [] index, double value ) {
	setValue( index, value, true );
    }

    public void setValue( int [] index, double value, boolean checkBounds ) {
	set[ getLinearIndexOf( index, checkBounds ) ] = value;
    }

    public void setValue( int linearIndex, double value ) {
	set[ linearIndex ] = value;
    }


    public Object toArray() {
	return toArray( (Object)null );
    }

    public Object toArray( Object array ) {
	
	switch( dim ) {
	    	    
	case 1:
	    return toArray( (double [])array );

	case 2:
	    return toArray( (double [][])array );

	case 3:
	    return toArray( (double [][][])array );

	case 4:
	    return toArray( (double [][][][])array );

	default:
	    throw new UnsupportedOperationException("dimension bigger then 4 are not supported");
	}	
    }

    public double [] toArray( double [] array ) {

	if( dim != 1 )
	    throw new IllegalArgumentException( "array has wrong dimension" );

	if( array == null || array.length != max[0] - min[0] + 1 )
	    array = new double [ max[0] - min[0] + 1 ];

	System.arraycopy( set, 0, array, 0, array.length );

	return array;
    }

    public double [][] toArray( double [][] array ) {

	if( dim != 2 )
	    throw new IllegalArgumentException( "array has wrong dimension" );

	final int xSize = max[0] - min[0] + 1;
	final int ySize = max[1] - min[1] + 1;

	if( array == null || 
	    array   .length != xSize ||
	    array[0].length != ySize )
	    array = new double [ xSize ][ ySize ];

	
	for( int i=0, j=0; i<xSize; i++, j+=ySize )
	    System.arraycopy( set, j, array[i], 0, ySize );

	return array;
    }

    public double [][][] toArray( double [][][] array ) {

	if( dim != 3 )
	    throw new IllegalArgumentException( "array has wrong dimension" );

	final int xSize = max[0] - min[0] + 1;
	final int ySize = max[1] - min[1] + 1;
	final int vSize = max[2] - min[2] + 1;

	if( array == null || 
	    array      .length != xSize ||
	    array[0]   .length != ySize ||
	    array[0][0].length != vSize )
	    array = new double [ xSize ][ ySize ][ vSize ];

	
	for( int i=0, k=0; i<xSize; i++ )
	    for( int j=0; j<ySize; j++, k+=vSize )
		System.arraycopy( set, k, array[i][j], 0, vSize );

	return array;
    }


    public double [][][][] toArray( double [][][][] array ) {

	if( dim != 4 )
	    throw new IllegalArgumentException( "array has wrong dimension" );

	final int xSize = max[0] - min[0] + 1;
	final int ySize = max[1] - min[1] + 1;
	final int vSize = max[2] - min[2] + 1;
	final int wSize = max[3] - min[3] + 1;

	if( array == null || 
	    array         .length != xSize ||
	    array[0]      .length != ySize ||
	    array[0][0]   .length != vSize ||
	    array[0][0][0].length != wSize )
	    array = new double [ xSize ][ ySize ][ vSize ][ wSize ];

	
	for( int i=0, l=0; i<xSize; i++ )
	    for( int j=0; j<ySize; j++ )
		for( int k=0; k<vSize; k++, l+=wSize )
		    System.arraycopy( set, l, array[i][j][k], 0, wSize );

	return array;
    }


    
}





