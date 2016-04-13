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

import java.io.Serializable;

public class MultiIndexedSet implements Serializable, Cloneable {

    public static final long serialVersionUID = 1L;

    final int [] min;
    final int [] max;

    final int [] size;

    final int dim;

    final int numOfVertices;

    public MultiIndexedSet( int [] max ) {
	this( new int[ max.length ], max );
    }

    public MultiIndexedSet( MultiIndexedSet multiIndexedSet ) {
	this( multiIndexedSet.min, multiIndexedSet.max );
    }

    public MultiIndexedSet( int [] min, int [] max ) {
	
	if( min.length != max.length )
	    throw new IllegalArgumentException( "length of min and max do not coincide" );
	
	this.min = (int[])min.clone();
	this.max = (int[])max.clone();

	dim = min.length;

	size = new int[ dim ];

	size[dim-1] = 1;
	for( int i=dim-1; i>0; i-- )
	    size[i-1] = (max[i] - min[i] + 1) * size[i];	

	int vol = max[0] - min[0] + 1;
	for( int i=1; i<dim; i++ )
	    vol *= max[i] - min[i] + 1;	
	numOfVertices = vol;
    }

    public int getNumOfVertices() {
	return numOfVertices;
    }

    public int getDim() {
	return dim;
    }    

    public boolean isOutOfBounds( int [] index ) {
	if( index.length != dim )
	    return true;

	for( int i=0; i<dim; i++ )
	    if( index[i] < min[i] || index[i] > max[i] )
		return true;

	return false;
    }

    public int getLinearIndexOf( int [] index ) {
	return getLinearIndexOf( index, true );
    }
    
    public int getLinearIndexOf( int [] index, boolean checkBounds ) {
	
	if( checkBounds && isOutOfBounds( index ) )
	    throw new IndexOutOfBoundsException("multi index is out of bounts" );
	
	int linearIndex = index[dim-1];
	
	for( int i=dim-2; i>=0; i-- )
	    linearIndex += size[i] * index[i];
	
	return linearIndex;
    }

    public boolean equals( MultiIndexedSet interval ) {

	if( this == interval )
	    return true;

	if( min.length != interval.min.length )
	    return false;

	for( int i=0; i<min.length; i++ ) {
	    if( min[i] != interval.min[i] || max[i] != interval.max[i] )
		return false;
	}
	
	return true;
    }

    public boolean equals( Object object ) {

	try {
	    return equals( (MultiIndexedSet)object );
	} catch( ClassCastException e ) {
	    return false;
	}
    }


}





