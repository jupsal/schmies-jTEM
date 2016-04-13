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


public class DiscretizedIntervalIterator extends DiscretizedInterval {

    public static final long serialVersionUID = 1L;

    final double [] value;
    final double [] delta;

    final int [] index;

    final int maxLinearIndex;

    int linearIndex = 0;

    boolean first = true;

    public DiscretizedIntervalIterator( Interval interval, int [] discr ) {
	this( interval.min, interval.max, discr );
    }

    public DiscretizedIntervalIterator( Interval interval, int discr ) {
	this( interval.min, interval.max, constant( discr, interval.getDim() ) );
    }

    public DiscretizedIntervalIterator( int [] discr ) {
	this( constant( 1., discr.length ), discr );
    }

    public DiscretizedIntervalIterator() {
	this( 1, 1 );
    }

    public DiscretizedIntervalIterator( int dim, int discr ) {
	this( constant( 1., dim ), constant( discr, dim ) );
    }

    public DiscretizedIntervalIterator( double [] max, int [] discr ) {

	this( new double[ max.length ], max, discr );
    }

    public DiscretizedIntervalIterator( double [] min, double [] max, int [] discr ) {

	super( min, max, discr );

	this.value = (double[])min.clone();

	index = new int[dim];	

	delta = new double[dim];

	for( int i=0; i<dim; i++ )
	    delta[i] = ( max[i] - min[i] ) / discr[i];

	maxLinearIndex = numOfVertices-1;
    }

    public double[] getValue() {
	return value;
    }

    public int [] getIndex() {
	return index;
    }

    public int getLinearIndex() {
	return linearIndex;
    }

    public boolean hasNext() {
	return linearIndex < maxLinearIndex;
    }

    public void restartIteration() {
	first = true;
	linearIndex = 0;
    }

    public Object next() {

	if( !(linearIndex < maxLinearIndex) )
	    throw new java.util.NoSuchElementException();

	if( !first )
	    iterate();
	else
	    first = false;

	return value.clone();
    }

    public void iterate() {

	if( !(linearIndex < maxLinearIndex) )
	    throw new java.util.NoSuchElementException();

	linearIndex++;
	first = false;

	for( int i=dim-1; i> 0; i-- ) {
	    index[i]++;	    
	    if( index[i] < discr[i] ) {
		value[i] += delta[i];
		return;
	    } else if( index[i] <= discr[i] ) { 
		value[i] = max[i];
		return;
	    }

	    index[i]=0;
	    value[i]=min[i];
	}

	index[0]++;
	value[0]+=delta[0];
    }

    public void remove() {
	throw new RuntimeException( "UnsupportedOperationException" ); //java 1 compatible
    }
}






