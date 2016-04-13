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

public class DiscretizedInterval implements Serializable, Cloneable {

    public static final long serialVersionUID = 1L;

    final double [] min;
    final double [] max;
    
    final int [] discr;

    final int dim;

    final int numOfVertices;

    static double [] constant( double value, int dim ) {
	double [] constant = new double[ dim ];
       
	for( int i=0; i<dim; i++ )
	    constant[i] = value;

	return constant;
    }

    static int [] constant( int value, int dim ) {
	int [] constant = new int[ dim ];
       
	for( int i=0; i<dim; i++ )
	    constant[i] = value;

	return constant;
    }

    public DiscretizedInterval( Interval interval, int [] discr ) {
	this( interval.min, interval.max, discr );
    }

    public DiscretizedInterval( Interval interval, int discr ) {
	this( interval.min, interval.max, constant( discr, interval.getDim() ) );
    }

    public DiscretizedInterval( int [] discr ) {
	this( constant( 1., discr.length ), discr );
    }

    public DiscretizedInterval() {
	this( 1, 1 );
    }

    public DiscretizedInterval( int dim, int discr ) {
	this( constant( 1., dim ), constant( discr, dim ) );
    }

    public DiscretizedInterval( double [] max, int [] discr ) {

	this( new double[ max.length ], max, discr );
    }

    public DiscretizedInterval( double [] min, double [] max, int [] discr ) {
	
	if( min.length != max.length || min.length != discr.length )
	    throw new IllegalArgumentException( "length of arrays do not coincide" );
	
	this.min   = (double[])min.clone();
	this.max   = (double[])max.clone();

	this.discr = (int[])discr.clone();

	dim = min.length;

	int vol = discr[0] + 1;
	for( int i=1; i<dim; i++ )
	    vol *= discr[i] + 1;	
	numOfVertices = vol;

    }

    public Interval getInterval() {
	return new Interval( min, max );
    }

    public void setInterval( Interval interval ) {
	if( interval.getDim() != dim )
	    throw new IllegalArgumentException( "interval has wrong dimension" );

	System.arraycopy( interval.min, 0, min, 0, dim );
	System.arraycopy( interval.max, 0, max, 0, dim );
    }


    public void getValue( int [] index, double [] value ) {

	for( int i=0; i<dim; i++ ) {
	    final int I = index[i];
	    final int D = discr[i];

	    if( I>D || I<-1 )
		throw new IndexOutOfBoundsException( "multi index out of bounds" );
	    
	    final double t = I / (double)D ;
	    
	    value[i] = (1-t) * min[i] + t * max[i];
	}
    }

    public double [] getValue( int [] index ) {	
	double [] value = new double[ dim ];
	getValue( index, value );
	return value;
    }

    public int getNumOfVertices() {
	return numOfVertices;
    }

    public int getDim() {
	return dim;
    }
}






