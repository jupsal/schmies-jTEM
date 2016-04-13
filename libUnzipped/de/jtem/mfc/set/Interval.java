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
import java.io.StreamTokenizer;
import java.io.StringReader;


public class Interval implements Serializable, Cloneable {

    private static final long serialVersionUID = 1L;

    private static final double EPS = 1e-12;

    double [] min;
    double [] max;

    private final void checkMinMax() {
	if( min.length != max.length )
	    throw new IllegalArgumentException( "length of min and max do not coincide" );

	for( int i=0; i<max.length; i++ )
	    if( min[i] > max[i] )
		throw new IllegalArgumentException( "min is bigger max" );
    }

    public Interval() {
	this( 1 );
    }

    public Interval( double a, double b ) {
	this( 1 );
	min[0] = a;
	max[0] = b;
	checkMinMax();
    }
    
    public Interval( int dim ) {
	
	min = new double[dim];
	max = new double[dim];
	
	for( int i=0; i<max.length; max[i++]=1 );
    }

    public Interval( Interval interval ) {
	min = (double [])interval.min.clone();
	max = (double [])interval.max.clone();	
    }

    public Interval( String s ) {
	assign( s );
    }

    public Interval( double [] min, double max [] ) {
	this.min = (double [])min.clone();
	this.max = (double [])max.clone();

	checkMinMax();
    }

    public int getDim() {
	return min.length;
	
    }

    public double getMin( int index ) {
	return min[index];
    }
    
    public void setMin( int index, double min ) {
	this.min[index] = min;
    }

    public double getMax( int index ) {
	return max[index];
    }

    public void setMax( int index, double max ) {
	this.max[index] = max;
    }

    public double getSize( int index ) {
	return max[index] - min[index];
    }

    public double getVolume() {
	double volume = max[0] - min[0];
	for( int i=1; i<min.length; i++ )
	    volume *=  max[i] - min[i];
	return volume;
    }


    public double getDiameter() {
	double lengthSqr = ( max[0] - min[0] ) * ( max[0] - min[0] );
	for( int i=1; i<min.length; i++ )
	    lengthSqr +=  ( max[i] - min[i] ) * ( max[i] - min[i] );
	return Math.sqrt( lengthSqr );
    }

    public void assignTensor( Interval interval1, Interval interval2 ) {

	if( min.length != interval1.min.length + interval2.min.length ) {
	    min = new double[ interval1.min.length + interval2.min.length ]; 
	    max = new double[ interval1.min.length + interval2.min.length ]; 
	}

	System.arraycopy( interval1.min, 0, min, 0,                    interval1.min.length );
	System.arraycopy( interval2.min, 0, min, interval1.min.length, interval2.min.length );

	System.arraycopy( interval1.max, 0, max, 0,                    interval1.max.length );
	System.arraycopy( interval2.max, 0, max, interval1.max.length, interval2.max.length );
    }

    public Interval tensor( Interval interval ) {
	Interval tensor = new Interval( min.length + interval.min.length );
	tensor.assignTensor( this, interval );
	return tensor;
    }

    public void assignBound( Interval interval1, Interval interval2 ) {
	if( interval1.min.length != interval2.min.length )
	    throw new IllegalArgumentException( "dimensions of intervals do not coincide" );

	if( min.length != interval1.min.length ) {
	    min = new double[ interval1.min.length ]; 
	    max = new double[ interval1.min.length ];
	}

	for( int i=0; i<min.length; i++ ) {
	    min[i] = Math.min( interval1.min[i], interval2.min[i] );
	    max[i] = Math.max( interval1.max[i], interval2.max[i] );
	}
    }

    public void assignBound( Interval interval ) {
	this.assignBound( this, interval );
    }

    public Interval bound( Interval interval ) {
	Interval bound = new Interval( min.length + interval.min.length );
	bound.assignBound( this, interval );
	return bound;
    }

    public void assignJoin( Interval interval1, Interval interval2 ) {

	if( !interval1.joinable( interval2 ) )
	    throw new ArithmeticException( "intervals are not joinable in a single interval" );

	assignBound( interval1, interval2 );
    }
    
    public void assignJoin( Interval interval ) {
	assignJoin( this, interval );
    }

    public Interval join( Interval interval ) {
	Interval join = new Interval( min.length + interval.min.length );
	join.assignJoin( this, interval );
	return join;
    } 

    public boolean joinable( Interval interval) {
	return joinable( interval, EPS );
    }

    public boolean joinable( Interval interval, double eps ) {
	if( min.length != interval.min.length )
	    return false;

	boolean wasAlreadyDifferent = false;

	for( int i=0; i<min.length; i++ ) {

	    if( Math.abs( min[i] - interval.min[i] ) > eps || 
		Math.abs( max[i] - interval.max[i] ) > eps    ) {

		if( wasAlreadyDifferent )
		    return false;
		
		if( !( min[i] - EPS < interval.min[i] && interval.min[i] < max[i] + EPS ||
		       min[i] - EPS < interval.max[i] && interval.max[i] < max[i] + EPS    ) ) {
		    return false;
		}

		wasAlreadyDifferent = true;
	    }
	}
	
	return true;
    }


    public boolean intersect( Interval interval) {
	return intersect( interval, EPS );
    }

    public boolean intersect( Interval interval, double eps ) {
	if( min.length != interval.min.length )
	    return false;

	for( int i=0; i<min.length; i++ ) {
		
	    if( min[i] - EPS > interval.max[i] || max[i] + EPS < interval.min[i] )
		return false;
	}
	
	return true;
    }

    public void assign( String s ) {

	double [] array = stringToArray( s, (double [])null );

	if( array.length % 2 != 0 )
	    throw new IllegalArgumentException();

	if( min == null || min.length != array.length / 2 ) {
	    min = new double[ array.length / 2 ];
	    max = new double[ array.length / 2 ];
	}

	for( int i=0, k=0; i<min.length; i++ ) {
	    min[i] = array[k++];
	    max[i] = array[k++];
	}
    }

    public void assign( Interval interval ) {
	assign( interval.min, interval.max );
    }

    public void assign( double [] min, double [] max ) {

	if( min.length != max.length )
	    throw new IllegalArgumentException( "length of min and max do not coincide" );
	
	if( min.length != this.min.length ){
	    this.min = (double [])min.clone();
	    this.max = (double [])max.clone();
	} else {
	    System.arraycopy( min, 0, this.min, 0, min.length );
	    System.arraycopy( max, 0, this.max, 0, max.length );
	}
    }

    public void assignByCenter( double [] center, double delta [] ) {

	if( center.length != delta.length )
	    throw new IllegalArgumentException( "length of arrays do not coincide" );

	assign( center, center );
	
	for( int i=0; i<center.length; i++ ) {
	    min[i] -= Math.abs( delta[i] );
	    max[i] += Math.abs( delta[i] );
	}
    }

    public void assignByCenter( double [] center, double relDelta ) {
	assign( center, center );

	for( int i=0; i<center.length; i++ ) {
	    min[i] -= Math.abs( center[i] * relDelta );
	    max[i] += Math.abs( center[i] * relDelta );
	}
    }

    public void getCenter( double [] center, int offset ) {
	
	for( int i=0; i<min.length; i++ )
	    center[ offset++ ] = ( min[i] + max[i] ) / 2;	
    }

    public double [] getCenter() {
	double [] center = new double[ min.length ];
	getCenter( center, 0 );
	return center;
    }

    public boolean equals( Interval interval, double eps ) {

	if( this == interval )
	    return true;

	if( min.length != interval.min.length )
	    return false;

	for( int i=0; i<min.length; i++ ) {
	    if( Math.abs( min[i] - interval.min[i] ) > eps )
		return false;
	    if( Math.abs( max[i] - interval.max[i] ) > eps )
		return false;
	}
	
	return true;
    }

    public boolean equals( Interval interval ) {
	return equals( interval, EPS );
    }

    public boolean equals( Object object ) {

	try {
	    return equals( (Interval)object, EPS );
	} catch( ClassCastException e ) {
	    return false;
	}
    }

    public String toString() {
	StringBuffer sb=new StringBuffer(100);
	sb.append("[");
	sb.append( min[0] );
	sb.append(",");
	sb.append( max[0] );
	for( int i=1; i<min.length; i++ ) {
	    sb.append( "]x[" );
	    sb.append( min[i] );
	    sb.append(",");
	    sb.append( max[i] );
	}
	sb.append("]");

	return sb.toString();
    }

    static double[] stringToArray( String s, double [] basket ) {

	if( basket == null || basket.length < s.length() / 2 + 1 )
	    basket = new double [s.length() / 2 + 1];
	
	try{
	    
	    StreamTokenizer st = new StreamTokenizer( new StringReader(s) );
	    
	    st.resetSyntax();	    
	    st.whitespaceChars('\u0000', '\u002C'); // ',', ';', ':',  '(', ')'
	    st.whitespaceChars('[', ']'); 
	    st.whitespaceChars('x', 'x'); 

	    st.parseNumbers();

	    int num = 0;
      
	    while( st.nextToken() != StreamTokenizer.TT_EOF ) {

		basket[num++] = st.nval;
	    }

	    final double [] array = new double[num];

	    System.arraycopy( basket, 0, array,0, num);

	    return array;

	}

	catch(java.io.IOException ex){ throw new Error(); }
    }
    
}





