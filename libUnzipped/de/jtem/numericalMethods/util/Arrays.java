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

package de.jtem.numericalMethods.util;


/**
 * @depricated
 This class provides utilities for arrays as static methods.
    <p />
    Most operation performed on the <i>target</i> array can imply
    a change of its length. Those methods create a new array and
    return this instead of the original which is only returned
    if it has already the right length.
    @author	Markus Schmies
 */

public final class Arrays {

    private Arrays() {}

    /** Just returns <i>target</i> if its length equals <i>size</i>.
	Otherwise, this includes <i>target</i> is null, it creates
	an array of length <i>size</i> and returns this instead. */
    public static final double [] size( final double [] target, final int size ) {
	return target.length == size ? target : new double[size];
    }

    /** Like {@link Arrays#size(double[],int)} but copies also the contents of
	the <i>target</i> if nessecary. */
    public static final double [] resize( final double [] target, final int size ) {
	if( target.length == size )
	    return target;

	double [] newTarget = new double [size];

	System.arraycopy( target, 0, newTarget, 0, Math.min( size, target.length ) );

	return newTarget;
    }

    /** Returns an array containing <i>numOfTimes</i> copies of the
	contens of <i>values</i>.
	The method uses <i>target</i> if it is not <code>null</code> and has the right size. */
    public static final double [] unrole( double [] target,
					  final double [] values, final int numOfTimes ) {
	int length = values.length * numOfTimes;

	if( target.length != values.length * numOfTimes )
	    target = new double[ length ];

	System.arraycopy( values, 0, target, 0, values.length );

	int i=values.length;
	for( ; 2*i < length; i *= 2 )
	    System.arraycopy( target, 0, target, i, i );

	System.arraycopy( target, 0, target, i, length - i );

	return target;
    }

    /** Retuns an array containing the contents of <i>con</i> and
	<i>cat</i>.
	The method uses <i>target</i> if it is not <code>null</code> and has the right size. */
    public static final double [] concat( double [] target,
					  final double [] con, final double [] cat ) {
	if( target.length != con.length + cat.length )
	    target = new double[ con.length + cat.length ];

	System.arraycopy( con, 0, target, 0,          con.length );
	System.arraycopy( cat, 0, target, con.length, cat.length );

	return target;
    }

    /** Fills <i>target</i> with <i>value</i>. */
    public static final void fill( final double [] target, final double value ) {
	fill( target, 0, target.length, value );
    }

    /** Fills <i>target</i> from <i>postion</i> to <i>postion+length-1</i>
	with <i>value</i>. */
    public static final void fill( final double [] target, final int position,
				   final int length, final double value ) {
	if( length < 1 ) return;

	target[position + 0]=value;
	if( length < 2 ) return;
	target[position + 1]=value;
	if( length < 3 ) return;
	target[position + 2]=value;
	if( length < 4 ) return;
	target[position + 3]=value;
	if( length < 5 ) return;
	target[position + 4]=value;

	int i=5;
	for( ; 2*i < length; i *= 2 )
	    System.arraycopy( target, position, target, position + i, i );

	System.arraycopy( target, position, target, position + i, length - i );
    }

    /** Returns the minimum of <i>target</i>. */
    public static final double getMin( final double [] target ) {
	return getMin( target, 0, target.length );
    }

    /** Returns the minimum of <i>target</i> from <i>postion</i> to <i>postion+length-1</i>. */
    public static final double getMin( final double [] target,
				       final int position, final int length ) {
	double min = target[position];
	for( int i=0; i<position+length; i++ )
	    if( min > target[i] )
		min = target[i];
	return min;
    }


    /** Returns the maximum of <i>target</i>. */
    public static final double getMax( final double [] target ) {
	return getMax( target, 0, target.length );
    }

    /** Returns the maximum of <i>target</i> from <i>postion</i> to <i>postion+length-1</i>. */
    public static final double getMax( final double [] target,
				       final int position, final int length ) {
	double max = target[position];
	for( int i=0; i<position+length; i++ )
	    if( max > target[i] )
		max = target[i];
	return max;
    }

    /** Just returns <i>target</i> if its length equals <i>size</i>.
	Otherwise, this includes <i>target</i> is null, it creates
	an array of length <i>size</i> and returns this instead. */
    public static final float [] size( final float [] target, final int size ) {
	return target.length == size ? target : new float[size];
    }

    /** Like {@link Arrays#size(float[],int)} but copies also the contents of
	the <i>target</i> if nessecary. */
    public static final float [] resize( final float [] target, final int size ) {
	if( target.length == size )
	    return target;

	float [] newTarget = new float [size];

	System.arraycopy( target, 0, newTarget, 0, Math.min( size, target.length ) );

	return newTarget;
    }

    /** Returns an array containing <i>numOfTimes</i> copies of the
	contens of <i>values</i>.
	The method uses <i>target</i> if it is not <code>null</code> and has the right size. */
    public static final float [] unrole( float [] target,
					  final float [] values, final int numOfTimes ) {
	int length = values.length * numOfTimes;

	if( target.length != values.length * numOfTimes )
	    target = new float[ length ];

	System.arraycopy( values, 0, target, 0, values.length );

	int i=values.length;
	for( ; 2*i < length; i *= 2 )
	    System.arraycopy( target, 0, target, i, i );

	System.arraycopy( target, 0, target, i, length - i );

	return target;
    }

    /** Retuns an array containing the contents of <i>con</i> and
	<i>cat</i>.
	The method uses <i>target</i> if it is not <code>null</code> and has the right size. */
    public static final float [] concat( float [] target,
					  final float [] con, final float [] cat ) {
	if( target.length != con.length + cat.length )
	    target = new float[ con.length + cat.length ];

	System.arraycopy( con, 0, target, 0,          con.length );
	System.arraycopy( cat, 0, target, con.length, cat.length );

	return target;
    }

    /** Fills <i>target</i> with <i>value</i>. */
    public static final void fill( final float [] target, final float value ) {
	fill( target, 0, target.length, value );
    }

    /** Fills <i>target</i> from <i>postion</i> to <i>postion+length-1</i>
	with <i>value</i>. */
    public static final void fill( final float [] target, final int position,
				   final int length, final float value ) {
	if( length < 1 ) return;

	target[position + 0]=value;
	if( length < 2 ) return;
	target[position + 1]=value;
	if( length < 3 ) return;
	target[position + 2]=value;
	if( length < 4 ) return;
	target[position + 3]=value;
	if( length < 5 ) return;
	target[position + 4]=value;

	int i=5;
	for( ; 2*i < length; i *= 2 )
	    System.arraycopy( target, position, target, position + i, i );

	System.arraycopy( target, position, target, position + i, length - i );
    }

    /** Returns the minimum of <i>target</i>. */
    public static final float getMin( final float [] target ) {
	return getMin( target, 0, target.length );
    }

    /** Returns the minimum of <i>target</i> from <i>postion</i> to <i>postion+length-1</i>. */
    public static final float getMin( final float [] target,
				       final int position, final int length ) {
	float min = target[position];
	for( int i=0; i<position+length; i++ )
	    if( min > target[i] )
		min = target[i];
	return min;
    }


    /** Returns the maximum of <i>target</i>. */
    public static final float getMax( final float [] target ) {
	return getMax( target, 0, target.length );
    }

    /** Returns the maximum of <i>target</i> from <i>postion</i> to <i>postion+length-1</i>. */
    public static final float getMax( final float [] target,
				       final int position, final int length ) {
	float max = target[position];
	for( int i=0; i<position+length; i++ )
	    if( max > target[i] )
		max = target[i];
	return max;
    }

    /** Just returns <i>target</i> if its length equals <i>size</i>.
	Otherwise, this includes <i>target</i> is null, it creates
	an array of length <i>size</i> and returns this instead. */
    public static final int [] size( final int [] target, final int size ) {
	return target.length == size ? target : new int[size];
    }

    /** Like {@link Arrays#size(int[],int)} but copies also the contents of
	the <i>target</i> if nessecary. */
    public static final int [] resize( final int [] target, final int size ) {
	if( target.length == size )
	    return target;

	int [] newTarget = new int [size];

	System.arraycopy( target, 0, newTarget, 0, Math.min( size, target.length ) );

	return newTarget;
    }

    /** Returns an array containing <i>numOfTimes</i> copies of the
	contens of <i>values</i>.
	The method uses <i>target</i> if it is not <code>null</code> and has the right size. */
    public static final int [] unrole( int [] target,
					  final int [] values, final int numOfTimes ) {
	int length = values.length * numOfTimes;

	if( target.length != values.length * numOfTimes )
	    target = new int[ length ];

	System.arraycopy( values, 0, target, 0, values.length );

	int i=values.length;
	for( ; 2*i < length; i *= 2 )
	    System.arraycopy( target, 0, target, i, i );

	System.arraycopy( target, 0, target, i, length - i );

	return target;
    }

    /** Retuns an array containing the contents of <i>con</i> and
	<i>cat</i>.
	The method uses <i>target</i> if it is not <code>null</code> and has the right size. */
    public static final int [] concat( int [] target,
					  final int [] con, final int [] cat ) {
	if( target.length != con.length + cat.length )
	    target = new int[ con.length + cat.length ];

	System.arraycopy( con, 0, target, 0,          con.length );
	System.arraycopy( cat, 0, target, con.length, cat.length );

	return target;
    }

    /** Fills <i>target</i> with <i>value</i>. */
    public static final void fill( final int [] target, final int value ) {
	fill( target, 0, target.length, value );
    }

    /** Fills <i>target</i> from <i>postion</i> to <i>postion+length-1</i>
	with <i>value</i>. */
    public static final void fill( final int [] target, final int position,
				   final int length, final int value ) {
	if( length < 1 ) return;

	target[position + 0]=value;
	if( length < 2 ) return;
	target[position + 1]=value;
	if( length < 3 ) return;
	target[position + 2]=value;
	if( length < 4 ) return;
	target[position + 3]=value;
	if( length < 5 ) return;
	target[position + 4]=value;

	int i=5;
	for( ; 2*i < length; i *= 2 )
	    System.arraycopy( target, position, target, position + i, i );

	System.arraycopy( target, position, target, position + i, length - i );
    }

    /** Returns the minimum of <i>target/<i>. */
    public static final int getMin( final int [] target ) {
	return getMin( target, 0, target.length );
    }

    /** Returns the minimum of <i>target</i> from <i>postion</i> to <i>postion+length-1</i>. */
    public static final int getMin( final int [] target,
				       final int position, final int length ) {
	int min = target[position];
	for( int i=0; i<position+length; i++ )
	    if( min > target[i] )
		min = target[i];
	return min;
    }


    /** Returns the maximum of <i>target</i>. */
    public static final int getMax( final int [] target ) {
	return getMax( target, 0, target.length );
    }

    /** Returns the maximum of <i>target</i> from <i>postion</i> to <i>postion+length-1</i>. */
    public static final int getMax( final int [] target,
				       final int position, final int length ) {
	int  max = target[position];
	for( int i=0; i<position+length; i++ )
	    if( max > target[i] )
		max = target[i];
	return max;
    }


}
