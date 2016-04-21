/**
This file is part of a jTEM project.
All jTEM projects are licensed under the FreeBSD license 
or 2-clause BSD license (see http://www.opensource.org/licenses/bsd-license.php). 

Copyright (c) 2002-2010, Technische Universität Berlin, jTEM
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

package de.jtem.blas;

import java.io.Serializable;

abstract class AbstractVector implements Serializable,Cloneable
{
    private static final long serialVersionUID = 1L;

    // are entries equal ?    a "=" b  <=>  |a-b|^2 < EPSILON
    protected static double EPSILON = 1.e-28;

    int length;

    /** return total number of entries  */
    public final int size()
    {
	return length;
    }

    /** return total number of entries  */
    public final int getNumEntries() {
	return length;
    }
    
    /** set number of entries  */
    public final void setNumEntries( int aNum) {
	resize( aNum );
    }

    /** throw <b>IllegalArgumentException</b> if <i>anIndex</i> exceeds number of entries */
    public final void checkIndex( int anIndex ) {
	if( anIndex < 0 || anIndex >= length  )
	    throw new IllegalArgumentException("index "+anIndex+" is not in [0,"+(length-1)+"]");
    }

    /** throw <b>IllegalArgumentException</b> if <i>aLength</i> does not equal <i>this</i>' length  */
    public final void checkLength( int aLength ) {
	if( aLength != length )
	    throw new IllegalArgumentException(aLength+" is not in equal to "+length);
    }
    
   /** throw <b>IllegalArgumentException</b> if <i>V</i>' length does not equal <i>this</i>' length  */
    public final void checkShape( final AbstractVector V )
	{
	    if(this.length!=V.length)
		throw new IllegalArgumentException("different vec sizes");
	}  

    /** return <tt>true</tt> if <i>V</i> and <i>this</i> have same shape, <tt>false</tt> o/w */
    public final boolean hasSameShape( final AbstractVector V ) {
	return V.length == length;
    }

    /** return <tt>true</tt> if <i>V</i> and <i>this</i> have same shape, <tt>false</tt> o/w */
    public final boolean hasSameShape( final double[] V ) {
	return length == V.length;
    }

    /** DEPRECATED ! */
    public boolean equals(Object o)
    {
	try
	    {
		AbstractVector v=(AbstractVector)o;
		if(v.length!=length) return false;
		return true;
	    }
	catch(ClassCastException ex)
	    {
		return false;
	    }
    }


    abstract void newSize( int newNumEntries, boolean copyValues );

    /** resize vector to <i>newNumEntries</i>, but ignore the values.*/
    public void newSize( int newNumEntries ) {
	newSize( newNumEntries, false );	
    }

    /** resize vector to <i>aVector</i>'s length, but ignore the values. */
    public final void newSize( AbstractVector aVector ) {
	newSize( aVector.length );
    }
    
    /** resize vector to <i>newNumEntries</i>
         w/o changing its (surviving) entries, possibly adding zero entries */
    public void resize( int newNumEntries ) {
	newSize( newNumEntries, true );
    }

    /** resize vector to <i>aVector</i>'s length
        w/o changing its (surviving) entries, possibly adding zero entries               */
    public final void resize( AbstractVector aVector ) {
	newSize( aVector.length, true );
    }

    /** resize vector to <i>newNumEntries</i>
         w/o changing its (surviving) entries, possibly adding zero entries */
    public void setSize( int newNumEntries ) {
        newSize(newNumEntries, true);
    }

    /** resize vector to <i>aVector</i>'s length
        w/o changing its (surviving) entries, possibly adding zero entries               */
    public final void setSize( AbstractVector aVector ) {
	newSize( aVector.length, true );
    }

    /** return the total sum of squares of entries' length */
    public abstract double normSqr();

    /** return <tt>true</tt> if <i>this</i> has only zero entries */
    public boolean isZero() {
      return (normSqr() <= length * EPSILON);
    }

    /** return a clone */
    public abstract Object clone();


    /** print this vector */
    public void print( String title ) {
	
	System.out.println( title );
	
	System.out.println( toString() );
    }
}

   









