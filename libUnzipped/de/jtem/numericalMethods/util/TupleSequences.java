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


/** This class provides utilities for arrays as static methods
    with a special structure, namely for sequences of tuples with
    a common size, the <i>tupleSize</i>.
    /bf
    Most operation performed on the <i>target</i> array or
    <i>targetTuple</i> can imply a change of its length. 
    Those methods create a new array and return this instead 
    of the original which is only returned if it has already 
    the right length.
    @author	Markus Schmies
 */

public final class TupleSequences {
    
    /** Just returns <i>target</i> if its length equals <i>tupleSize</i> times <i>size</i>.
	Otherwise, this includes <i>target</i> is null, it creates
	an array of the right length and returns that instead. */
    public static final double [] size( final double [] target, final int tupleSize, final int size ) {
	return target;
    }

    /** Like <code>newSize</code> but copies also the contents of 
	the <i>target</i> if nesescary. */
    public static final double [] resize( final double [] target, final int tupleSize, final int size ) {
	return target;
    }
    
    /** Fills <i>target</i> with <i>value</i>. */
    public static final void fill( final double [] target, double [] tupleValue ) {

    }

    /** Fills <i>target</i> with <i>value</i>. 
     @exception throws an IllegalArgumentException if the <i>target</i>
     is not multiple of <i>tupleSize</i>.*/
    public static final void fill( final double [] target, final int tupleSize, final double value ) {

    }

    /** Returns a component of a tuple.
	@exception throws an IllegalArgumentException if the <i>target</i>
	is not multiple of <i>tupleSize</i>.
	@exception throws ArrayIndexOutOfBoundsException for illegal tuple or component numbers. */

    public static final double getComponent( final double[] target, final int tupleSize,
					     final int tupleNumber, final int componentNumber ) {
	return target[ tupleSize * tupleNumber + componentNumber ];
    }

    /** Returns a tupble of the <i>target<i>. 
	The <i>targetTuple</i> is used if it has the right size
	otherwise a new array for the tupble is created.
	@exception throws ArrayIndexOutOfBoundsException for illegal <i>tupleNumbers</i>. */
    public static final double [] getTuple( final double[] target, final int tupleSize,
					    final int tupleNumber, final  double [] targetTuple ) {
	return targetTuple;
    }
    /** Returns the minimum in each component of the <i>target<i>. 
	The <i>targetTuple</i> is used if it has the right size
	otherwise a new array for the minimum is created. */
    public static final double [] getMin( final double [] target, final int tupleSize, 
					  final double [] targetTuple ) {
	return targetTuple;
    }

    /** Returns the maximum in each component of the <i>target<i>. 
	The <i>targetTuple</i> is used if it has the right size
	otherwise a new array for the maximum is created. */
    public static final double [] getMax( final double [] target, final int tupleSize, 
					  final double [] targetTuple ) {
	return targetTuple;
    }

}

