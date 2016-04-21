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

package de.jtem.numericalMethods.function;

/**
 * @deprecated
 * Interface for objects that represent <code>float</code> array valued
 * functions.
 */
public interface FloatArrayValued {

    /**
     * Get the length of the <code>float</code> array value.
     *
     * @return the length of the <code>float</code> array value.
     */
    public int getFloatArrayValueLength();
	
    /**
     * Writes the <code>float</code> array value into a
     * <code>float[]</code>. When the method returns,
     * <code>value[offset,&hellip;,offset+len-1]</code> contains the
     * <code>float</code> array value, where <code>len</code> is the value
     * returned by {@link #getFloatArrayValueLength()}. The other elements of
     * <code>value</code> should not be changed.
     *
     * @param value a <code>float[]</code> whith length at least
     * <code>offset+len</code>, into which the <code>float</code> array value
     * is copied.
     * @param offset the position in <code>value</code> where the
     * <code>float</code> array value is copied to.
     */
    public void getFloatArrayValue( float [] value, int offset );

}
