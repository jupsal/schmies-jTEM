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

package de.jtem.mfc.matrix;

import de.jtem.mfc.field.Complex;

/** This interface provides complex functions with values
 *  in complex 2-by-2 matrices.
 *  
 * @author schmies
 */

public interface Complex2By2Function {

    /** Field.Complex identity matrix function */
    public static Complex2By2Function ID = new Complex2By2Function() {
	    public Complex2By2 eval( Complex z ) { 
		return new Complex2By2();
	    }
	    public void eval( Complex z, Complex2By2 result ) {
		result.assignIdentity();
	    }
	};

    /** Field.Complex zero matrix function */
    public static Complex2By2Function ZERO = new Complex2By2Function() {
	    public Complex2By2 eval( Complex z ) {
		return new Complex2By2( 0, 0, 0, 0, 0, 0, 0, 0 );
	    }
	    public void eval( Complex z, Complex2By2 result ) {
		result.assignZero();
	    }
	};

 
    /* Evalaluates this at <code>z</code> and returns result in a new Complex2By2. */
    public Complex2By2 eval( Complex z );

    /* Evalaluates this at <code>z</code> and returns result in<code>result</code> . */
    public void eval( Complex z, Complex2By2 result );

}
