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

package de.jtem.numericalMethods.geometry.geodesic;

/**
 * @author schmies
 *
 * To change the template for this generated type comment go to
 * Window - Preferences - Java - Code Generation - Code and Comments
 */
public class DoubleTriple {

	DoubleTriple() {
		super();
	}
	

	static void add( double [] c, double [] a, double [] b ) { //TODO: rename or remove
		c[0]=a[0]+b[0];
		c[1]=a[1]+b[1];
		c[2]=a[2]+b[2];
	}
	
	static void sub( double [] c, double [] a, double [] b ) { //TODO: rename or remove
		c[0]=a[0]-b[0];
		c[1]=a[1]-b[1];
		c[2]=a[2]-b[2];
	}
	
	static void scale( double [] c, double [] a, double f ) { //TODO: rename or remove
		c[0]=a[0]*f;
		c[1]=a[1]*f;
		c[2]=a[2]*f;
	}
	
	static void copy( double [] dst, double [] src ) { //TODO: rename or remove
		dst[0]=src[0];
		dst[1]=src[1];
		dst[2]=src[2];
	}
	
	static void zero( double [] bary ) { //TODO: rename or remove
		bary[0] = bary[1] =bary[2] =0;
	}
	
	

}
