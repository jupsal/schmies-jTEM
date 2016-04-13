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

package de.jtem.riemann.schottky;

import de.jtem.mfc.field.Complex;
import junit.framework.TestCase;

public class SchottkyGroupElementTest extends TestCase {

  SchottkyGroupElement element = new SchottkyGroupElement();

  Complex z, w, d;

  public void setUp() {

    element.assign( new Complex(  3.4, -5.6 ),
                    new Complex( -4.2,  1.3 ), new Complex( 0.034, -0.0013 ) );

    z = new Complex(  -1.1, 4.3 );
    w = new Complex(  3.2, -4.2 );
    d = new Complex();
  }


  public void testDiff() {
    element.diff( z, w, d );

    Complex D = element.applyTo(z).minus( element.applyTo(w) );

    assertEquals( 0, d.dist( D ), 1e-14 );
  }

  public void testDiffPow() {

    for( int i=1; i<16; i++ ) {

      Complex D = element.applyTo(z).pow( i ).minus( element.applyTo(w).pow( i ) );

      element.diff( z, w, d );
      element.diffPow( z, w, i, d, d );

      assertEquals( "pow of " + i, 0, d.dist( D ) / D.abs(), 1e-14 );

      element.diffPow( z, w, i, d );

      assertEquals( "pow of " + i, 0, d.dist( D ) / D.abs(), 1e-14 );

    }
  }

  public void testInverseOfCSqr() {

    Complex inverseOfCSqr = new Complex();

    element.getInverseOfCSqr( inverseOfCSqr );

    assertEquals( 0, inverseOfCSqr.dist( element.getC().sqr().invert() ), 1e-14 );
  }
}
