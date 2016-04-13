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

package de.jtem.riemann.theta;

import de.jtem.blas.ComplexMatrix;
import de.jtem.blas.ComplexVector;
import de.jtem.mfc.field.Complex;
import junit.framework.TestCase;

/**
 * We test the theta functions via its modular transformation property.
 * Theta allows to switch off the optimazation which is brought by the
 * siegel reduction. The choosen period matrix generates a reduceable lattices,
 * we evaluate theta, its first and second derivative with and without
 * the use of the siegelreduction. The result should all differ by the
 * same factor which should be a 8th root of unity.
 * @author Markus Schmies
 * @version 0.9
 */
public class TestTheta extends TestCase {

    int n;

    ComplexMatrix B;

    Theta theta;
    Theta siegelTheta;

    ComplexVector Z, X, Y;

    public void setUp() {

          n = 2;

          B = new ComplexMatrix( 2 );

          // lettice of the following period matrix is reducable
          B.set( 0, 0, -5.10972365633887, Math.PI );
          B.set( 0, 1, -4.24199777781055, 0       );

          B.set( 1, 0, -4.24199777781055, 0       );
          B.set( 1, 1, -5.78346380443502, Math.PI );

          theta       = new Theta( B, 1e-16, false );  // without lattice reduction
          siegelTheta = new Theta( B, 1e-16, true );   // with    lattice reduction

          theta.setUniformApproximation(false);
          
          Z = new ComplexVector( n );
          X = new ComplexVector( n );
          Y = new ComplexVector( n );

          X.set( 1, 1, 0 );
          Y.set( 0, 1, 0 );

    }

    /** test derivative of theta function using the modular property*/
    public void testTheta() {

      Complex k = siegelTheta.theta(Z).divide( theta.theta(Z) );

      for( int i=0; i<15; i++ ) {

        Z.assignRandom();
        Z.assignTimes( 10 );
        
        Complex f =  siegelTheta.theta(Z).divide( theta.theta(Z) );

        assertEquals( "real part", k.re, f.re, 1e-13 );
        assertEquals( "imag part", k.im, f.im, 1e-13 );
      }
    }

    /** test derivative of theta function using the modular property*/
    public void testDTheta() {

      Complex k = siegelTheta.theta(Z).divide( theta.theta(Z) );

      for( int i=0; i<15; i++ ) {

        Z.assignRandom();
        Z.assignTimes( 10 );

        Complex f =  siegelTheta.dTheta(Z,X).divide( theta.dTheta(Z,X) );

        assertEquals( "real part", k.re, f.re, 1e-13 );
        assertEquals( "imag part", k.im, f.im, 1e-13 );
      }
    }

    /** test second derivative of theta function using the modular property*/
    public void testDDTheta() {

      Complex k = siegelTheta.theta(Z).divide( theta.theta(Z) );

      for( int i=0; i<15; i++ ) {

        Z.assignRandom();
        Z.assignTimes( 10 );

        Complex f =  siegelTheta.ddTheta(Z,X,Y).divide( theta.ddTheta(Z,X,Y) );

        assertEquals( "real part", k.re, f.re, 1e-13 );
        assertEquals( "imag part", k.im, f.im, 1e-13 );
      }
    }
}
















