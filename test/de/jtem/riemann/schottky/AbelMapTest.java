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

import de.jtem.blas.ComplexVector;
import de.jtem.mfc.field.Complex;
import junit.framework.TestCase;

public class AbelMapTest extends TestCase {

  /**
   * Computes spiral with constant radius
   */
  static ComplexVector spiral( Complex center,
                               double radius,
                               double startAngle,
                               double endAngle ) {

    final double delta = endAngle - startAngle;

    if( Math.abs( delta ) < 1e-12 )
      return new ComplexVector();

    final int n = (int)(Math.abs( delta ) / 2 / Math.PI * 10 + 1);

    ComplexVector spiral = new ComplexVector( n + 1 );

    for( int i=0; i<=n; i++ ) {
      final double angle = startAngle + i * delta / n;

      spiral.set( i,
                  center.getRe() + radius * Math.cos( angle ) ,
                  center.getIm() + radius * Math.sin( angle ) );
    }

    return spiral;
  }

  /**
   * Test analytic integral of differential with index k by integrating
   * along a growing spiral around kth a-circle and checks periods.
   * @param g genus of helicoid to test.
   */
  public void testHelicoid( int g ) {

    Schottky schottky = TestSchottky.getSchottkyOfHe( g );

    AbelMap abel = new AbelMap( schottky );

    final int n = schottky.getNumGenerators(); // should be g

    for( int i=0; i<n; i++ ) {

      for( int k=0; k < 4; k++ ) {

        ComplexVector path = spiral(schottky.center[i][1],
                                    schottky.radius[i] * 1.2, 0,
                                    k * Math.PI / 2);

        Complex refIntegral = abel.eval(i, path, 1e-10 );

        for( int l=1; l<4; l++ ) {

          refIntegral.im += 2 * Math.PI;

          ComplexVector path_ = spiral(schottky.center[i][1],
                        schottky.radius[i] * 1.2, 0,
                        (k + 4 * l) * Math.PI / 2);

          if( k>0 )
            path_.set( path_.size()-1, path.get(path.size()-1)); //assure ends coinside

          Complex integral = abel.eval(i, path_, 1e-10);

          assertEquals("real part of integral", refIntegral.re, integral.re, 1e-14 );
          assertEquals("imag part of integral", refIntegral.im, integral.im, 1e-14 );
        }
      }
    }
  }

  public void testHe1() {
    testHelicoid( 1 );
  }

  public void testHe2() {
    testHelicoid(2);
  }

  public void testHe3() {
    testHelicoid(3);
  }

  public void testHe4() {
    testHelicoid(4);
  }

}
