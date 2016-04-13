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
import de.jtem.blas.IntegerMatrix;
import de.jtem.mfc.field.Complex;
import junit.framework.TestCase;

/**
 * Testing directly by evaluating thata functions.
 * @author Schmies
 * @version 1.0
 */
public class TestModularPropertySupport extends TestCase {

  int n = 2;

  SiegelReduction sr = new SiegelReduction( n );

  ComplexMatrix B = new ComplexMatrix( 2 );

  ModularPropertySupport mps;

  ModularTransformation modular;
  Theta theta, transformedTheta;

  ComplexVector Z;


  public void setUp() {

    n = 2;

    sr = new SiegelReduction(n);

    B = new ComplexMatrix(2);

    B.set(0, 0, -5.10972365633887, Math.PI);
    B.set(0, 1, -4.24199777781055, 0);

    B.set(1, 0, -4.24199777781055, 0);
    B.set(1, 1, -5.78346380443502, Math.PI);

    sr.setPeriodMatrix(B);

    mps = new ModularPropertySupport(n);

    modular = sr.getModularTransformation();

    mps.setModularTransformation(modular);

    mps.setPeriodMatrix(B);

    theta = new Theta(B, 1e-16, false, false);
    transformedTheta = new Theta(mps.getTransformedPeriodMatrix(), 1e-16, false, false);

    Z = new ComplexVector(n);

  }

  public void testTheta() {

    Z.assignZero();
    mps.setZ(Z);
    Complex k = transformedTheta.theta(mps.getTransformedZ()).
        times(mps.getFactor().exp()).divide(theta.theta(Z));

    for (int i = 0; i < 10; i++) {
      Z.assignRandom();
      Z.assignTimes(10);
      mps.setZ(Z);

      Complex f = transformedTheta.theta(mps.getTransformedZ()).
          times(mps.getFactor().exp()).divide(theta.theta(Z));

      assertEquals("real part", k.re, f.re, 1e-13);
      assertEquals("imag part", k.im, f.im, 1e-13);

    }
  }
  
  public void  testPaperExample() {
  	
  	//static public void main( String []arg ) {
  	int n = 2;

  	ModularPropertySupport mps = new ModularPropertySupport(n);
  	
  	ComplexMatrix B = new ComplexMatrix(2);

  	// example with reducable lattice
  	B.set(0, 0, -111.207, 0);
  	B.set(0, 1, -96.616, 0 );
  	B.set(1, 0, -96.616, 0 );
  	B.set(1, 1, -111.207, 0);

  	mps.setPeriodMatrix( B );

  	// compare with data of modular transformation
  	IntegerMatrix a = new IntegerMatrix( new int[][] { { 0, 0 }, { 0, 0 } } );
  	IntegerMatrix b = new IntegerMatrix( new int[][] { { 8, 7 }, { 7, 6 } } );
  	IntegerMatrix c = new IntegerMatrix( new int[][] { {6, -7 }, { -7,8 } } );
  	IntegerMatrix d = new IntegerMatrix( new int[][] { { 0, 0 }, { 0, 0 } } );

  	ModularTransformation modular = new ModularTransformation( a,b,c,d);
  	
  	mps.setModularTransformation(modular);
  	
  	System.out.println( mps.getTransformedPeriodMatrix().divide( 2*Math.PI) );
  	
  	B.assignDivide( new de.jtem.mfc.field.Complex( 0, 2*Math.PI) );
  	
  	System.out.println( a.times(B).plus(b).times( c.times(B).plus(d).invert()  ));
  	
  }
  

}















