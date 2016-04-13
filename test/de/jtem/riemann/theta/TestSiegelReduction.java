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
import de.jtem.blas.IntegerMatrix;
import junit.framework.TestCase;

public class TestSiegelReduction extends TestCase {

    public void test2DExample() {
    //static public void main( String []arg ) {
      int n = 2;

      SiegelReduction sr = new SiegelReduction(n);

      ComplexMatrix B = new ComplexMatrix(2);

      // example with reducable lattice
      B.set(0, 0, -5.109723656338871, Math.PI);
      B.set(0, 1, -4.2419977778105515, 0);
      B.set(1, 0, -4.2419977778105515, 0);
      B.set(1, 1, -5.78346380443502, Math.PI);

      sr.setPeriodMatrix( B );

      // compare with data of modular transformation
      IntegerMatrix a = new IntegerMatrix( new int[][] { { 1, 0 }, { 0, 1 } } );
      IntegerMatrix b = new IntegerMatrix( new int[][] { { 0, 0 }, { 0, 0 } } );
      IntegerMatrix c = new IntegerMatrix( new int[][] { {-1, 1 }, { 1,-1 } } );
      IntegerMatrix d = new IntegerMatrix( new int[][] { { 1, 0 }, { 0, 1 } } );

      assertTrue( "a", a.equals( sr.modular.a ) );
      assertTrue( "b", b.equals( sr.modular.b ) );
      assertTrue( "c", c.equals( sr.modular.c ) );
      assertTrue( "d", d.equals( sr.modular.d ) );

      // compare with data of reduced period matrix
      ComplexMatrix reduced = sr.getReducedPeriodMatrix();

      assertEquals( "reduced period matrix: real of entry 00", -8.893837400587667, reduced.re[0][0], 1e-14 );
      assertEquals( "reduced period matrix: real of entry 11", -8.893837400587667, reduced.re[1][1], 1e-14 );
      assertEquals( "reduced period matrix: real of entry 01", -0.700547044613776, reduced.re[0][1], 1e-14 );
      assertEquals( "reduced period matrix: real of entry 10", -0.700547044613776, reduced.re[1][0], 1e-14 );

      assertEquals( "reduced period matrix: imag of entry 00",  0.878558945495512, reduced.im[0][0], 1e-14 );
      assertEquals( "reduced period matrix: imag of entry 11", -0.878558945495512, reduced.im[1][1], 1e-14 );
      assertEquals( "reduced period matrix: imag of entry 01",  3.141592653589791, reduced.im[0][1], 1e-14 );
      assertEquals( "reduced period matrix: imag of entry 10",  3.141592653589791, reduced.im[1][0], 1e-14 );
    }
    
    

    public void test2DExample2() {
    	//static public void main( String []arg ) {
    	int n = 2;

    	SiegelReduction sr = new SiegelReduction(n);

    	ComplexMatrix B = new ComplexMatrix(2);

    	// example with reducable lattice
    	B.set(0, 0, -111.207, 0);
    	B.set(0, 1, -96.616, 0 );
    	B.set(1, 0, -96.616, 0 );
    	B.set(1, 1, -111.207, 0);

    	sr.setPeriodMatrix( B );

    	// compare with data of modular transformation
    	IntegerMatrix a = new IntegerMatrix( new int[][] { { -1, 1 }, { 0, 1 } } );
    	IntegerMatrix b = new IntegerMatrix( new int[][] { { 0, 0 }, { 0, 0 } } );
    	IntegerMatrix c = new IntegerMatrix( new int[][] { {0, 0 }, { 0, 0 } } );
    	IntegerMatrix d = new IntegerMatrix( new int[][] { { -1, 0 }, { 1, 1 } } );

    	System.out.println( sr.modular );
    	
    	System.out.println( sr.getReducedPeriodMatrix().divide( 2*Math.PI ));
    	
    	assertTrue( "a", a.equals( sr.modular.a ) );
    	assertTrue( "b", b.equals( sr.modular.b ) );
    	assertTrue( "c", c.equals( sr.modular.c ) );
    	assertTrue( "d", d.equals( sr.modular.d ) );

    	// compare with data of reduced period matrix
    	ComplexMatrix reduced = sr.getReducedPeriodMatrix();

//    	assertEquals( "reduced period matrix: real of entry 00", -8.893837400587667, reduced.re[0][0], 1e-14 );
//    	assertEquals( "reduced period matrix: real of entry 11", -8.893837400587667, reduced.re[1][1], 1e-14 );
//    	assertEquals( "reduced period matrix: real of entry 01", -0.700547044613776, reduced.re[0][1], 1e-14 );
//    	assertEquals( "reduced period matrix: real of entry 10", -0.700547044613776, reduced.re[1][0], 1e-14 );
//
//    	assertEquals( "reduced period matrix: imag of entry 00",  0.878558945495512, reduced.im[0][0], 1e-14 );
//    	assertEquals( "reduced period matrix: imag of entry 11", -0.878558945495512, reduced.im[1][1], 1e-14 );
//    	assertEquals( "reduced period matrix: imag of entry 01",  3.141592653589791, reduced.im[0][1], 1e-14 );
//    	assertEquals( "reduced period matrix: imag of entry 10",  3.141592653589791, reduced.im[1][0], 1e-14 );
    }
}
















