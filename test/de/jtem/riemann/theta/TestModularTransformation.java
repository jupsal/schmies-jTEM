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

import de.jtem.blas.IntegerMatrix;
import junit.framework.TestCase;

/**
 * The testing is performed by producing the different generators
 * randomly and applieng them to the different generating methods.
 */

public class TestModularTransformation extends TestCase  {

    public void test() {

      for( int j=0; j<10; j++ ) {

        ModularTransformation m = new ModularTransformation(3);

        IntegerMatrix a = new IntegerMatrix(3);

        assertTrue( m.respectsSymplecticStructure() );

        for (int i = 0; i < 10; i++) {

          m.applyGenerator(getRandomInvertableIntegerMatrix(3), null);

          assertTrue( m.respectsSymplecticStructure() );

          m.applyGenerator(getRandomSymmetricIntegerMatrix(3));

          assertTrue( m.respectsSymplecticStructure() );

          m.applyGenerator();

          assertTrue( m.respectsSymplecticStructure() );

          m.applyGenerator(getRandomSymmetricIntegerMatrix(3));

          assertTrue( m.respectsSymplecticStructure() );

          m.applyGenerator(getRandomInvertableIntegerMatrix(3), null);

          assertTrue( m.respectsSymplecticStructure() );

          m.applyGenerator(getRandomSymmetricIntegerMatrix(3));

          assertTrue( m.respectsSymplecticStructure() );

          m.applyGenerator();
        }
      }
    }


    static IntegerMatrix getRandomInvertableIntegerMatrix(int n) {

      IntegerMatrix a = new IntegerMatrix(n);

      a.assignRandom(4);
      while (a.determinant() != 1 || a.minus(IntegerMatrix.id(n)).isZero()) {
        a.assignRandom(4);
      }

      return a;
    }

    static IntegerMatrix getRandomSymmetricIntegerMatrix(int n) {

      IntegerMatrix a = new IntegerMatrix(n);

      a.assignRandom(2); ;
      a.assignPlus(a.transpose());

      while (a.isZero()) {
        a.assignRandom(2);
        a.assignPlus(a.transpose());
      }

      return a;
    }

}
















