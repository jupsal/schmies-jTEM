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

package de.jtem.numericalMethods.calculus.minimizing;

import de.jtem.numericalMethods.calculus.function.RealFunctionOfOneVariable;
import de.jtem.numericalMethods.calculus.function.RealFunctionOfSeveralVariables;

/**
 * This class represents the routine of calculating a minimum of the multidimensional function across
 * one line.
 * @author Markus Schmies, Vitali Lieder
 * @version 1.0
 */
class RealFunctionOnLine implements RealFunctionOfOneVariable, java.io.Serializable {

  private static final long serialVersionUID = 1L;

  /**
   * The multi dimensional function f
   */
  final RealFunctionOfSeveralVariables f;

  /**
   * line across we move
   */
  final Line line;

  RealFunctionOnLine( Line line, RealFunctionOfSeveralVariables f) {

    this.line = line;

    this.f = f;
  }

  RealFunctionOnLine(double[] point, double[] direction,
                     RealFunctionOfSeveralVariables f) {

    this( new Line( point, direction), f );
  }

  /**
   * Return the value of function f at the line in x direction.
   * @param x direction at the line.
   * @return the value of
   */
  public final double eval(double x) {
    line.getPoint(x, line.otherPoint);
    return f.eval(line.otherPoint);
  }

  /**
   * Moves the point of the line where function f has a minimum along the line
   * and return the value of it. Actually doing this by making braking on interval (-1,0,1) and then searching
   * minimum in the bracketing interval with brent's methods.
   * @param tol tolerance of calculation
   * @return the value of the minimum on line
   */
  final double brent(double tol) {

    final double[] abc = new double[3];
    final double[] valuesAtABC = new double[3];
    final double[] result = new double[2];

    abc[0] = -1.0;
    abc[1] = 0.0;
    abc[2] = 1.0;

    Braket.search(abc, valuesAtABC, this);

    Brent.search(abc, result, this, tol);

    final double xmin = result[0];

    for (int j = 0; j < line.n; j++) {
      line.point[j] += (line.direction[j] *= xmin);
    }

    return result[1];
  }


}
