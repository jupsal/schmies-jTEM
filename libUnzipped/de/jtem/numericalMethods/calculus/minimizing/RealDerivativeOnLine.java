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
import de.jtem.numericalMethods.calculus.function.RealFunctionOfSeveralVariablesWithGradient;

/**
 * This class represents an routine of calculating an
 * @author Markus Schmies, Vitali Lieder
 * @version 1.0
 */
class RealDerivativeOnLine implements RealFunctionOfOneVariable, java.io.Serializable {

  private static final long serialVersionUID = 1L;

  final RealFunctionOfSeveralVariablesWithGradient f;

  final Line line;

  final double [] gradAtX;

  RealDerivativeOnLine( Line line, RealFunctionOfSeveralVariablesWithGradient f) {
    this.line = line;
    this.f = f;
    gradAtX = new double[ line.n ];
  }

  RealDerivativeOnLine(double[] point, double[] direction,
                     RealFunctionOfSeveralVariablesWithGradient f) {
    this( new Line( point, direction), f );
  }

  public final double eval(double x) {
    line.getPoint(x, line.otherPoint);
    f.eval( line.otherPoint, gradAtX );
    double derivative = 0;
    for( int i=0; i<line.n; i++ ) {
      derivative += gradAtX[i] * line.direction[i];
    }
    return derivative;
  }
}
