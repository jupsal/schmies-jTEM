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

import de.jtem.numericalMethods.calculus.function.RealFunctionOfSeveralVariables;

/**
 * This abstract class represents a routine to find a minimum of given mulitdimensional function f
 * on line.
 * @author Markus Schmies, Vitali Lieder
 * @version 1.0
 */
public abstract class MinimizingOnLine implements java.io.Serializable {

  private static final long serialVersionUID = 1L;

  final Line line;

  final RealFunctionOnLine g;

  final double[] abc = new double[3];
  final double[] valuesAtABC = new double[3];
  final double[] result = new double[2];

  /**
   * Constructor with line along which minimum will be find and given function.
   * @param line given line along which minimum will be searched.
   * @param f given multidimensional function.
   */
  public MinimizingOnLine( Line line, RealFunctionOfSeveralVariables f) {

    this.line = line;

    this.g = new RealFunctionOnLine( line, f );
  }

  /**
   * This constructor build a new line with given point and direction.
   * @param point the initial point od line.
   * @param direction direction of line.
   * @param f given multidimensional function.
   */
  public MinimizingOnLine( double[] point, double[] direction,
                     RealFunctionOfSeveralVariables f) {

    this( new Line( point, direction), f );
  }

  /**
   * Routine to find minimum along given line fo given function f.
   * @param tol the precision.
   * @return the value of minimum on line
   */
  public abstract double search( double tol );
}
