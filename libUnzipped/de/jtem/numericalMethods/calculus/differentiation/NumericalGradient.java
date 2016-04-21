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

package de.jtem.numericalMethods.calculus.differentiation;

import de.jtem.numericalMethods.calculus.function.RealFunctionOfSeveralVariables;
import de.jtem.numericalMethods.calculus.function.RealFunctionOfSeveralVariablesWithGradient;


/**
 * Real function of several variables with numerical gradient.
 * 
 * @author schmies
 */
final public class NumericalGradient implements RealFunctionOfSeveralVariablesWithGradient, java.io.Serializable  {

  private static final long serialVersionUID = 1L;

  RealFunctionOfSeveralVariables f;
  int n;
  double h;

  public NumericalGradient(final RealFunctionOfSeveralVariables theF) {
    this ( theF, 1e-8);
  }

  public NumericalGradient(final RealFunctionOfSeveralVariables theF,
                           final double theH) {
    f = theF;
    n = theF.getNumberOfVariables();
    h = theH;
  }

  /**
   * Get the value of h.
   * @return Value of h.
   */
  public double getH() {
    return h;
  }

  /**
   * Set the value of h.
   * @param v  Value to assign to h.
   */
  public void setH(double v) {
    this.h = v;
  }

  /**
   * Returns number of variables which must be a positive number.
   * @return number of variables.
   */
  public int getNumberOfVariables() {
    return n;
  }

  /**
   * Evaluates function at x.
   * @param x
   * @return value of function at x.
   */
  public double eval(final double[] x) {
    return f.eval(x);
  }

  /**
   * Evaluates function at x and saves the gradient at point x.
   * @param x
   * @param gradient gradient on output
   * @return value of function at x.
   */
  public double eval(double[] x, double[] grad) {
    return computeGradient( f, x, grad, h );
  }
  
  /**
   * Computes gradient of f at x using asymmetric differences with size h. 
   * Evaluates: ( f(x + h e_i) - f(x) ) / h ) e_i .
   * @param x
   * @param gradient gradient on output
   * @param h 
   * @return value of function at x.
   */
  public static double computeGradient( RealFunctionOfSeveralVariables f, 
  		double[] x, double[] grad, double h ) {

  	final int n=f.getNumberOfVariables();
  	
  	double fp = f.eval( x );

  	for (int i = 0; i < n; i++) {
  		final double tmp = x[i];
  		x[i] += h;
  		final double h_ =x[i] - tmp;
  		grad[i] = ( f.eval( x ) -fp ) /h_;
  		x[i] = tmp;
  	}

  	return fp;
  }

}




