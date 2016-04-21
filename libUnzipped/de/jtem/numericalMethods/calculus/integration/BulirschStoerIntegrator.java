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

package de.jtem.numericalMethods.calculus.integration;


import de.jtem.numericalMethods.calculus.function.RealFunctionOfOneVariable;
import de.jtem.numericalMethods.calculus.function.RealVectorValuedFunctionOfOneVariable;
import de.jtem.numericalMethods.calculus.odeSolving.BulirschStoer;

/**
 * Integrator based on the adaptive ode-solver
 * {@link de.jtem.numericalMethods.calculus.odeSolving.BulirschStoer} .
 * For usage consult example class: {@link ExampleIntegration}.
 *
 * @see de.jtem.numericalMethods.calculus.odeSolving.BulirschStoer
 * @see ExampleIntegration
 * @author Markus Schmies
 * @version 1.0
 */
public final class BulirschStoerIntegrator
    extends OdeSolverBasedIntegrator {

  private static final long serialVersionUID = 1L;

  double initialStepSize = 0.1;

  double eps = 1e-7;

  public BulirschStoerIntegrator() {
    super(new BulirschStoer());
  }

  public BulirschStoerIntegrator(RealFunctionOfOneVariable f) {
    super(new BulirschStoer(1));
    setFunction(f);
  }

  public BulirschStoerIntegrator(RealVectorValuedFunctionOfOneVariable f) {
    super(new BulirschStoer(f.getDimensionOfTargetSpace()));
    setFunction(f);
  }

  public double getInitialStepSize() {
    return initialStepSize;
  }

  public void setInitialStepSize(double initialStepSize) {
    if (this.initialStepSize == initialStepSize) {
      return;
    }

    this.initialStepSize = initialStepSize;
  }

  public double getEps() {
    return eps;
  }

  public void setEps(double eps) {

    ( (BulirschStoer) odeSolver).setEps(eps);

    if (this.eps == eps) {
      return;
    }

    this.eps = eps;
  }

  /**
   * Integrates real function on interval.
   * @param f function
   * @param a start of interval
   * @param b end of interval
   * @return integral of f over [a,b] with numerical accuracy of 1e-12
   */
  public static double integrate(RealFunctionOfOneVariable f,
                                 double a, double b) {

    return integrate(f, a, b, 1e-12);
  }

  /**
   * Integrates real function on interval with prescribed accuracy.
   * @param f function
   * @param a start of interval
   * @param b end of interval
   * @param acc tolerance for result
   * @return integral of f over [a,b] with provided tol
   */
  public static double integrate(RealFunctionOfOneVariable f,
                                 double a, double b,
                                 double acc) {

    BulirschStoerIntegrator integrator = new BulirschStoerIntegrator(f);

    integrator.setEps(acc);

    return integrator.integrate(a, b);
  }

  /**
   * Integrates set of real function on interval.
   * @param f set of real function
   * @param a start of interval
   * @param b end of interval
   * @param integrals over [a,b] with numerical accuracy of 1e-12 (on output)
   */
  public static void integrate(RealVectorValuedFunctionOfOneVariable f,
                               double[] values,
                               double a, double b) {

    integrate(f, a, b, values, 1e-12);
  }

  /**
   * Integrates set of real functions on interval with prescribed accuracy.
   * @param f set of real function
   * @param a start of interval
   * @param b end of interval
   * @param values of integrals over [a,b] (on output)
   * @param acc accuracy for result
   */
  public static void integrate(RealVectorValuedFunctionOfOneVariable f,
                               double a, double b,
                               double[] values,
                               double acc) {

    BulirschStoerIntegrator integrator
        = new BulirschStoerIntegrator(f);

    integrator.setEps(acc);

    integrator.integrate(a, b, values);
  }

}