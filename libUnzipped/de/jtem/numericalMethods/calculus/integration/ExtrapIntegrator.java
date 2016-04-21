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
import de.jtem.numericalMethods.calculus.odeSolving.Extrap;

/**
 * Integrator based on the adaptive ode-solver
 * {@link de.jtem.numericalMethods.calculus.odeSolving.Extrap} .
 * For usage consult example class: {@link ExampleIntegration}.
 *
 * @see de.jtem.numericalMethods.calculus.odeSolving.Extrap
 * @see ExampleIntegration
 * @author Markus Schmies
 * @version 1.0
 */
public final class ExtrapIntegrator
    extends OdeSolverBasedIntegrator {

  private static final long serialVersionUID = 1L;

  double initialStepSize = 0.1;

  double absTol = 1e-7;
  double relTol = 1e-7;

  public ExtrapIntegrator() {
    super(new Extrap());
  }

  public ExtrapIntegrator( RealFunctionOfOneVariable f ) {
    super(new Extrap(1));
    setFunction( f );
  }

  public ExtrapIntegrator( RealVectorValuedFunctionOfOneVariable f ) {
    super(new Extrap(f.getDimensionOfTargetSpace()));
    setFunction( f );
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

  public void setEps(double tol) {

  	( (Extrap) odeSolver).setAbsTol(tol);
  	( (Extrap) odeSolver).setRelTol(tol);
  	
  	this.absTol = tol;
  	this.relTol = tol;
  }
  
  public double getAbsTol() {
    return absTol;
  }

  public void setAbsTol(double absTol) {

    ( (Extrap) odeSolver).setAbsTol(absTol);

    if (this.absTol == absTol) {
      return;
    }

    this.absTol = absTol;
  }

  public double getRelTol() {
    return relTol;
  }

  public void setRelTol(double relTol) {

    ( (Extrap) odeSolver).setRelTol(relTol);

    if (this.relTol == relTol) {
      return;
    }

    this.relTol = relTol;
  }

  private boolean justStarted = false;

  public void startAt(double a) {

    super.startAt(a);

    justStarted = true;
  }

  public void integrateTo(double b) {

    if (justStarted) {
      ( (Extrap) odeSolver).setCurrentStepSize(initialStepSize);
    }

    super.integrateTo(b);

    justStarted = false;
  }

  /**
   * Integrates real function on interval.
   * @param f function
   * @param a start of interval
   * @param b end of interval
   * @return integral of f over [a,b] with numerical accuracy of 1e-14
   */
  public static double integrate(RealFunctionOfOneVariable f,
                                 double a, double b) {

    return integrate(f, a, b, 1e-14);
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

    ExtrapIntegrator integrator = new ExtrapIntegrator(f);

    integrator.setAbsTol(acc);
    integrator.setRelTol(acc);

    return integrator.integrate(a, b);
  }

  /**
   * Integrates set of real function on interval.
   * @param f set of real function
   * @param a start of interval
   * @param b end of interval
   * @param integrals over [a,b] with numerical accuracy of 1e-14 (on output)
   */
  public static void integrate(RealVectorValuedFunctionOfOneVariable f,
                               double[] values,
                               double a, double b) {

    integrate(values, f, a, b, 1e-14);
  }

  /**
   * Integrates set of real functions on interval with prescribed accuracy.
   * @param values of integrals over [a,b] (on output)
   * @param f set of real function
   * @param a start of interval
   * @param b end of interval
   * @param acc accuracy for result
   */
  public static void integrate(double[] values,
                               RealVectorValuedFunctionOfOneVariable f, double a,
                               double b,
                               double acc) {

    ExtrapIntegrator integrator
        = new ExtrapIntegrator( f );

    integrator.setAbsTol(acc);
    integrator.setRelTol(acc);

    integrator.integrate(a, b, values);
  }
}