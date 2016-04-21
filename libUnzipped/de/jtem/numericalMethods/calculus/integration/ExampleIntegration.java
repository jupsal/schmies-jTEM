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

/**
 * This classs illustrates the usage of the numerical integrators,
 * that are based on adaptive ode solvers.
 * <p>
 * We integrate first
 * </p><p>
 * f(t) = sin(t) ,
 * </p><p>
 * which is defined using the interface {@link RealFunctionOfOneVariable}:
 * </p>
 * <pre>
 *  // creating function f(t) = sin(t)
 *  RealFunctionOfOneVariable f = new RealFunctionOfOneVariable() {
 *     public double eval(double t) {
 *       return Math.sin(t);
 *     }
 *   };
 * </pre>
 * To integrate the function you can call a static function on the
 * several integrator classes, e.g.:
 * <pre>
 * ExtrapIntegrator.integrate(g, 0, Math.PI, resultOfExtrap, acc);
 * </pre>
 * or create an instance of {@link ExtrapIntegrator} and set the
 * relevant data as properties.
 * </p><p>
 * Further we integrate the vector valued function
 * </p><p>
 * g(t) = ( cos(t), sin(t) )
 * </p><p>
 * which is defined using the interface {@link RealFunctionSetOfOneVariable}:
 * </p>
 * <pre>
 *  // creating function g(t) = ( sin(t), cos(t) )
 *  RealFunctionSetOfOneVariable g = new RealFunctionSetOfOneVariable() {
 *    public int getNumberOfFunctions() {
 *      return 2;
 *    }
 *    public void eval(double t, double[] values, int offset) {
 *      values[offset    ] = Math.sin(t);
 *      values[offset + 1] = Math.cos(t);
 *    }
 *  };
 * </pre>
 * The handling for the vector valued function is simular to the one
 * dimensional case:
 * <pre>
 *    double[][] result = new double[2];
 *    ExtrapIntegrator.integrate(g, 0, Math.PI, result[0], acc);
 * </pre>
 *
 * @see RealFunctionOfOneVariable
 * @see RealVectorValuedFunctionOfOneVariable
 * @see ExtrapIntegrator
 * @see RungeKuttaFehlbergIntegrator
 * @see BulirschStoerIntegrator
 *
 * @author Markus Schmies
 *
 * @version 0.9
 */
public class ExampleIntegration {

  /**
   * Integrates the above given functions.
   * The application uses three different
   * adaptive integrators based on the ode solvers,
   * Extrap, RungeKuttaFehlberg, and BulirschStoer and
   * compares the results for different accuracies.
   * @param arg are not used
   */
  static public void main(String[] arg) {

    // creating function f(t) = sin(t)
    RealFunctionOfOneVariable f = new RealFunctionOfOneVariable() {
      public double eval(double t) {
        return Math.sin(t);
      }
    };

    System.out.println(
        "integrating f(t) = sin(t) on interval [0,PI]; result is 2");

    // integrate with different accuracies
    for (int i = 3; i < 15; i++ ) {
      double acc = Math.pow( 10, -i );
      System.out.println("setting accuracy to " + acc);

      double resultOfExtrap
          = ExtrapIntegrator.integrate(f, 0, Math.PI, acc);
      double resultOfBulrischStoer
          = BulirschStoerIntegrator.integrate(f, 0, Math.PI, acc);
      double resultOfRungeKuttaFehlberg
          = RungeKuttaFehlbergIntegrator.integrate(f, 0, Math.PI, acc);

      System.out.println("result of Extrap:             " + resultOfExtrap);
      System.out.println("result of BulirschStoer:      " + resultOfBulrischStoer);
      System.out.println("result of RungeKuttaFehlberg: " + resultOfRungeKuttaFehlberg);
      System.out.println("___");
    }

    // creating function g(t) = ( sin(t), cos(t) )
    RealVectorValuedFunctionOfOneVariable g = new RealVectorValuedFunctionOfOneVariable() {
      public int getDimensionOfTargetSpace() {
        return 2;
      }

      public void eval(double t, double[] values ) {
        values[0] = Math.sin(t);
        values[1] = Math.cos(t);
      }
    };

    System.out.println
        ("integrating g(t) = (sin(t),cos(t) on interval [0,PI]; result is (2,0)");

    double[][] result = new double[3][2];

    // integrate with different accuracies
    for (int i = 3; i < 15; i++  ) {
      double acc = Math.pow( 10, -i );

      System.out.println("setting accuracy to " + acc);

      ExtrapIntegrator
          .integrate(result[0], g, 0, Math.PI, acc);
      BulirschStoerIntegrator
          .integrate(g, 0, Math.PI, result[1], acc);
      RungeKuttaFehlbergIntegrator
          .integrate(g, 0, Math.PI, result[2], acc);

      System.out.println("result of Extrap:             ( " +
                         result[0][0] + " , " + result[0][1] + " ) ");
      System.out.println("result of BulirschStoer:      ( " +
                         result[1][0] + " , " + result[1][1] + " ) ");
      System.out.println("result of RungeKuttaFehlberg: ( " +
                         result[2][0] + " , " + result[2][1] + " ) ");
      System.out.println("___");

    }

  }

}