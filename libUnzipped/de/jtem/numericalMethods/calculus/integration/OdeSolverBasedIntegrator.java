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

import java.io.Serializable;

import de.jtem.numericalMethods.calculus.function.RealFunctionOfOneVariable;
import de.jtem.numericalMethods.calculus.function.RealVectorValuedFunctionOfOneVariable;
import de.jtem.numericalMethods.calculus.odeSolving.ODE;
import de.jtem.numericalMethods.calculus.odeSolving.OdeSolver;

abstract class OdeSolverBasedIntegrator
	implements Serializable, RealVectorValuedFunctionIntegrator {

	private static final long serialVersionUID = 1L;

	private int numberOfFunctions;

	OdeSolver odeSolver;

	WraperFunction wf;

	double[] I;

	OdeSolverBasedIntegrator(OdeSolver odeSolver) {

		this.odeSolver = odeSolver;
	}

	public OdeSolver getOdeSolver() {
		return odeSolver;
	}

	private void changeNumberOfFunctions(int numberOfFunctions) {

		I = new double[numberOfFunctions];

		this.numberOfFunctions = numberOfFunctions;
	}

	public void setFunction(RealFunctionOfOneVariable f) {
		wf = new OdeSolverBasedIntegrator.WraperFunction.RealFunction(f);

		changeNumberOfFunctions(1);
	}

	public void setFunction(RealVectorValuedFunctionOfOneVariable f) {
		wf = new OdeSolverBasedIntegrator.WraperFunction.RealFunctionSet(f);

		changeNumberOfFunctions(f.getDimensionOfTargetSpace());
	}

	double startX;
	
	public void startAt(double a) {
		startX = a;

		for (int i = 0; i < I.length; i++) {
			I[i] = 0;
		}
	}

	public void integrateTo(double b) {

		if (wf == null) {
			throw new RuntimeException("no function to integrate");
		}

		odeSolver.setNumOfEquations(numberOfFunctions);
		odeSolver.odex(wf, I, startX, b);
		startX = b;
	}

	public void getIntegral(double[] integral, int offset) {

		for (int i = 0, j = offset; i < I.length; i++, j++) {
			integral[j] = I[i];
		}
	}

	public void integrate(double start, double end, double[] integral) {
		integrate(start, end, integral, 0);
	}

	public void integrate(
		double start,
		double end,
		double[] integral,
		int offset) {

		startAt(start);
		integrateTo(end);
		getIntegral(integral, offset);
	}

// TODO: remove 14.11.04
//	public double integrate(double start, double end, int index) {
//		startAt(start);
//		integrateTo(end);
//		return I[1 + index];
//	}

	public double integrate(double start, double end) {
		startAt(start);
		integrateTo(end);
		return I[0];
	}

	static abstract class WraperFunction
		implements ODE, Serializable, Cloneable {

		private static final long serialVersionUID = 1L;

		static class RealFunction extends WraperFunction {

			final RealFunctionOfOneVariable realFunctionOfOneVariable;

			RealFunction(RealFunctionOfOneVariable f) {

				realFunctionOfOneVariable = f;
			}

			public int getNumberOfEquations() {
				return 1;
			}
			
			public void eval( double t, double [] x, double y[] ) {
				y[0] = realFunctionOfOneVariable.eval(t);
			}
			
		}

		static class RealFunctionSet extends WraperFunction {

			final RealVectorValuedFunctionOfOneVariable realFunctionSetOfOneVariable;

			RealFunctionSet(RealVectorValuedFunctionOfOneVariable f) {

				realFunctionSetOfOneVariable = f;
			}

			public int getNumberOfEquations() {
				return realFunctionSetOfOneVariable.getDimensionOfTargetSpace();
			}
			
			public void eval( double t, double [] x, double y[] ) {
				realFunctionSetOfOneVariable.eval(t, y);
			}
			
		}

	}

}