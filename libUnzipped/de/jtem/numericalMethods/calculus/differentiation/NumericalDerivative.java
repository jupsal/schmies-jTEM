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

import java.util.HashMap;

import de.jtem.numericalMethods.calculus.function.RealFunctionOfOneVariable;
import de.jtem.numericalMethods.calculus.function.RealFunctionOfSeveralVariables;
import de.jtem.numericalMethods.calculus.function.RealFunctionOfSeveralVariablesWithGradient;
import de.jtem.numericalMethods.calculus.function.RealVectorValuedFunctionOfSeveralVariables;
import de.jtem.numericalMethods.calculus.function.RealVectorValuedFunctionOfSeveralVariablesWithJacobien;

/**
 * Computing derivative numericly.
 * 
 * @author schmies
 */
public class NumericalDerivative {

	
	NumericalDerivative() {
	}
	
	final static int DEFALUT_MAX_TABLE_LENGTH = 10;
	
	public static RealVectorValuedFunctionOfSeveralVariablesWithJacobien 
		createDerivativeNumerically(
			final RealVectorValuedFunctionOfSeveralVariables F, final double h ) {
		return createDerivativeNumerically( F, h, DEFALUT_MAX_TABLE_LENGTH, false );
	}
	
	public static RealVectorValuedFunctionOfSeveralVariablesWithJacobien 
		createDerivativeNumerically(
				final RealVectorValuedFunctionOfSeveralVariables F, final double h, final int maxTableLength, final boolean hashValues ) {
		return new   RealVectorValuedFunctionOfSeveralVariablesWithJacobien() {

			public void eval(double[] x, double[] values, int offset, double[][] jacobien) {
				F.eval( x, values, offset );
				NumericalDerivative.computeJacobienByRidders(F,x,jacobien,null,h, maxTableLength, hashValues);
			}

			public int getDimensionOfTargetSpace() {
				return F.getDimensionOfTargetSpace();
			}

			public int getNumberOfVariables() {
				return F.getNumberOfVariables();
			}

			public void eval(double[] x, double[] values, int offset) {
				F.eval( x, values, offset );
			}
			
		};
	}
	
	static class CoordinateVariableProjection implements RealFunctionOfOneVariable {

		final RealVectorValuedFunctionOfSeveralVariables F;
		
		final int numberOfVariables, dimensionOfTargetSpace;
		int coordinate = 0, variable = 0;
		
		final double [] x;
		final double [] values;
		
		final HashMap hashMap;
		
		CoordinateVariableProjection( 
				RealVectorValuedFunctionOfSeveralVariables F, 
				double [] x, boolean hashValues ) {
			
			this.F = F;
			this.x = x;
			
			hashMap = hashValues ? new HashMap() : null;
			
			numberOfVariables = F.getNumberOfVariables();
			dimensionOfTargetSpace = F.getDimensionOfTargetSpace();
			
			values = new double[dimensionOfTargetSpace];
		}
		
		public double eval(double t) {
			Double T = null;
			
			if( hashMap != null ) {
				T = new Double(t);			
				double[] hashedValues = (double[]) hashMap.get(T);
				if (hashedValues != null) {
					return hashedValues[coordinate];
				} 
			}
			
			double xOfVariable = x[variable];
			x[variable] += t;
			F.eval(x, values, 0);
			x[variable] = xOfVariable;
			
			if( hashMap != null ) {
				hashMap.put( T, values.clone() );
			}
			
			return values[coordinate];
		}
	}

	public static void computeJacobienByRidders( 
			RealVectorValuedFunctionOfSeveralVariables F, double [] x,
			double [][] jacobien, double [][] error, double h ) {
		computeJacobienByRidders( F, x, jacobien, error, h, DEFALUT_MAX_TABLE_LENGTH, false );
	}
	
	public static void computeJacobienByRidders( 
			RealVectorValuedFunctionOfSeveralVariables F, double [] x,
			double [][] jacobien, double [][] error, double h, int maxTableLength, boolean hashValues ) {
		
		final CoordinateVariableProjection cvp 
			= new CoordinateVariableProjection( F, x, hashValues );
		
		final double [] err = new double[1];
	
		for( cvp.variable=0; 
				cvp.variable<cvp.numberOfVariables; cvp.variable++ ) {
			    
			if(hashValues)
				cvp.hashMap.clear();
			    
			for( cvp.coordinate=0; 
				cvp.coordinate<cvp.dimensionOfTargetSpace; cvp.coordinate++ ) {
			
				jacobien[cvp.coordinate][cvp.variable] = Ridders.compute(cvp,0,h,err);
				
				if( error != null )
					error[cvp.coordinate][cvp.variable] = err[0];
			}
		}
	}
	
	

	static class VariableProjection implements RealFunctionOfOneVariable {

		final RealFunctionOfSeveralVariables F;
		
		final int numberOfVariables;
		int variable = 0;
		
		final double [] x;
		
		VariableProjection( 
				RealFunctionOfSeveralVariables F, 
				double [] x ) {
			
			this.F = F;
			this.x = x;
			
			numberOfVariables = F.getNumberOfVariables();
		}
		
		public double eval(double t) {
		
			double xOfVariable = x[variable];
			x[variable] += t;
			final double value = F.eval(x);
			x[variable] = xOfVariable;
			
			return value;
		}
	}
	

	public static RealFunctionOfSeveralVariablesWithGradient
	createDerivativeNumerically(
			final RealFunctionOfSeveralVariables F, final double h ) {
		return createDerivativeNumerically( F, h, DEFALUT_MAX_TABLE_LENGTH );
	}
	
	public static RealFunctionOfSeveralVariablesWithGradient 
	createDerivativeNumerically(
			final RealFunctionOfSeveralVariables F, final double h, final int maxTableLength ) {
		return new   RealFunctionOfSeveralVariablesWithGradient() {

			public double eval(double[] x, double [] grad ) {
				NumericalDerivative.computeGradientByRidders(F,x,grad,null,h, maxTableLength);
				return F.eval( x );
			}

			public int getNumberOfVariables() {
				return F.getNumberOfVariables();
			}

			public double eval(double[] x) {
				return F.eval( x );
			}
		};
	}
	
	public static void computeGradientByRidders( 
			RealFunctionOfSeveralVariables F, double [] x,
			double [] grad, double [] error, double h ) {
		computeGradientByRidders( F, x, grad, error, h, DEFALUT_MAX_TABLE_LENGTH );
	}
	
	public static void computeGradientByRidders( 
			RealFunctionOfSeveralVariables F, double [] x,
			double [] grad, double [] error, double h, int maxTableLength ) {
		
		final VariableProjection cvp 
		= new VariableProjection( F, x );
		
		final double [] err = new double[1];
		
		for( cvp.variable=0; 
		cvp.variable<cvp.numberOfVariables; cvp.variable++ ) {
			
			
			grad[cvp.variable] = Ridders.compute(cvp,0,h,err);
			
			if( error != null )
				error[cvp.variable] = err[0];
		}
	}
}
