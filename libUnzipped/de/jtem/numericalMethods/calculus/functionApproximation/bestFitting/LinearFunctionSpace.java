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

package de.jtem.numericalMethods.calculus.functionApproximation.bestFitting;


import de.jtem.numericalMethods.calculus.function.RealFunctionOfSeveralVariablesWithGradient;
import de.jtem.numericalMethods.calculus.minimizing.ConjugateGradient;

/**
 * @author schmies
 *
 * To change the template for this generated type comment go to
 * Window - Preferences - Java - Code Generation - Code and Comments
 */
public class LinearFunctionSpace {

	BasisFunction [][] f;
	
	public final int n, m, dim;

	double [][] c;

	LinearFunctionSpace( BasisFunction [][] f ) {
		this.f =f;
		
		n = f[0][0].getNumberOfVariables();
		m = f[0].length;
		dim =f.length;
		
		c =new double[dim][m];
	}
		
	public int getDimension() {
		return m*dim;
	}
	
	void eval( double [] x, double [] parameter, double [] u ) {
		eval( x, parameter, u, c );
	}
	void eval( double [] x, double [] parameter, double [] u, double [][] c ) {

		for(int j=0; j<m; j++ ) {
			u[j]=0;
		}
		
		for( int i=0, k=0; i<dim; i++ ) {
			for(int j=0; j<m; j++,k++ ) {
				final double f_ij= f[i][j].eval(x);
				c[i][j]= f_ij;
				u[j] += parameter[k] * f_ij;
				
			}
		}
			
	}
	
	class Energy implements RealFunctionOfSeveralVariablesWithGradient {

		final Sample [] samples;
		
		final double [] y;
		
		final int numberOfDimensionBeingFixed;
		
		Energy( Sample [] samples,  int numberOfDimensionBeingFixed) {
			this.samples = samples;
			this.numberOfDimensionBeingFixed = numberOfDimensionBeingFixed;
			y =new double[m];
		}
		Energy( Sample [] samples ) {
			this( samples, 0 );
		}
		
		final double absSqr( double [] u ) {
			double uu=0;
			for(int i=0; i<m; i++)
				uu+=u[i]*u[i];
			return uu;
		}
	
		/* (non-Javadoc)
		 * @see de.jtem.numericalMethods.calculus.minimizing.RealFunctionOfSeveralVariablesWithGradient#eval(double[], double[])
		 */
		public double eval( final double [] parameter, final double [] gradient ) {
			double energy = 0;
			if( gradient != null ) {
				for( int j=0; j<m*dim; j++ ) {
					gradient[j]=0;
				}
			}
			for( int i=0; i<samples.length; i++ ) {
				final Sample sample = samples[i];
				
				LinearFunctionSpace.this.eval( sample.x, parameter, y, c );
				
				for( int j=0; j<m; j++) {
					y[j] -= sample.y[j];
				}
				if( gradient != null ) {
					for( int j=numberOfDimensionBeingFixed, 
							l=m*numberOfDimensionBeingFixed; j<dim; j++ ) {
						for( int k=0; k<m; k++, l++) {
							gradient[l] += c[j][k] * y[k] * sample.weight;
						}
					}
				}
				energy += absSqr(y) * sample.weight;
			}
			return energy/2;
		}

		public double eval(double[] parameter) {
			return eval( parameter, null );
		}

		public int getNumberOfVariables() {
			return getDimension();
		}
		
	}
	
	/**
	 * 
	 * @param f
	 * @param samples
	 * @param tol
	 * @param maxNumberOfSteps
	 */
	public LinearFunctionSpaceElement searchForBestFit( double [] startParameter, Sample [] samples, double tol, int maxNumberOfSteps ) {
		
		LinearFunctionSpaceElement best = new  LinearFunctionSpaceElement(this);
		
		Energy energy = new Energy( samples );
		
		ConjugateGradient.search( best.coefficient, tol, energy, maxNumberOfSteps, false, null );
		
		return best;
	}
	
	public LinearFunctionSpaceElement searchForBestFit( 
			LinearFunctionSpaceElement start, Sample [] samples, double tol, int maxNumberOfSteps ) {
		return searchForBestFit( start, samples, tol, maxNumberOfSteps, 0 );
	}
	public LinearFunctionSpaceElement searchForBestFit( 
			LinearFunctionSpaceElement start, Sample [] samples, double tol, int maxNumberOfSteps, int numberOfDimensionsBeingFixed ) {
		
		LinearFunctionSpaceElement best = new  LinearFunctionSpaceElement(this, start.getCoefficients());
		
		Energy energy = new Energy( samples, numberOfDimensionsBeingFixed  );
		
		ConjugateGradient.search( best.coefficient, tol, energy, maxNumberOfSteps, false, null );
		
		return best;
	}
	
	public LinearFunctionSpaceElement searchForBestFit( Sample [] samples, double tol, int maxNumberOfSteps ) {
		return searchForBestFit(new  LinearFunctionSpaceElement(this), samples, tol, maxNumberOfSteps );
	}
	
	
}
