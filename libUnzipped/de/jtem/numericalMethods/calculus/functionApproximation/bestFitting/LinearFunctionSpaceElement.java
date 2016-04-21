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

/**
 * @author schmies
 *
 */
public class LinearFunctionSpaceElement {

		final LinearFunctionSpace space;
		final double [] coefficient;
		
		public LinearFunctionSpaceElement( LinearFunctionSpace space ) {
			
			this.space = space;
			this.coefficient =  new double[ space.getDimension() ];
		}
		public LinearFunctionSpaceElement( LinearFunctionSpace space,
		final double [] parameter ) {
			
			this.space = space;
			this.coefficient = (double[])parameter.clone();
		}
		
		public void eval( double [] x, double [] y ) {
			space.eval( x, coefficient, y );
		}
		
		public double [] getCoefficients() {
			return (double[])coefficient.clone();
		}
		
		public LinearFunctionSpace getApproximationFunctionSpace() {
			return space;
		}
		
		final static void minus( double [] a, double [] b, double [] c ) {
			for(int i=0; i<a.length; i++)
				c[i] = a[i] - b[i];
		}
		
		final static double absSqr( final double [] u ) {
			double uu=0;
			for(int i=0; i<u.length; i++)
				uu+=u[i]*u[i];
			return uu;
		}
		
		public Sample [] getDeviationSamples( Sample [] samples ) {
				final int nos=samples.length;
				Sample [] deviationSample = new Sample[ nos ];
				double [] y = new double[space.m];
				for( int i=0; i<nos; i++ ) {
					eval( samples[i].x, y );
					minus( samples[i].y, y, y);
					deviationSample[i] = new Sample(samples[i].x, y, samples[i].weight );
				}
				return deviationSample;
		}
		
		public double getAverageDeviation( Sample [] samples ) {
			final int nos=samples.length;
			if( nos == 0 )
					return 0;
			double [] y = new double[space.m];
			double dev=0;
			for( int i=0; i<nos; i++ ) {
				eval( samples[i].x, y );
				minus( y, samples[i].y, y);
				dev += Math.sqrt(absSqr( y ));
			}
			return dev / nos;
		}
		
		public double getMaximalDeviation( Sample [] samples ) {
			final int nos=samples.length;
			if( nos == 0 )
				return 0;
			double [] y = new double[space.m];
			double max=0;
			for( int i=0; i<nos; i++ ) {
				eval( samples[i].x, y );
				minus( y, samples[i].y, y);
				double dev = Math.sqrt(absSqr( y ));
				if( max<dev )
						max = dev;
			}
			return max;
		}
		
		public StringBuffer componentString( int index ) {
			StringBuffer sb=new StringBuffer(300);
			boolean first = true;
			for( int i=0; i<space.dim; i++) {
				double coeff = coefficient[i*space.m + index];
				if( coeff == 0 || space.f[i][index].isZero() )
					continue;
				if( first ) {
					first = false;
				} else {
					sb.append( " +");
				}
				sb.append( " ( "+coefficient[i*space.m + index]+") * ( " + space.f[i][index].toString() + " )");
			}
			
			if( first ) sb.append( " 0 ");
			
			return sb;
		}
		
		public String toString()  {
			StringBuffer sb=new StringBuffer(300);
			sb.append("{\n");
			for( int i=0; i<space.m; i++ ) {
				sb.append( "\t y["+i+"] = "+componentString(i));
				sb.append(";\n");
			}
			sb.append("}\n");
			return sb.toString();
		}
		
		public LinearFunctionSpaceElement getPartialDerivative( int k ) {
			return new LinearFunctionSpaceElement( 
					LinearFunctionSpaceFactory.createPartialDerivative(space,k), coefficient );
		}
}
