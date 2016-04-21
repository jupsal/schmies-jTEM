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

package de.jtem.numericalMethods.calculus.rootFinding;

/**
 * Example illustrates usage of {@link NewtonRaphson}.
 * @author schmies
 */
public class ExampleNewtonRaphson {

	public static void main(String[] argv) {

		NewtonRaphson.RealFunctionWithDerivative rf = new NewtonRaphson.RealFunctionWithDerivative() {
			public void eval(double x, double[] f, int offsetF, double[] df, int offsetDF) {
					f[offsetF] = Math.sin(x);
					df[offsetDF] = Math.cos(x);
			}
		};

		double[] rootValues = new double[6];

		NewtonRaphson.search(rf, 0.6, rootValues);

		System.out.println(" root at x = " + rootValues[0]);
		System.out.println(" sin ( x ) = " + rootValues[1]);
		System.out.println(" sin'( x ) = " + rootValues[2]);

		NewtonRaphson.ComplexFunctionWithDerivative cf = new NewtonRaphson.ComplexFunctionWithDerivative() {
			public void eval(double x, double y, double[] f, int offsetF, double[] df, int offsetDF) {
				double expY = Math.exp(y);
				  	    double cosX = Math.cos(x);
				  	    double sinX = Math.sin(x);
				  
				  	    f[ offsetF       ] = sinX * ( expY + 1/expY ) / 2;
				  	    f[ offsetF + 1 ] = cosX * ( expY - 1/expY ) / 2;
				  
				  	    df[ offsetDF ] = cosX * ( expY + 1/expY ) / 2;
				  	    df[ offsetDF + 1] = sinX * ( expY - 1/expY ) / 2;
			}
		};

		NewtonRaphson.search(cf, 0.5, 0.3, rootValues);

		System.out.println(" root at z = " + rootValues[0] + " + i " + rootValues[1]);
		System.out.println(" sin ( z ) = " + rootValues[2] + " + i " + rootValues[3]);
		System.out.println(" sin'( z ) = " + rootValues[4] + " + i " + rootValues[5]);

	}
}
