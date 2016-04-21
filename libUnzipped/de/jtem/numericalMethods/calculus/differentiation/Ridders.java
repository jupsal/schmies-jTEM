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

import de.jtem.numericalMethods.calculus.function.RealFunctionOfOneVariable;

/**
 * Ridders' method of polynomial extrapolation to compute 
 * derivative numerically.
 * 
 * @author schmies
 *
 */
public class Ridders {


	Ridders() {
	}



	final static double CON = 1.4;
	final static double CON2 =(CON*CON);
	final static double SAFE = 2.0;

	/**
	 * Computes the derivative of f at x by Ridders' method of
	 * polynomial extrapolation. The value h the estimated initial stepsize,
	 * which need not to be small, but should be an increment such that f
	 * substantially chanes. It uses a maximal table length of 10
	 * @param f function to deferentiate
	 * @param x position
	 * @param h initial stepsize
	 * @param error on output, first entry, may be null
	 * @return the derivative of f at x by Ridders' method
	 */
	public static double compute( RealFunctionOfOneVariable f, double x, double h, double [] error ) {
		return compute( f, x, h, 10, error, 0 );
	}
	
	/**
	 * Computes the derivative of f at x by Ridders' method of
	 * polynomial extrapolation. The value h the estimated initial stepsize,
	 * which need not to be small, but should be an increment such that f
	 * substantially chanes.
	 * @param f function to deferentiate
	 * @param x position
	 * @param h initial stepsize
	 * @param maxTableLength maximal length of extrapolation table
	 * @param error on output, may be null
	 * @param offset for the error value in the associated array
	 * @return the derivative of f at x by Ridders' method
	 */
	public static double compute( RealFunctionOfOneVariable f, double x, double h, int maxTableLength, double [] error, int errorOffset ) {
		
		double [][] table = new double[maxTableLength][maxTableLength];

		if (h == 0.0) 
			throw new IllegalArgumentException("h must be nonzero.");
		
		double df=Double.NaN, hh=h;
		
		table[0][0]=(f.eval(x+hh)-f.eval(x-hh))/(2.0*hh);
		
		double err=Double.MAX_VALUE;
		
		for (int i =1;i<maxTableLength;i++) {
			hh /= CON;
			table[0][i]=(f.eval(x+hh)-f.eval(x-hh))/(2.0*hh);
			double fac=CON2;
			for ( int j=1;j<i;j++) {
				table[j][i]=(table[j-1][i]*fac-table[j-1][i-1])/(fac-1.0);
				fac=CON2*fac;
				double errt=Math.max( Math.abs(table[j][i]-table[j-1][i]),Math.abs(table[j][i]-table[j-1][i-1]));
				if (errt <= err) {
					err=errt;
					df=table[j][i];
				}
			}
			if (Math.abs(table[i][i]-table[i-1][i-1]) >= SAFE*(err)) {
				if( error != null ) 
					error[0] = err;
				return df;
			}
		}
		if( error != null ) 
			error[errorOffset] = err;
		return df;
	}
}
