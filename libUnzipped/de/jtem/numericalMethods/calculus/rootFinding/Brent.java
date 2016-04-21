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

import de.jtem.numericalMethods.calculus.function.RealFunctionOfOneVariable;

/**
 * Searches the root of a one dimensional real function using Brent's
 * algorithm.
 * 
 * @author schmies
 */
public class Brent {

	final static int ITMAX = 100;
	final static double EPS = 1e-15;

	/**
	 * Searchs the root of function f in brackets [left,right] with
	 * prescribed precision. 
	 * left and right must really bracket a root, e.g. f(left)*f(right)<0.
	 * @param f real function
	 * @param left bracket
	 * @param right bracket
	 * @param precision of root
	 * @return a root of x
	 */
	public static double search( RealFunctionOfOneVariable f, 
			final double left, double right, double precision ) {
		
		double a=left;
		double b=right;
		double c=right;

		double fa=f.eval(a);
		double fb=f.eval(b);

		if ((fa > 0.0 && fb > 0.0) || (fa < 0.0 && fb < 0.0))
			throw new IllegalArgumentException( "root is not bracketed");
		double fc=fb;
		
		double d=0, e=0;
		for (int iter=1;iter<=ITMAX;iter++) {
			if ((fb > 0.0 && fc > 0.0) || (fb < 0.0 && fc < 0.0)) {
				c=a;
				fc=fa;
				e=d=b-a;
			}
			if (Math.abs(fc) < Math.abs(fb)) {
				a=b;
				b=c;
				c=a;
				fa=fb;
				fb=fc;
				fc=fa;
			}
			double tol1=2.0*EPS*Math.abs(b)+0.5*precision;
			double xm=0.5*(c-b);
			if (Math.abs(xm) <= tol1 || fb == 0.0) return b;
			if (Math.abs(e) >= tol1 && Math.abs(fa) > Math.abs(fb)) {
				final double s=fb/fa;
				double p, q, r;
				if (a == c) {
					p=2.0*xm*s;
					q=1.0-s;
				} else {
					q=fa/fc;
					r=fb/fc;
					p=s*(2.0*xm*q*(q-r)-(b-a)*(r-1.0));
					q=(q-1.0)*(r-1.0)*(s-1.0);
				}
				if (p > 0.0) q = -q;
				p=Math.abs(p);
				final double min1=3.0*xm*q-Math.abs(tol1*q);
				final double min2=Math.abs(e*q);
				if (2.0*p < (min1 < min2 ? min1 : min2)) {
					e=d;
					d=p/q;
				} else {
					d=xm;
					e=d;
				}
			} else {
				d=xm;
				e=d;
			}
			a=b;
			fa=fb;
			if (Math.abs(d) > tol1)
				b += d;
			else
				b += xm >= 0 ? Math.abs(tol1) : -Math.abs(tol1);
			fb=f.eval(b);
		}
		throw new RuntimeException( "exceeded maximal number of iterations");
	}

}
