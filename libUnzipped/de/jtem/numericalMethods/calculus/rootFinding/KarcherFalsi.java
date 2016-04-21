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
 * KarcherFalsi is an improved regula falsi for finding bracketed roots.
 * The main weakness of the regula falsi is that usually only one endpoint
 * of the bracketing interval is improved and that therefore during many
 * iterations interval halfing would have been better. 
 * KarcherFalsi counts how often the improved approximation of the root is
 * on the same side of the interval. A weight proportional to this count
 * moves the next regula falsi approximation towards the other endpoint.
 * To avoid overshooting the weight depends quadratically on the ratio
 * of the two most recent approximations where the function has the same
 * sign.
 * As few function values as possible are computed.
 * 
 * @author karcher
 *
 */
public class KarcherFalsi {
//
//	%% ....... Part 2,  The improved regula falsi ............. <<<============<<
//	%% THEORY for this program:
//	   % Standard regula falsi replaces a function (with different signs at the
//	   % endpoints of an intervall) by its secant, then takes the zero of the 
//	   % secant as the improving iteration. This has two weaknesses:
//	   % (1) Usually only *one* endpoint is improved under the iteration -- this
//	   % onesided convergence happens for example if the 2nd derivative does
//	   % NOT change sign.
//	   % (2) Often intervall-halfing is dramatically better than regula falsi,
//	   % for example for f(x) = exp(x) - 1, starting on [a,b] = [-12, +10].
//	   %
//	   % THIS  PROGRAM  avoids these weaknesses by computing a weight from the
//	   % ratio between the last function value and the other most recent value
//	   % with the same sign. With this weight the zero of the secant is moved
//	   % towards the endpoint which was NOT changed on the last iteration.
//	   %
//	   % Program details. COUNT keeps track of the bad cases, it says how often
//	   % the same side was improved in a row. INDEX is either 1 or 2;
//	   % arg(index) is the endpoint which was improved in the last round, and
//	   % arg(3-index) was not so recently improved.
//	   % WEIGHT  is *quadratic* in the ratio of function values: 
//	   %      (most recent value / previous value with same sign)^2, 
//	   % so that a significantly improved value results in only a *small* shift
//	   % of the zero of the secant. Also, weight < = 1, so that the zero of
//	   % the secant is at most shifted to the midpoint between this zero and
//	   % the unimproved endpoint.
//
//
//	if sign(leftValue) * sign(rightValue)  <= 0,                   % otherwise do nothing
//	   index = 1;     count = 1;    iteration_steps = 0;    % initialize
//	   
//	% ..... Start with one unmodified linear interpolation, detect endpoint-zeros:
//	if leftValue == 0, u = arg(1);       % needed for the case leftValue = rightValue = 0
//	else
//	 u=( arg(1)*rightValue - arg(2)*leftValue )/(rightValue - leftValue); % linear interpolation
//	end % if .. else
//	 v=feval(fct,u);  x0 = u;         % This is  OUTPUT  if abs(v) < error     
//
//	% .................. MAIN Iteration
//	while abs(v) > error                 %  IF abs(v) < error THEN  output = x0      
//	         % NOTE: v = 0 does NOT occur in while-loop, (termination if fct=0) 
//	     if sign(v) == sign(val(index)), index=index; count=count+1;   % BAD case
//	        %% arg(index) was improved in the last round, arg(3-index) was not
//		 else index=3-index; count=1;                                  % GOOD case
//	     end % if .. else
//	  weight = count*min(0.25,(v/val(index) )^2) ;            %%   The weight (##) 
//	  arg(index)=u; val(index)=v;                             % new endpoints   
//	  u=( arg(1)*rightValue - arg(2)*leftValue )/(rightValue - leftValue);  % zero of secant
//
//	    %% .......... The improved guess with the help of weight: .......... (##)
//	       u = ( u + weight*arg(3-index) )/(1 + weight); 
//	       v = feval(fct,u);   x0=u;  
//
//	    if gflag     % DEMO-graphics, this is SLOW if fct is slow:
//	     s = [0:20]/20;   x = arg(2)*s+arg(1)*(1-s);    y = feval(fct,x); 
//		 iteration_steps = iteration_steps+1;
//	     plot(x',y', [arg,u;arg,u], [0,0,0; val,v], arg',[0;0]), drawnow, pause(4)
//		end % if gflag
//	end % while             END of iteration, error bound achieved. Prepare output:
//	     x0 = u; value = v;
//	else  x0 = []; value = []; message = 'There are no zeros in this interval'                    
//	end % if sign(leftValue) * sign(rightValue)  <= 0   
//	%% ....... End of Part 2 

	final static int ITMAX = 10000;
	final static double EPS = 1e-15;

	public static int stepCount;
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
	
	static double sign( double v ) { return Math.signum(v); }
	static double abs( double v ) { return Math.abs(v); }

	public static double search(RealFunctionOfOneVariable f, double left,
			double right, double precision) {

		// a and b are the endponts of the bracketing interval
		// and stay not orderd during the algorithm
		double a = left;
		double b = right;
		double f_a = f.eval(a);
		double f_b = f.eval(b);
		double u = 0.0;
		double v = 0.0;

		if (sign(f_a) * sign(f_b) > 0) //TODO: discuss this 
			throw new IllegalArgumentException("root is not bracketed");

		int count = 1;
		int iteration_steps = 0;

		if (f_a==0.0)  {
			u = a;
			}
		else{
			u = (a * f_b - b * f_a) / (f_b - f_a);
			if (Math.abs(f_a)< 1e-13){ // Deals only with very small f_a
				double d = (b-a)/1024.0;
				a = a - d;
				b = a + d;
				f_a = f.eval(a);
				f_b = f.eval(b);
				v = f.eval(u);
				if (sign(f_a) * sign(v) > 0) {
                     if  (sign(f_a) * sign(f_b) > 0) //TODO: discuss this 
					 throw new IllegalArgumentException("root is not bracketed");
				     else u = (a * f_b - b * f_a) / (f_b - f_a); 
                     }
				else{
					u = (a * v - u * f_a) / (v - f_a); 
				}
			}  // small f_a is treated.
			

		//		 .................. MAIN Iteration
		while (abs(a - b) > 2 * precision ) {
			v = f.eval(u);
			if (v == 0.0){
				a = u; b = u;
			}
			else {
			if (sign(v) == sign(f_a)) { // BAD case, same endpoint a is improved
				count = count + 1;
			} 
			else { // GOOD case, count was not increased, the other endpoint b improved
				double tmp1 = f_a;
				f_a = f_b;
				f_b = tmp1;
				double tmp2 = a;
				a = b;
				b = tmp2;
				count = 1;
			}
			final double weight = count * Math.min(0.25, (v / f_a) * (v / f_a)); // The weight (##) 
			a = u;
			f_a = v;

			u = (a * f_b - b * f_a) / (f_b - f_a); // zero of secant

			//.......... The improved guess with the help of weight: .......... (##)
			u = (u + weight * b) / (1 + weight);
	
			iteration_steps++;
			} // else
		} // while             END of iteration, error bound achieved.
		
		} // if (!(f_a==0))
		stepCount = iteration_steps;
		return u;
	}

	
}
