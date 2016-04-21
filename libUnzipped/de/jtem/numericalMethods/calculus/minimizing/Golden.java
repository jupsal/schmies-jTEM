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

package de.jtem.numericalMethods.calculus.minimizing;

import de.jtem.numericalMethods.calculus.function.RealFunctionOfOneVariable;




/**This class represents the routine to find a local minimum of the function f in a bracking triplet(a,b,c).
 * The triplet must be of the constuction, where middle point of tiplet must lies between first and third point
 * ( Such f(b) < f(a),f(c) ).
 * @author Markus Schmies, Vitali Lieder
 * @version 1.0
 */
public final class Golden
    implements java.io.Serializable {

  private static final long serialVersionUID = 1L;

  static final double R = 0.61803399;
  static final double C = (1.0 - R);


  /**
   * Search the minimum of function f in the interval (ax,cx) with golden section search
   * and write it in the X (such X must be a double array of dim. 2).
   * @param ax first component of triplet.
   * @param bx second component of triplet.
   * @param cx third component of triplet.
   * @param X X[0] postion and X[1] value of minimum on output.
   * @param f given function
   * @param tol tolerance of calculation
   */
  public static final void search(double ax, double bx, double cx, double[] X,
                                  RealFunctionOfOneVariable f, double tol) {

    double f0, f1, f2, f3, x0, x1, x2, x3, xmin;

    x0 = ax;
    x3 = cx;
    if (Math.abs(cx - bx) > Math.abs(bx - ax)) {
      x1 = bx;
      x2 = bx + C * (cx - bx);
    }
    else {
      x2 = bx;
      x1 = bx - C * (bx - ax);
    }


    f1 = f.eval(x1);
    f2 = f.eval(x2);

    while (Math.abs(x3 - x0) > tol * (Math.abs(x1) + Math.abs(x2))) {
      if (f2 < f1) {
        x0 = x1;
        x1 = x2;
        x2 = R * x1 + C * x3; //SHFT(x0,x1,x2,R*x1+C*x3)

        f0 = f1;
        f1 = f2;
        f2 = f.eval(x2); //SHFT(f0,f1,f2,(*f)(x2))
      }
      else {
        x3 = x2;
        x2 = x1;
        x1 = R * x2 + C * x0; //SHFT(x3,x2,x1,R*x2+C*x0)

        f3 = f2;
        f2 = f1;
        f1 = f.eval(x1); //SHFT(f3,f2,f1,(*f)(x1))
      }
    }
    if (f1 < f2) {
      X[0] = x1;
      X[1] = f1;
    }
    else {
      X[0] = x2;
      X[1] = f2;
    }
  }

}






