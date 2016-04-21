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


/**
 * This class represents the routine to find a braket minimum of the function f.
 * After searching the new intervals with triplet will be saved in given array,
 * so that a < b < c or a > b > c and f(a),f(c) > f(b).
 * @author Markus Schmies, Vitali Lieder
 * @version 1.0
 */
public final class Braket implements java.io.Serializable {

    private static final long serialVersionUID = 1L;

    static final double GOLD   = 1.618034;
    static final double GLIMIT = 100.0;
    static final double TINY   = 1.0e-20;

    /**
     * Default constructor of the class.
     */
    public Braket(){}

    private static final double sign( double a, double b ) {
        return b > 0 ? Math.abs( a ) : - Math.abs( a );
    }

    /**
     * Given distinct initial interval a and b ( t[0] and t[1] ) and function f.
     * Routine searches in downhill direction the braket minimum of function f and save the new points
     * in triplet t ( a,b,c in t[0],t[1],t[2] ), also the values of function f are saved in fOfT.
     * @param t the initial point a (t[0]) and b(t[1]). The braket point c will be saved in t[2].
     * @param fOfT value of the function f at t.
     * @param f given function.
     */
    public static void search(double[] t, double[] fOfT,
                              RealFunctionOfOneVariable f) {
      search(t, fOfT, f, null);
    }


    /**
     * Given distinct initial interval a and b ( t[0] and t[1] ) and function f.
     * Routine searches in downhill direction the braket minimum of function f and save the new points
     * in triplet t ( a,b,c in t[0],t[1],t[2] ), also the values of function f are saved in fOfT.
     * @param t the initial point a (t[0]) and b(t[1]). The braket point c will be saved in t[2].
     * @param fOfT value of the function f at t.
     * @param f given function.
     * @param info some degug output.
     */
    public static void search(double[] t, double[] fOfT,
                              RealFunctionOfOneVariable f, final Info info) {

      double ulim, u, r, q, fu, dum = 0;

      double ax = t[0];
      double bx = t[1];

      double fa = f.eval(ax);

      double fb = f.eval(bx);

      double inValue = fb;

      if (fb > fa) {
        dum = ax;
        ax = bx;
        bx = dum;
        dum = fa;
        fa = fb;
        fb = dum;
      }

      double cx = (bx) + GOLD * (bx - ax);

      double fc = f.eval(cx);

      while (fb > fc) {
        r = (bx - ax) * (fb - fc);
        q = (bx - cx) * (fb - fa);
        u = (bx) - ( (bx - cx) * q - (bx - ax) * r) /
            (2.0 * sign(Math.max(Math.abs(q - r), TINY), q - r));
        ulim = (bx) + GLIMIT * (cx - bx);
        if ( (bx - u) * (u - cx) > 0.0) {

          fu = f.eval(u);

          if (fu < fc) {
            ax = bx;
            bx = u;
            fa = fb;
            fb = fu;
            break;
          }
          else if (fu > fb) {
            cx = u;
            fc = fu;
            break;
          }

          u = (cx) + GOLD * (cx - bx);

          fu = f.eval(u);

        }
        else if ( (cx - u) * (u - ulim) > 0.0) {

          fu = f.eval(u);

          if (fu < fc) {
            bx = cx;
            cx = u;
            u = cx + GOLD * (cx - bx);

            fb = fc;
            fc = fu;
            fu = f.eval(u);
          }
        }
        else if ( (u - ulim) * (ulim - cx) >= 0.0) {

          u = ulim;

          fu = f.eval(u);

        }
        else {

          u = (cx) + GOLD * (cx - bx);

          fu = f.eval(u);
        }
        ax = bx;
        bx = cx;
        cx = u;
        fa = fb;
        fb = fc;
        fc = fu;
      }

      t[0] = ax;
      fOfT[0] = fa;
      t[1] = bx;
      fOfT[1] = fb;
      t[2] = cx;
      fOfT[2] = fc;

      if ( ( fb > inValue ) && ( info !=null ) ){
        info.setMessage("proc Braket failed to decrease center value! ");
        info.printDebug();
      }
    }
}





