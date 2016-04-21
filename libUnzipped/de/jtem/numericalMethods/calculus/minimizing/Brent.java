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
 * The triplet must be of the constuction, where middle point of tiplet lies between first and third point
 * ( Such f(b) < f(a),f(c) ) .
 * @author Markus Schmies, Vitali Lieder
 * @version 1.0
 */
public final class Brent implements java.io.Serializable {

    private static final long serialVersionUID = 1L;

    static final double CGOLD = 0.3819660;
    static final double ZEPS  = 1.0e-17;

    /**
     * Maximum number of evalutaions in search methods.
     */
    static int ITMAX = 100;

    /**
     * Default constructor of the class.
     */
    public Brent(){}


    /**
       * Get the value of ITMAX.
       * @return Value of ITMAX.
       */
    public static double getITMAX() {return ITMAX;}

    /**
       * Set the value of ITMAX.
       * @param v  Value to assign to ITMAX.
       */
    public static void setITMAX(int  v) {ITMAX = v;}

    private static final double sign( double a, double b ) {
        return b > 0 ? Math.abs( a ) : - Math.abs( a );
    }

    /**
     * Search the minimum of function f in the interval (t[0],t[2]) with parabolic interpolation
     * and Brent's method and write it in the X (such X must be a double array of dim. 2).
     * @param t bracketing triplet of function f with t[0]<t[1]<t[2] or t[0]>t[1]>t[2]
     * and f(t[0]),f(t[2])>f(t[1]).
     * @param X the minimum point of function f in the interval (t[0],t[2]).X[0] is the position
     *  and X[1] is the value of minimum.
     * @param f given function.
     * @param tol the fractional precission.
     */
    public static final void search(double[] t, double[] X, RealFunctionOfOneVariable f,
                                    double tol) {

      search(t[0], t[1], t[2], X, f, tol);
    }

    /**
    * Search the minimum of function f in the interval (ax,cx) with parabolic interpolation
    * and Brent's method and write it in the X (such X must be a double array of dim. 2).
    * @param ax first point of bracketing triplet
    * @param bx second point of bracketing triplet
    * @param cx third point of bracketing triplet
    * @param X the minimum point of function f in the interval (ax,cx).X[0] is the position
    *  and X[1] is the value of minimum.
    * @param f given function.
    * @param tol the fractional precission.
    */
   public static final void search(double ax, double bx, double cx, double[] X,
                                   RealFunctionOfOneVariable f, double tol) {
     search(ax,bx,cx,X,f,tol,null);
   }


    /**
     * Search the minimum of function f in the interval (ax,cx) with parabolic interpolation
     * and Brent's method and write it in the X (such X must be a double array of dim. 2).
     * @param ax first point of bracketing triplet
     * @param bx second point of bracketing triplet
     * @param cx third point of bracketing triplet
     * @param X the minimum point of function f in the interval (ax,cx).X[0] is the position
     *  and X[1] is the value of minimum.
     * @param f given function.
     * @param tol the fractional precission.
     * @param info some debug information.
     */
    public static final void search(double ax, double bx, double cx, double[] X,
                                    RealFunctionOfOneVariable f, double tol,final Info info) {

      int iter;

      double d = 0.0, etmp, fu, fv, fw, fx, p, q, r;
      double tol1, tol2, u, xm;
      double e = 0.0;

      double a = (ax < cx ? ax : cx);
      double b = (ax > cx ? ax : cx);

      double x, w, v;

      x = w = v = bx;

      fw = fv = fx = f.eval(x);

      double inValue = fw;

      if(info!=null) info.setMaxIter(ITMAX);

      for (iter = 1; iter <= ITMAX; iter++) {
        xm = 0.5 * (a + b);
        tol2 = 2.0 * (tol1 = tol * Math.abs(x) + ZEPS);
        if (Math.abs(x - xm) <= (tol2 - 0.5 * (b - a))) {
          X[0] = x;
          X[1] = fx;
          if(info!=null){
            info.setCurrentIter(iter);
          }
          return;
        }
        if (Math.abs(e) > tol1) {
          r = (x - w) * (fx - fv);
          q = (x - v) * (fx - fw);
          p = (x - v) * q - (x - w) * r;
          q = 2.0 * (q - r);
          if (q > 0.0)
            p = -p;
          q = Math.abs(q);
          etmp = e;
          e = d;
          if (Math.abs(p) >= Math.abs(0.5 * q * etmp) ||
              p < q * (a - x) || p >= q * (b - x))
            d = CGOLD * (e = (x >= xm ? a - x : b - x));
          else {
            d = p / q;
            u = x + d;
            if (u - a < tol2 || b - u < tol2)
              d = sign(tol1, xm - x);
          }
        }
        else {
          d = CGOLD * (e = (x >= xm ? a - x : b - x));
        }

        u = (Math.abs(d) >= tol1 ? x + d : x + sign(tol1, d));

        fu = f.eval(u);

        if (fu <= fx) {
          if (u >= x)
            a = x;
          else
            b = x;
          v = w;
          w = x;
          x = u; //SHFT(v, w, x, u);
          fv = fw;
          fw = fx;
          fx = fu; //SHFT(fv, fw, fx, fu);
        }
        else {
          if (u < x)
            a = u;
          else
            b = u;
          if (fu <= fw || w == x) {
            v = w;
            w = u;
            fv = fw;
            fw = fu;
          }
          else if (fu <= fv || v == x || v == w) {
            v = u;
            fv = fu;
          }
        }
      }

      X[0] = x;
      X[1] = fx;

   // save some debug information
      if( info!=null){
        String str = "Too many iteration in BRENT\n";
        if(fx > inValue)
          str+=" proc Brent failed to decrease center value! " + ax +
                           " " + bx + " " + cx;
        info.setMessage(str);
        info.printDebug();
      }
    }


}
