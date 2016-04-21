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

public final class DBrent
    implements java.io.Serializable {

  private static final long serialVersionUID = 1L;

  static final double CGOLD = 0.3819660;
  static final double ZEPS = 1.0e-17;

  static int ITMAX = 100;

  /**
   * Get the value of ITMAX.
   * @return Value of ITMAX.
   */
  public static int getITMAX() {
    return ITMAX;
  }

  /**
   * Set the value of ITMAX.
   * @param v  Value to assign to ITMAX.
   */
  public static void setITMAX(int v) {
    DBrent.ITMAX = v;
  }

  private static final double sign(double a, double b) {
    return b > 0 ? Math.abs(a) : -Math.abs(a);
  }

  /**
   * Search the minimum of function f in the interval (abc[0],abc[2]) with her first derivative
   *  and write it in the X (such X must be a double array of dim. 2).
   * @param abc bracketing triplet of function f with abc[0]<abc[1]<t[2] or abc[0]>abc[1]>abc[2]
   * and f(abc[0]),f(abc[2])>f(abc[1]).
   * @param X the minimum point of function f in the interval (abc[0],abc[2]).X[0] is the position
   *  and X[1] is the value of minimum.
   * @param f given function.
   * @param df first derivative of f.
   * @param tol the fractional precission.
   */
  public static final void search(double[] abc,
                                  double[] X,
                                  RealFunctionOfOneVariable f,
                                  RealFunctionOfOneVariable df,
                                  double tol) {
    search(abc[0], abc[1], abc[2], X, f, df, tol,null);
  }

  /**
  * Search the minimum of function f in the interval ( ax,cx ) with her first derivative
  *  and write it in the X (such X must be a double array of dim. 2).
  * @param ax first component of bracketing triplet.
  * @param bx second component of bracketing triplet.
  * @param cx third component of bracketing triplet.
  * @param X the minimum point of function f in the interval (abc[0],abc[2]).X[0] is the position
  *  and X[1] is the value of minimum.
  * @param f given function.
  * @param df first derivative of f.
  * @param tol the fractional precission.
  */
 public static final void search(double ax, double bx, double cx,
                                 double[] X,
                                 RealFunctionOfOneVariable f,
                                 RealFunctionOfOneVariable df,
                                 double tol) {
   search(ax,bx,cx,X,f,df,tol,null);
 }

 /**
   * Search the minimum of function f in the interval ( ax,cx ) with her first derivative
   *  and write it in the X (such X must be a double array of dim. 2).
   * @param ax first component of bracketing triplet.
   * @param bx second component of bracketing triplet.
   * @param cx third component of bracketing triplet.
   * @param X the minimum point of function f in the interval (abc[0],abc[2]).X[0] is the position
   *  and X[1] is the value of minimum.
   * @param f given function.
   * @param df first derivative of f.
   * @param tol the fractional precission.
   * @param info some debug information.
   */
  public static final void search(double ax, double bx, double cx,
                                  double[] X,
                                  RealFunctionOfOneVariable f,
                                  RealFunctionOfOneVariable df,
                                  double tol,final Info info) {

    boolean ok1, ok2;

    double a, b, d = 0.0, d1, d2, du, dv, dw, dx, e = 0.0;
    double fu, fv, fw, fx, olde, tol1, tol2, u, u1, u2, v, w, x, xm;

    a = (ax < cx ? ax : cx);
    b = (ax > cx ? ax : cx);
    x = w = v = bx;

    fw = fv = fx = f.eval(x);
    dw = dv = dx = df.eval(x);

    if(info!=null) info.setMaxIter(ITMAX);

    for (int iter = 1; iter <= ITMAX; iter++) {
      xm = 0.5 * (a + b);
      tol1 = tol * Math.abs(x) + ZEPS;
      tol2 = 2.0 * tol1;
      if (Math.abs(x - xm) <= (tol2 - 0.5 * (b - a))) {
        X[0] = x;
        X[1] = fx;
        if(info!=null) info.setCurrentIter(iter);
        return;
      }
      if (Math.abs(e) > tol1) {
        d1 = 2.0 * (b - a);
        d2 = d1;
        if (dw != dx)
          d1 = (w - x) * dx / (dx - dw);
        if (dv != dx)
          d2 = (v - x) * dx / (dx - dv);
        u1 = x + d1;
        u2 = x + d2;
        ok1 = (a - u1) * (u1 - b) > 0.0 && dx * d1 <= 0.0;
        ok2 = (a - u2) * (u2 - b) > 0.0 && dx * d2 <= 0.0;
        olde = e;
        e = d;
        if (ok1 || ok2) {
          if (ok1 && ok2)
            d = (Math.abs(d1) < Math.abs(d2) ? d1 : d2);
          else if (ok1)
            d = d1;
          else
            d = d2;
          if (Math.abs(d) <= Math.abs(0.5 * olde)) {
            u = x + d;
            if (u - a < tol2 || b - u < tol2)
              d = sign(tol1, xm - x);
          }
          else
            d = 0.5 * (e = (dx >= 0.0 ? a - x : b - x));
        }
        else
          d = 0.5 * (e = (dx >= 0.0 ? a - x : b - x));
      }
      else
        d = 0.5 * (e = (dx >= 0.0 ? a - x : b - x));

      if (Math.abs(d) >= tol1) {
        u = x + d;
        fu = f.eval(u);
      }
      else {
        u = x + sign(tol1, d);
        fu = f.eval(u);

        if (fu > fx) {
          X[0] = x;
          X[1] = fx;
          if(info!=null) info.setCurrentIter(iter);
          return;
        }
      }

      du = df.eval(u);
      if (fu <= fx) {
        if (u >= x)
          a = x;
        else
          b = x;
        v = w;
        fv = fw;
        dv = dw; //MOV3(v, fv, dv, w, fw, dw);
        w = x;
        fw = fx;
        dw = dx; //MOV3(w, fw, dw, x, fx, dx);
        x = u;
        fx = fu;
        dx = du; //MOV3(x, fx, dx, u, fu, du);
      }
      else {
        if (u < x)
          a = u;
        else
          b = u;
        if (fu <= fw || w == x) {
          v = w;
          fv = fw;
          dv = dw; //MOV3(v, fv, dv, w, fw, dw);
          w = u;
          fw = fu;
          dw = du; //MOV3(w, fw, dw, u, fu, du);
        }
        else if (fu < fv || v == x || v == w) {
          v = u;
          fv = fu;
          dv = du; //MOV3(v, fv, dv, u, fu, du);
        }
      }
    }

    if( info!=null){
      info.setCurrentIter(ITMAX);
      info.setMessage(info.getMessage()+"\n"+"Too many iterations in routine DBRENT\n");
      info.printDebug();
    }

    X[0] = x;
    X[1] = fx;
  }

}

