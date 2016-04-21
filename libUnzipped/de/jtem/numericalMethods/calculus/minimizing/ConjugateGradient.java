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

import de.jtem.numericalMethods.calculus.function.RealFunctionOfSeveralVariablesWithGradient;


/**
 * This class represents the routine to find a minimum of multidimensional function f with conjugate
 * gradient methods in mutidimensional.
 * In short: we start at initial point p then move to the minimum of line from point pi and direction
 * of local gradient as many times as needed.
 * @author Markus Schmies, Vitali Lieder
 * @version 1.0
 */
public final class ConjugateGradient
    implements java.io.Serializable {

  private static final long serialVersionUID = 1L;

  static final double CGOLD = 0.3819660;

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
    ITMAX = v;
  }

  final static double EPS = 1e-10;

  static private boolean useDBrent = false;

  /**
   * Get the value of useDBrent.
   * @return Value of useDBrent.
   */
  static public boolean getUseDBrent() {
    return useDBrent;
  }

  /**
   * Set the value of useDBrent.
   * @param v  Value to assign to useDBrent.
   */
  static public void setUseDBrent(boolean v) {
    useDBrent = v;
  }

  /**
   * Search the minimum of function f  with conjugate gradient methods and return the value of the minimum.
   * @param p starting point.
   * @param ftol precision tolerance.
   * @param f given function.
   * @return the minimum of function f.
   */
  public static double search(double[] p,
                              double ftol,
                              RealFunctionOfSeveralVariablesWithGradient f) {
    return search(p, ftol, f, ITMAX, useDBrent, null);
  }

  /**
   * Search the minimum of function f  with conjugate gradient methods and return the value of the minimum.
   * @param p starting point.
   * @param ftol precision tolerance.
   * @param f given function.
   * @param itMax maximum number of evalutaion for finding a minimum.
   * @param useDBrent wether to use brent's methods with first derivative to find a minimum along line or brent's method.
   * @param info some debugging informations.
   * @return the minimum of function f.
   */
  public static double search(double[] p,
                              double ftol,
                              RealFunctionOfSeveralVariablesWithGradient f,
                              int itMax, boolean useDBrent, final Info info) {

    int n = p.length;

    final double[] h = new double[n];
    final double[] g = new double[n];
    final double[] xi = new double[n];

    final MinimizingOnLine minimizingOnLine;

    if (useDBrent) {
      minimizingOnLine = new DBrentOnLine(p, xi, f);
    }
    else {
      minimizingOnLine = new BrentOnLine(p, xi, f);
    }

    double fp = f.eval(p, xi);
    double fret = fp;

    for (int j = 0; j < n; j++) {
      g[j] = -xi[j];
      xi[j] = h[j] = g[j];
    }

    if(info != null) info.setMaxIter(itMax);

    for (int its = 0; its < itMax; its++) {

      fret = minimizingOnLine.search(2.0e-8);

      if (2.0 * Math.abs(fret - fp) <=
          ftol * (Math.abs(fret) + Math.abs(fp) + EPS)){
        if(info!=null) {
          info.setMessage("ConjugateGradient " + its + ": " + 2.0 * Math.abs(fret - fp) +
                          " (" + ftol * (Math.abs(fret) + Math.abs(fp) + EPS) +
                   ")");
          info.setCurrentIter(its);
          info.printDebug();
        }
        return fret;
      }

      fp = f.eval(p, xi);

      double dgg = 0, gg = 0;

      for (int j = 0; j < n; j++) {
        gg += g[j] * g[j];
        /* dgg += xi[j]*xi[j]; */
        dgg += (xi[j] + g[j]) * xi[j];
      }
      if (gg == 0.0){
        if(info!=null) {
          info.setCurrentIter(its);
          info.printDebug();
        }
        return fret;
      }

      double gam = dgg / gg;
      for (int j = 0; j < n; j++) {
        g[j] = -xi[j];
        xi[j] = h[j] = g[j] + gam * h[j];
      }
    }

    if( info!=null){
      info.setCurrentIter(itMax);
      info.setMessage(info.getMessage()+"\n"+"Too many iterations in ConjugateGradient\n");
    }

    return fret;
  }

}
