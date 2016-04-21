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

package de.jtem.numericalMethods.calculus.integration;

/**
 * Collection of Newton-Cotes summation formulas for different orders.
 *
 * @author Markus Schmies
 * @version 1.0
 */
public final class NewtonCotes {

  NewtonCotes() {
  }

  /**
   * Returns sum of values in table
   * @param value table
   * @return sum of values in table
   */
  public static double sum(double[] value) {

    int n = value.length;

    double res = 0;

    for (int i = 0; i < n; i++) {
      res += value[i];

    }
    return res;
  }

  /**
   * Returns average of values in table
   * @param value table
   * @return average of values in table
   */
  public static double average(double[] value) {

    return sum(value) / value.length;
  }

  /**
   * Computes Newton-Cotes formula using Euler's rule.
   * Order 3 method; integrates polynoms up to order 1 exactly.
   * @param value table of legnth 1 + n with n>0, equally spaced interval of length one
   * @return integral of tabled (equally spaced, interval length 1) function
   */
  public static double eulerSum(double[] value) {

    return (2 * sum(value) - (value[0] + value[value.length - 1])) / 2 /
        (value.length - 1);
  }

  /**
   * Computes Newton-Cotes formula using Simpson's rule.
   * Order 5 method; integrates polynoms up to order 3 exactly.
   * @param value table of legnth 1 + 2n with n>0, equally spaced interval of length one
   * @return integral of tabled (equally spaced, interval length 1) function
   */
  public static double simpsonSum(double[] value) {

    if (value.length % 2 != 1) {
      throw new IllegalArgumentException("length of array have to be odd");
    }

    int n = value.length - 1;

    double mod0Sum = value[n];
    double mod1Sum = 0;

    for (int i = 0; i < n; ) {
      mod0Sum += value[i++];
      mod1Sum += value[i++];
    }

    mod0Sum = 2 * mod0Sum - value[0] - value[n];

    return (mod0Sum + 4 * mod1Sum) / 6 / (n / 2);
  }


  /**
   * Computes Newton-Cotes formula using Kepler's rule.
   * Order 5 method; integrates polynoms up to order 3 exactly.
   * @param value table of legnth 1 + 3n with n>0, equally spaced interval of length one
   * @return integral of tabled (equally spaced, interval length 1) function
   */
  public static double keplerSum(double[] value) {

    if (value.length % 3 != 1) {
      throw new IllegalArgumentException("length of array modulo 3 has to be 1");
    }

    int n = value.length - 1;

    double mod0Sum = value[n];
    double mod1Sum = 0;
    double mod2Sum = 0;

    for (int i = 0; i < n; ) {
      mod0Sum += value[i++];
      mod1Sum += value[i++];
      mod2Sum += value[i++];
    }

    mod0Sum = 2 * mod0Sum - value[0] - value[n];

    return (mod0Sum + 3 * mod1Sum + 3 * mod2Sum) / 8 / (n / 3);
  }

  /**
   * Computes Newton-Cotes formula using Milne's rule.
   * Order 7 method; integrates polynoms up to order 5 exactly.
   * @param value table of legnth 1 + 4n with n>0, equally spaced interval of length one
   * @return integral of tabled (equally spaced, interval length 1) function
   */
  public static double milneSum(double[] value) {

    if (value.length % 4 != 1) {
      throw new IllegalArgumentException("length of array modulo 4 has to be 1");
    }

    int n = value.length - 1;

    double mod0Sum = value[n];
    double mod1Sum = 0;
    double mod2Sum = 0;
    double mod3Sum = 0;

    for (int i = 0; i < n; ) {
      mod0Sum += value[i++];
      mod1Sum += value[i++];
      mod2Sum += value[i++];
      mod3Sum += value[i++];
    }

    mod0Sum = 2 * mod0Sum - value[0] - value[n];

    return (7 * mod0Sum + 32 * mod1Sum + 12 * mod2Sum + 32 * mod3Sum) / 90 /
        (n / 4);
  }

  /**
   * Computes Newton-Cotes formula using an order 7 rule different to Milne's.
   * Order 7 method; integrates polynoms up to order 5 exactly.
   * @param value table of legnth 1 + 5n with n>0, equally spaced interval of length one
   * @return integral of tabled (equally spaced, interval length 1) function
   */
  public static double order7Sum(double[] value) {

    if (value.length % 5 != 1) {
      throw new IllegalArgumentException("length of array modulo 5 has to be 1");
    }

    int n = value.length - 1;

    double mod0Sum = value[n];
    double mod1Sum = 0;
    double mod2Sum = 0;
    double mod3Sum = 0;
    double mod4Sum = 0;

    for (int i = 0; i < n; ) {
      mod0Sum += value[i++];
      mod1Sum += value[i++];
      mod2Sum += value[i++];
      mod3Sum += value[i++];
      mod4Sum += value[i++];
    }

    mod0Sum = 2 * mod0Sum - value[0] - value[n];

    return (19 * mod0Sum + 75 * mod1Sum + 50 * mod2Sum + 50 * mod3Sum +
            75 * mod4Sum) / 288 / (n / 5);
  }

  /**
   * Computes Newton-Cotes formula using Weddle's rule.
   * Order 9 method; integrates polynoms up to order 7 exactly.
   * @param value table of legnth 1 + 4n with n>0, equally spaced interval of length one
   * @return integral of tabled (equally spaced, interval length 1) function
   */
  public static double weddleSum(double[] value) {

    if (value.length % 6 != 1) {
      throw new IllegalArgumentException("length of array modulo 6 has to be 1");
    }

    int n = value.length - 1;

    double mod0Sum = value[n];
    double mod1Sum = 0;
    double mod2Sum = 0;
    double mod3Sum = 0;
    double mod4Sum = 0;
    double mod5Sum = 0;

    for (int i = 0; i < n; ) {
      mod0Sum += value[i++];
      mod1Sum += value[i++];
      mod2Sum += value[i++];
      mod3Sum += value[i++];
      mod4Sum += value[i++];
      mod5Sum += value[i++];
    }

    mod0Sum = 2 * mod0Sum - value[0] - value[n];

    return (41 * mod0Sum + 216 * mod1Sum + 27 * mod2Sum + 272 * mod3Sum +
            27 * mod4Sum + 216 * mod5Sum) / 840 / (n / 6);
  }

  /**
    * Computes Newton-Cotes formula using an order 9 rule different to Weddle's.
    * Order 9 method; integrates polynoms up to order 7 exactly.
    * @param value table of legnth 1 + 5n with n>0, equally spaced interval of length one
    * @return integral of tabled (equally spaced, interval length 1) function
    */
  public static double order9Sum(double[] value) {

    if (value.length % 7 != 1) {
      throw new IllegalArgumentException("length of array modulo 7 has to be 1");
    }

    int n = value.length - 1;

    double mod0Sum = value[n];
    double mod1Sum = 0;
    double mod2Sum = 0;
    double mod3Sum = 0;
    double mod4Sum = 0;
    double mod5Sum = 0;
    double mod6Sum = 0;

    for (int i = 0; i < n; ) {
      mod0Sum += value[i++];
      mod1Sum += value[i++];
      mod2Sum += value[i++];
      mod3Sum += value[i++];
      mod4Sum += value[i++];
      mod5Sum += value[i++];
      mod6Sum += value[i++];
    }

    mod0Sum = 2 * mod0Sum - value[0] - value[n];

    return (751 * mod0Sum + 3577 * mod1Sum + 1323 * mod2Sum
            + 2989 * mod3Sum + 2989 * mod4Sum + 1323 * mod5Sum + 3577 * mod6Sum) /
        17280 / (n / 7);
  }

  /**
    * Computes Newton-Cotes formula using the rule with the highest possible order
    * suitable for the length of the provided table.
    * @param value table of legnth (>1), equally spaced interval of length one
    * @return integral of tabled (equally spaced, interval length 1) function
    */
  public static double highestOrderSum(double[] value) {

    if (value.length % 7 == 1) {
      return order9Sum(value);
    }
    else if (value.length % 6 == 1) {
      return weddleSum(value);
    }
    else if (value.length % 5 == 1) {
      return order7Sum(value);
    }
    else if (value.length % 4 == 1) {
      return milneSum(value);
    }
    else if (value.length % 3 == 1) {
      return keplerSum(value);
    }
    else if (value.length % 2 == 1) {
      return simpsonSum(value);
    }

    return eulerSum(value);
  }
}