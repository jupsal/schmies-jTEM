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

package de.jtem.riemann.schottky;

import de.jtem.blas.ComplexMatrix;
import de.jtem.blas.ComplexVector;
import de.jtem.mfc.field.Complex;
import junit.framework.TestCase;

public class TestSchottky extends TestCase {

  static double[][] schottkyDataForHe = new double [][] {
      { // He1
         1.47361, 0,
        -1.47361, 0, -0.0265972, 0 },

      { // He2
         1.899744838484059,  1.3383276978981307,
        -1.899744838484059, -1.3383276978981307, 8.776111981589E-5, -1.0564238051974E-4,

         1.899744838484059, -1.3383276978981307,
        -1.899744838484059,  1.3383276978981307, 8.776111981589E-5,  1.0564238051974E-4 },

      { // He3
          1.349015,  0,
         -1.349015,  0, -2.82698e-2, 0,
          3.056639,  1.378618,
         -3.056639, -1.378618, 3.75145e-6, -2.581535e-6,
          3.056639, -1.378618,
         -3.056639,  1.378618, 3.75145e-6,  2.581535e-6 },

      { // He4
         1.7717633047238306,  1.2907336759162422,
        -1.7717633047238306, -1.2907336759162422, 1.2228783186228273E-4, -1.3731780537305815E-4,

         1.7717633047238306, -1.2907336759162422,
        -1.7717633047238306,  1.2907336759162422, 1.2228783186228273E-4,  1.3731780537305815E-4,

         4.5750218512914410,  1.5333453003187825,
        -4.5750218512914410, -1.5333453003187825, 3.7557520421193814E-8, -2.6313463901539743E-8,

         4.5750218512914410, -1.5333453003187825,
        -4.5750218512914410,  1.5333453003187825, 3.7557520421193814E-8,  2.6313463901539743E-8 }
      };

      /** Period matrices for helicoid computed via schottky uniformization
       * with precesion 10^-13. A comparism with results computed via
       * ramified covering shows that these values are atleas 9 digits precise. */
    static ComplexMatrix [] periodMatrixForHe = new ComplexMatrix[] {

        new ComplexMatrix( new double [][]
            {{ -3.6269493318929333 }},
                           new double [][]
            {{  3.141592653589793  }} ),

        new ComplexMatrix( new double [][]
            {{-8.893841400451032, -0.7005977257529801},
             {-0.7005977257529801, -8.893841400451029}},
                           new double [][]
            {{-0.8785450734734299, 3.1415926535897927},
             {3.1415926535897927, 0.8785450734734307}} ),
        new ComplexMatrix( new double [][]
            {{-3.5659698343995747, -1.4871742106353656, -1.487174210635264},
             {-1.4871742106353656, -12.37620011550392, -1.7534784593908188},
             {-1.487174210635264, -1.7534784593908188, -12.376200115503982}},
                           new double [][]
            {{3.14159265358979, 0.7518734356222462, -0.7518734356221779},
             {0.7518734356222462, -0.45967872251396846, 3.141592653589757},
             {-0.7518734356221779, 3.141592653589757, 0.45967872251397424}} ),
        new ComplexMatrix( new double [][]
            {{-8.602326913647333, -0.6335291647252136, -1.8078866234012052, -0.9345739224309066},
             {-0.6335291647252136, -8.602326913647364, -0.9345739224308964, -1.8078866234012057},
             {-1.8078866234012052, -0.9345739224308964, -16.895842115969558, -2.184347339559907},
             {-0.9345739224309066, -1.8078866234012057, -2.184347339559907, -16.895842115969558}},
                           new double [][]
            {{-0.8444415657810661, -3.1415926535897714, -0.6658900839052415, -1.503325882255731},
             {-3.1415926535897714, 0.844441565781059, 1.5033258822557398, 0.6658900839052563},
             {-0.6658900839052415, 1.5033258822557398, -0.6120701299917983, -3.1415926535897856},
             {-1.503325882255731, 0.6658900839052563, -3.1415926535897856, 0.6120701299918495}} ) };


  static double [] epsForHe = new double[] { 0, 1e-6, 1e-9, 1e-5, 1e-2 };

  public static SchottkyEstimates getSchottkyOfHe( int g ) {
    return  new SchottkyEstimates
        ( schottkyDataForHe[g-1], epsForHe[g] );
  }


  public static double distSqrMod2PiI(Complex u, Complex v) {

    final double diff = u.getIm() - v.getIm();

    final int N = (int) Math.floor( (u.getIm() - v.getIm()) / 2 / Math.PI + 0.5);

    final double re = u.getRe() - v.getRe();
    final double im = u.getIm() - v.getIm() - N * 2 * Math.PI;

    return re * re + im * im;
  }

  public static double distSqrMod2PiI(ComplexMatrix u, ComplexMatrix v) {
    double distSqr = 0;

    for (int i = 0; i < u.getNumRows(); i++)
      for (int j = 0; j < u.getNumCols(); j++)
        distSqr += distSqrMod2PiI(u.get(i, j), v.get(i, j));

    return distSqr;
  }

  public static double distMaxMod2PiI(ComplexMatrix u, ComplexMatrix v) {
    double max = 0;

    for (int i = 0; i < u.getNumRows(); i++)
      for (int j = 0; j < u.getNumCols(); j++) {
        final double distSqr = distSqrMod2PiI(u.get(i, j), v.get(i, j));
        if( distSqr > max )
          max = distSqr;
      }

    return Math.sqrt( max );
  }

  public static boolean equalsMod2PiI(Complex u, Complex v, double eps) {

    if (Math.abs(u.getRe() - v.getRe()) > eps)
      return false;

    final double diff = u.getIm() - v.getIm();

    final int N = (int) Math.floor( (u.getIm() - v.getIm()) / 2 / Math.PI + 0.5);

    return Math.abs(u.getIm() - v.getIm() - N * 2 * Math.PI) <= eps;
  }

  public static boolean equalsMod2PiI(ComplexMatrix u, ComplexMatrix v,
                                      double eps) {
    if (v.getNumCols() != u.getNumCols() || v.getNumRows() != u.getNumRows())
      return false;

    for (int i = 0; i < u.getNumRows(); i++)
      for (int j = 0; j < u.getNumCols(); j++)
        if (!equalsMod2PiI(u.get(i, j), v.get(i, j), eps))
          return false;

    return true;
  }

  public void testHelicoid( int g ) {

    final int i = g-1;

    SchottkyEstimates schottky = getSchottkyOfHe( g );

    // check period matrix

    assertEquals( "period Matrix acc=1e-4", 0.0,
                  distMaxMod2PiI( schottky.getPeriodMatrix( 1e-4 ) ,
                                  schottky.getPeriodMatrix( 1e-10 ) ),1e-4 );
    assertEquals( "period Matrix acc=1e-6", 0.0,
                  distMaxMod2PiI( schottky.getPeriodMatrix( 1e-6 ) ,
                                  schottky.getPeriodMatrix( 1e-10 ) ),1e-6 );
    assertEquals( "period Matrix acc=1e-8", 0.0,
                  distMaxMod2PiI( schottky.getPeriodMatrix( 1e-8 ) ,
                                  schottky.getPeriodMatrix( 1e-10 ) ),1e-8 );



    assertEquals( "period matrix", 0.0,
                  distMaxMod2PiI( schottky.getPeriodMatrix(),
                                  periodMatrixForHe[i] ), epsForHe[g] );

    Complex test    = new Complex();
    Complex control = new Complex();


    for( int n=0; n<g; n++ ) {

      assertEquals( "V", 0.0, Math.sqrt(distSqrMod2PiI( schottky.V( n, 1e-6 ), schottky.V( n, 1e-8 ) ) ),1e-6 );

    }

    schottky.testVEstimates();
    schottky.testPeriodMatrixEstimates();

    // testing numerical integrals of a-circles
    testAPeriods( schottky );

  }
//
//  public void testHe1() {
//    testHelicoid( 1 );
//  }

  public void testHe2() {
    testHelicoid(2);
  }
//
//  public void testHe3() {
//    testHelicoid(3);
//  }

//  public void testHe4() {
//    testHelicoid(4);
//  }

  /** tests whether numerical integration of a-circles gives right
   * normalized results: 0, 2 pi i, - 2 pi i.
   */
  public void testAPeriods( Schottky schottky ) {

    AbelMap abel = new AbelMap(schottky);

    final int n = schottky.getNumGenerators(); // should be g

    final double acc = schottky.getAccuracy();

    for (int i = 0; i < n; i++) {
      for (int j = 0; j < 2; j++) {

        ComplexVector path = AbelMapTest.spiral
            (schottky.center[i][j],
             schottky.radius[i] * 1.2, 0, 2 * Math.PI);

        for (int k = 0; k < n; k++) {
          Complex integral
              = abel.numEval(k, path, acc / 10 );

          if( k != i ) {// might lose one digit because of numerical integration
            assertEquals( "real part of differential: " + k +" circle " + i + "," + j,
                          0, integral.re, acc *10 );
            assertEquals( "imag part of differential: " + k +" circle " + i + "," + j,
                          0, integral.im, acc *10 );
          } else {
            assertEquals( "real part of differential: " + k +" circle " + i + "," + j,
                          0,                               integral.re, acc *10 );
            assertEquals( "imag part of differential: " + k +" circle " + i + "," + j,
                          (j == 0 ? -1 : 1) * 2 * Math.PI, integral.im, acc *10 );
            }
        }
      }
    }

  }


  final static class SchottkyEstimates extends SchottkyAnalysis {

    SchottkyEstimates( Schottky schottky ) {
      super( schottky.getUniformizationData(), schottky.getAccuracy() );
    }

     SchottkyEstimates( double [] uniformizationData, double accuracy ) {
      super( uniformizationData, accuracy );
    }

    double logFactorForPeriodMatrix;
    double factorForLeftIsM;
    double factorForLeftIsNotM;
    double accuracy;

    int n, m;

    final private void B(final double estimateOfParent, final SchottkyGroupElement element) {

      if (element.updateID != updateID)
        updateElement(element);

      double estimate = Double.MAX_VALUE;

      final long noe = getNumOfElementsOfCosetWithWordLength(element.
          wordLength);

            element.diff(fixpoint[n][1], fixpoint[n][0], tmp);

            double dist = Math.abs( tmp.re ) + Math.abs( tmp.im );

      if (element.left != m) {

        final double localLogFactor = ( 1/dist( element, fixpoint[m][0] ) +
                                        1/dist( element, fixpoint[m][1] ) );

        assertTrue( "local log factor: (wl=" +
                    ((int) element.wordLength)+ ")  " +
                    localLogFactor + " > " + logFactorForPeriodMatrix,
                    localLogFactor < logFactorForPeriodMatrix );


        estimate = dist * localLogFactor;

        tmp.assignCrossRatio(fixpoint[m][0], element.imageOfB[n],
                             fixpoint[m][1], element.imageOfA[n] );
        tmp.assignLog();

        assertTrue( "period matrix estimate failed: (wl=" +
                    ((int) element.wordLength)+ ")  " +
                    tmp.abs() + " > " + estimate,
                    tmp.abs() < estimate );

        estimate = dist * logFactorForPeriodMatrix;

        assertTrue( "period matrix estimate failed",
                    estimate < estimateOfParent * theta1 * theta1 );
      }

      if( element.left != m  && factorForLeftIsNotM * dist * noe < accuracy ||
          element.left == m  && factorForLeftIsM    * dist * noe < accuracy )
        return;


      if (element.child == null)
        createLeftChilds(element);

      final SchottkyGroupElement[] child = element.child;

      final int numOfChilds = child.length;

      for (int i = 0; i < numOfChilds; i++)
        B( estimate, child[i]);
    }

      void testPeriodMatrixEstimates() {
        testPeriodMatrixEstimates( getAccuracy() );
      }

      void testPeriodMatrixEstimates( double accuracy ) {

        if( numGenerators == 1 ) { // treat genus one case seperatly
          return;
        }

        this.accuracy = accuracy;

        final double v = theta1 * theta1;

        final double sumForLeftIsM    = (2*numGenerators-2) * r( v );
        final double sumForLeftIsNotM = (2*numGenerators-4) * r( v ) + rPlus( v ) + rMinus( v );

        for( m = 0; m<numGenerators; m++ ) {

          factorForLeftIsM = factorForLeftIsNotM = logFactorForPeriodMatrix
              = 1/periodMatrix.k1( m, fixpoint[m][0] )
              + 1/periodMatrix.k1( m, fixpoint[m][1] );

          factorForLeftIsM    *= sumForLeftIsM;
          factorForLeftIsNotM *= sumForLeftIsNotM;

          for( n=0; n<numGenerators; n++ ) {

            if( m<=n ) {

              for (int i = 0; i < numGenerators; i++)
                if( i != n ) {
                  B( Double.MAX_VALUE, generator[i]);
                  B( Double.MAX_VALUE, generatorInv[i]);
                }
            }
          }
        }
      }


      double factor; // set in V

      private final void V( double distOfParent, SchottkyGroupElement element) {

        if (element.updateID != updateID)
          updateElement(element);

        final long noe = getNumOfElementsOfCosetWithWordLength(element.wordLength);

        element.diff(fixpoint[n][1], fixpoint[n][0], tmp);

        final double dist = tmp.abs();

        final double error = factor * noe * (Math.abs(tmp.re) + Math.abs(tmp.im));

        assertTrue( "V + " + distOfParent + " "+ distOfParent * theta1 * theta1 + "  " + dist, distOfParent * theta1 * theta1 > dist );

        if (error < accuracy)
          return;

        if (element.child == null)
          createLeftChilds(element);

        final SchottkyGroupElement[] child = element.child;

        final int numOfChilds = child.length;

        for (int i = 0; i < numOfChilds; i++)
          V( dist, child[i]);

      }

      final void testVEstimates() {
        for( int i=0; i<numGenerators; i++ )
          testVEstimates( i, getAccuracy() );
      }

      final void testVEstimates( int n, double accuracy) {

        this.n = n;
        this.accuracy = accuracy;
        this.factor = 1 / (theta1 * theta1 * (2 * numGenerators - 1));

        for (int i = 0; i < numGenerators; i++)
          if (i != n) {
            V( Double.MAX_VALUE, generator[i]);
            V( Double.MAX_VALUE, generatorInv[i]);
          }
      }



  }


}
