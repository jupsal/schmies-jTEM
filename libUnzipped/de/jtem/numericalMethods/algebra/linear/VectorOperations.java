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

package de.jtem.numericalMethods.algebra.linear;

/** This class handles the basic operations with vectors.
 ** A vector is a one dimensional array (double, int ).
 ** All methods assume that the passed data have sensible sizes,
 ** which, due to performance reasons, is never checked. So be careful!
 **
 ** @author	schmies
 */

public final class VectorOperations {

    private VectorOperations() {}

    /* ************************ */
    /* ******** times ********* */
    /* ************************ */




    /** Performs dot product.
     * @param   a  The first multiplicand.
     * @param   b  The second multiplicand.
     * @return  The product.
     */
    public final static int times( final int [] a, final int [] b ) {

        final int n = a.length;

	int sum = 0;

	for(int k=0; k<n; k++) {
	    sum += a[k]*b[k];
	}

	return sum;
    }

     /** Performs dot product
     * @param   VRe  The real part of first multiplicand.
     * @param   VIm  The imag. part of first multiplicand.
     * @param   WRe  The real part second multiplicand.
     * @param   WIm  The imag. part second multiplicand.
     * @return  The product.
     */
    public final static double dot( final double [] VRe, final double [] VIm,
			     final double [] WRe, final double [] WIm)
    {
	final int length = VRe.length;
	double rr=0.0;
	for(int i=0; i<length; i++) rr+=WRe[i]*VRe[i]+WIm[i]*VIm[i];

	return rr;
    }

    /** Performs dot product.
     * @param   a  The first multiplicand.
     * @param   b  The second multiplicand.
     * @return  The product.
     */
    public final static double times( final double [] a, final double [] b ) {

        final int n       = a.length;

	double sum = 0;

	for(int k=0; k<n; k++) {
	    sum += a[k]*b[k];
	}

	return sum;
    }


    /** Performs dot product.
     * @param   a  The first multiplicand.
     * @param   b  The second multiplicand.
     * @return  The product.
     */
    public final static double times( final int [] a, final double [] b ) {

        final int n       = a.length;

	double sum = 0.0;

	for(int k=0; k<n; k++) {
	    sum += a[k]*b[k];
	}

	return sum;
    }


    /** Performs dot product.
     * @param   a  The first multiplicand.
     * @param   b  The second multiplicand.
     * @return  The product.
     */
    public final static double times( final double [] a, final int [] b ) {

        final int n       = a.length;

	double sum = 0;

	for(int k=0; k<n; k++) {
	    sum += a[k]*b[k];
	}

	return sum;
    }


    /** Scales vector.
     * @param a Vector to scale.
     * @param b Scalar.
     * @param c Product.
     */
    public final static void times( final int [] a, final int b, final int [] c ) {

        final int n = a.length;

	for(int k=0; k<n; k++) {
	    c[k] = a[k]*b;
	}
    }

    /** Scales vector.
     * @param a Vector to scale.
     * @param b Scalar.
     * @param c Product.
      */
    public final static void times( final double [] a, final double b, final double [] c ) {

        final int n = a.length;

	for(int k=0; k<n; k++) {
	    c[k] = a[k]*b;
	}
    }


    /** Scales vector.
     * @param a Vector to scale.
     * @param b Scalar.
     * @param c Product.
      */
    public final static void times( final int [] a, final double b, final double [] c ) {

        final int n = a.length;

	for(int k=0; k<n; k++) {
	    c[k] = a[k]*b;
	}
    }
    /** Scales vector.
     * @param VRe  The real part of Vector to scale.
     * @param VIm  The imag. part of Vector to scale.
     * @param bRe  The real part of Scalar.
     * @param bIm  The imag part of Scalar.
     * @param WRe  The real part of Product.
     * @param WIm  The imag. part of Product.
     */
    public final static void times( final double [] VRe, final double [] VIm,
			      final double bRe, final double bIm,
			      final double [] WRe, final double [] WIm)
    {
	final int length = VRe.length ;

	for(int i=0; i<length; i++) {
		final double rr = VRe[i];
		final double ii = VIm[i];

		WRe[i]=rr*bRe-ii*bIm;
		WIm[i]=rr*bIm+ii*bRe;
	    }
    }

    /** Scales vector.
     * @param VRe  The real part of Vector to scale.
     * @param VIm  The imag. part of Vector to scale.
     * @param bRe  The real part of Scalar.
     * @param bIm  The imag part of Scalar.
     * @param WRe  The real part of Product.
     * @param WRe  The imag. part of Product.
     */
    public final static void divide( final double [] VRe, final double [] VIm,
			      final double bRe, final double bIm,
			      final double [] WRe, final double [] WIm)
    {
	final int length = VRe.length;
    	final double nn = bRe*bRe + bIm*bIm;

	for(int i=0; i<length; i++) {
	    final double rr=VRe[i];
	    final double ii=VIm[i];

	    WRe[i]=(rr*bRe+ii*bIm)/nn;
	    WIm[i]=(ii*bRe-rr*bIm)/nn;
	}
    }


    /** Scales vector.
     * @param a  The Vector to scale.
     * @param b  The Scalar.
     * @param c  The Product.
     */
    public final static void divide( final int [] a, final int b, final int [] c ) {

        final int n = a.length;

	for(int k=0; k<n; k++) {
	    c[k] = a[k]/b;
	}
    }
    /** Scales vector: a / b = c .
     * @param b  The Vector to scale.
     * @param a  The Scalar.
     * @param c  The Product.
     */
    public final static void divide( final int a, final int []  b, final int [] c ) {

        final int n = b.length;

	for(int k=0; k<n; k++) {
	    c[k] = a/b[k];
	}
    }

   /** Scales vector: a / b = c .
     * @param b  The Vector to scale.
     * @param a  The Scalar.
     * @param c  The Product.
     */
    public final static void divide( final double a, final int []  b, final double [] c ) {

        final int n = b.length;

	for(int k=0; k<n; k++) {
	    c[k] = a/b[k];
	}
    }

    /** Scales vector: a / b = c .
     * @param b  The Vector to scale.
     * @param a  The Scalar.
     * @param c  The Product.
     */
    public final static void divide( final double a, final double [] b, final double [] c ) {

        final int n = b.length;

	for(int k=0; k<n; k++) {
	    c[k] = a/b[k];
	}
    }

    /** Scales vector: a / (re+i*im) = (resRe+i*resIm) .
     * @param b  The Vector to scale.
     */
    public final static void divide( final double a, final double [] re, final double [] im, 
				     final double [] resRe, final double [] resIm ) {

        final int n = re.length;

	for(int k=0; k<n; k++) {
	    final double Re = re[k];
	    final double Im = im[k];
	    final double nn = Re*Re+Im*Im;
	    if( nn == 0. )
		{
		    resRe[k] = a /nn;
		    resIm[k] = 0.;
		}
	    else
		{
		    resRe[k] =  a * Re /nn;
		    resIm[k] = -a * Im /nn;
		}
	}
    }

    /** Scales vector a / b = c.
     * @param a  The Vector to scale.
     * @param b  The Scalar.
     * @param c  The Product.
     */
    public final static void divide( final double [] a, final double b, final double [] c ) {

        final int n = a.length;

	for(int k=0; k<n; k++) {
	    c[k] = a[k]/b;
	}
    }

    /** scales vector a /b = c.
     * @param a  The Vector to scale.
     * @param b  The Scalar.
     * @param c  The Product.
     */
    public final static void divide( final int [] a, final double b, final double [] c ) {

        final int n = a.length;

	for(int k=0; k<n; k++) {
	    c[k] = a[k]/b;
	}
    }

    /* ************************ */
    /* ******** plus ********** */
    /* ************************ */


    /** Performs a + b = c.
     * @param a  The first summend.
     * @param b  The second summend.
     * @param c  The Product.
     */
    public final static void plus( final int [] a, final int [] b, final int [] c ) {

        final int n = a.length;

	for(int k=0; k<n; k++) {
	    c[k] = a[k] + b[k];
	}
    }

    /** Performs a + b = c.
     * @param a  The first summend.
     * @param b  The second summend.
     * @param c  The Product.
     */
    public final static void plus( final double [] a, final double [] b, final double [] c ) {

        final int n = a.length;

	for(int k=0; k<n; k++) {
	    c[k] = a[k] + b[k];
	}
    }


    /** Performs a + b = c.
     * @param a  The first summend.
     * @param b  The second summend.
     * @param c  The Product.
     */
    public final static void plus( final int [] a, final double [] b, final double [] c ) {

        final int n = a.length;

	for(int k=0; k<n; k++) {
	    c[k] = a[k] + b[k];
	}
    }


    /** Performs a + b = c.
     * @param a  The first summend.
     * @param b  The second summend.
     * @param c  The Product.
     */
    public final static void plus( final double [] a, final int [] b, final double [] c ) {

        final int n = a.length;

	for(int k=0; k<n; k++) {
	    c[k] = a[k] + b[k];
	}
    }

    /** Performs a + b = c.
     * @param a  The first summend.
     * @param b  The second summend.
     * @param c  The Product.
     */
    public final static void plus( final int [] a, final int b, final int [] c ) {

	final int n = a.length;

	for(int k=0; k<n; k++)
	    c[k]=a[k]+b;
    }
    /** Performs a + b = c .
     * @param a  The first summend.
     * @param b  The second summend.
     * @param c  The Product.
     */
    public final static void plus( final double [] a, final double b, final double [] c ) {

	final int n = a.length;

	for(int k=0; k<n; k++)
	    c[k]=a[k]+b;
    }

    /* ************************ */
    /* ******** minus ********* */
    /* ************************ */


    /** Performs a - b = c .
     * @param  b   The subtrahend.
     * @param  c   The product.
     */
    public final static void minus( final int [] a, final int [] b, final int [] c ) {

        final int n = a.length;

	for(int k=0; k<n; k++) {
	    c[k] = a[k] - b[k];
	}
    }

    /** Performs a - b = c.
     * @param  b   The subtrahend.
     * @param  c   The difference.
     */
    public final static void minus( final double [] a, final double [] b, final double [] c ) {

        final int n = a.length;

	for(int k=0; k<n; k++) {
	    c[k] = a[k] - b[k];
	}
    }


    /** Performs a - b = c.
     * @param  b   The subtrahend.
     * @param  c   The product.
     */
   public final static void minus( final int [] a, final double [] b, final double [] c ) {

        final int n = a.length;

	for(int k=0; k<n; k++) {
	    c[k] = a[k] - b[k];
	}
    }


    /** Performs a - b = c .
     * @param  b   The subtrahend.
     * @param  c   The product.
     */
    public final static void minus( final double [] a, final int [] b, final double [] c ) {

        final int n = a.length;

	for(int k=0; k<n; k++) {
	    c[k] = a[k] - b[k];
	}
    }

    /** Performs a - b = c .
     * @param  b   The subtrahend.
     * @param  c   The product.
     */
    public final static void minus( final int [] a, final int b, final int [] c ) {

	final int n = a.length;

	for(int k=0; k<n; k++)
	    c[k]=a[k]-b;
    }
    /** Performs a - b = c .
     * @param  b   The subtrahend.
     * @param  c   The product.
     */
    public final static void minus( final int a, final int [] b, final int [] c ) {

	final int n = b.length;

	for(int k=0; k<n; k++)
	    c[k]=a-b[k];
    }
    /** Performs a - b = c .
     * @param  b   The subtrahend.
     * @param  c   The product.
     */
    public final static void minus( final double [] a, final double b, final double [] c ) {

	final int n = a.length;

	for(int k=0; k<n; k++)
	    c[k]=a[k]-b;
    }
    /** Performs a - b = c.
     * @param  b   The subtrahend.
     * @param  c   The product.
     */
    public final static void minus( final double a, final double [] b, final double [] c ) {

	final int n = b.length;

	for(int k=0; k<n; k++)
	    c[k]=a-b[k];
    }
    /** Assign vector <code>b</code> to vector <code>a</code>.**/
    final static void assign( final int [] a, final int [] b ) {
	System.arraycopy( a, 0, b, 0, a.length );
    }


    /** Assign vector <code>b</code> with vector <code>a</code>. **/
    public final static void assign( final double [] a, final double [] b ) {
	System.arraycopy( a, 0, b, 0, a.length );
    }

    /** Assign vector <code>b</code> with vector <code>a</code>. **/
    public final static void assign( final int [] a, final double [] b ) {
	final int n = a.length;
	for(int j=0; j<n; j++)
	    b[j] = a[j];
    }

    /** Assign vector <code>a</code> to scalar <code>v</code>. **/
    public final static void assign( final int [] a, final int v ) {
	final int n = a.length;
	for(int j=0; j<n; j++)
	    a[j] = v;
    }

    /** Assign vector <code>a</code> to scalar <code>v</code>. **/
    public final static void assign( final double [] a, final double v ) {
	final int n = a.length;
	for(int j=0; j<n; j++)
	    a[j] = v;
    }

    /** Assign vector <code>a</code> to zero. **/
    public final static void assignZero( final double [] a ){
	final int n = a.length;
	for(int j=0; j<n; j++)
	    a[j] = 0.;
    }

    /** Assign vector <code>a</code> to zero. **/
    public final static void assignZero( final int [] a ){
	final int n = a.length;
	for(int j=0; j<n; j++)
	    a[j] = 0;
    }

    /** Returns a copy of vector <code>a</code>.
     * @param a   The vector to copy.
     * @return  The copy-vector.
     */
    final static int [] copy( final int [] a ) {
	return (int[])a.clone();
    }

    /** Returns a clone of vector <code>a</code>.
     * @param a   The vector to copy.
     * @return  The copy-vector.
     */
    public final static double [] copy( final double [] a ) {
	return (double[])a.clone();
    }

    /** Rounds vector <code>a</code> and stores result in <code>b</code> .
     ** a and b may coinside.
     * @param  a  The vector to round.
     * @param  b  The rounded vector.
     */
    public final static void round( final double [] a, final double [] b ) {

	final int n = a.length;

	for(int j=0; j<n; j++)
	    b[j] = Math.floor( a[j] + 0.5 );
    }

    /** Rounds vector <code>a</code> and stores result in <code>b</code> .
     ** a and b may coinside.
     * @param  a  The vector to round.
     * @param  b  The rounded vector.
     */
    public final static void floor( final double [] a, final double [] b ) {

	final int n = a.length;

	for(int j=0; j<n; j++)
	    b[j] = Math.floor( a[j] );
    }

    /** Rounds vector <code>a</code> and stores result in <code>b</code> .
     ** a and b may coinside.
     * @param  a  The vector to round.
     * @param  b  The rounded vector.
     */
    public final static void round( final double [] a, final int [] b ) {

	final int n = a.length;

	for(int j=0; j<n; j++)
	    b[j] = (int)Math.floor( a[j] + 0.5 );
    }

    /** Rounds vector <code>a</code> and stores result in <code>b</code> .
     ** a and b may coinside.
     * @param  a  The vector to round.
     * @param  b  The rounded vector.
     */
    public final static void floor( final double [] a, final int [] b ) {

	final int n = a.length;

	for(int j=0; j<n; j++)
	    b[j] = (int)Math.floor( a[j] );
    }

    /** Negates vector <code>a</code> and stores result in <code>b</code>.
     ** a and b may coinside.
     * @param  a   The vector to negate.
     * @param  b   The negated vector.
     */
    public final static void neg( final int [] a, final int [] b ) {

	final int n = a.length;

	for(int j=0; j<n; j++)
	    b[j] = -a[j];
    }

    /** Negates vector <code>a</code> and stores result in <code>b</code>.
     ** a and b may coinside.
     * @param  a   The vector to negate.
     * @param  b   The negated vector.
     */
    public final static void neg( final double [] a, final double [] b ) {

	final int n = a.length;

	for(int j=0; j<n; j++)
	    b[j] = -a[j];
    }

    /** Negates vector <code>a</code> and stores result in <code>b</code>.
     ** a and b may coinside.
     * @param  a   The vector to negate.
     * @param  b   The negated vector.
     */
    public final static void neg( final int [] a, final double [] b ) {

	final int n = a.length;

	for(int j=0; j<n; j++)
	    b[j] = -a[j];
    }

    /* Return square of norm of <i>re</i>. */
    public static double normSqr( final double [] re )
    {
	int l = re.length;
	double n = 0.;
	for (int i = 0; i < l; i++) n += re[i] * re[i];
	return n;
    }

    /* Return square of norm of <i>re</i>. */
    public static double normSqr( final int [] re )
    {
	int l = re.length;
	double n = 0.;
	for (int i = 0; i < l; i++) n += re[i] * re[i];
	return n;
    }

    /** Throw <b>IllegalArgumentException</b> if <i>v</i>' length does not equal <i>w</i>' length.  */
    public final static  void checkShape( final double [] v, final double [] w )
    {
	if(w.length!=v.length)
	    throw new IllegalArgumentException("different vec sizes");
    }

    /** Throw <b>IllegalArgumentException</b> if <i>v</i>' length does not equal <i>w</i>' length.  */
    public final static  void checkShape( final int [] v, final int [] w )
	{
	    if(w.length!=v.length)
		throw new IllegalArgumentException("different vec sizes");
	}

    /** Throw <b>IllegalArgumentException</b> if <i>v</i>' length does not equal <i>w</i>' length.  */
    public final static  void checkShape( final double [] v, final int [] w )
	{
	    if(w.length!=v.length)
		throw new IllegalArgumentException("different vec sizes");
	}

    /** Throw <b>IllegalArgumentException</b> if <i>v</i>' length does not equal <i>w</i>' length.  */
    public final static  void checkShape( final int [] v, final double [] w )
	{
	    if(w.length!=v.length)
		throw new IllegalArgumentException("different vec sizes");
	}

    /** Assigns random all entries of v.*/
    public static final void random( double [] v) {
	final int l = v.length;

        for (int i = 0; i < l; i++)
            v[i] = 2. * Math.random() - 1.;
    }

     /** Assigns random all entries of v.
      * all entries are from 0 to range.
      */
    public static final void random( int [] v, int range) {
	final int l = v.length;

        for (int i = 0; i < l; i++)
            v[i] = (int) Math.random()*range;
    }

}//end of class





