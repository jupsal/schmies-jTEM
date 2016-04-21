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

/** This class handles the basic operations with matrices .
 ** A matrix is a two dimensional array. Rows relate to the first index,
 ** columns to the second. Thus a matrix is an array of row vectors which
 ** have all the same length.
 ** All methods assume that the passed data have sensible sizes, 
 ** which, due to performance reasons, is never checked. So be careful!
 **   
 ** @author	schmies
 */

public final class MatrixOperations {

    private MatrixOperations() {}


    /* ************************ */
    /* ******** times ********* */
    /* ************************ */



    /**
     * Performs matrix multiplication: A * B = C.
     * If ARe and/or AIm and/or BRe and/or BIm equals CRe or CIm a temporary
     * copy is assigned, which should be avoided if run times are critical. 
     * <p>
     * @param   ARe  Real part of first multiplicand.
     * @param   AIm  Imag. part of first multiplicand.
     * @param   BRe  Real part of second multiplicand.
     * @param   BIm  Imag. part of second multiplicand.
     * @param   CRe  Real part of product.
     * @param   CIm  Imag. part of product.
     */
    public final static void times( final double[][] ARe, final double[][] AIm, final double[][] BRe, final double[][] BIm, final double[][] CRe, final double[][] CIm )
    {
        double[][] tmpARe, tmpAIm, tmpBRe, tmpBIm;
        tmpARe = (ARe == CRe || ARe == CIm ) ? copy(ARe) : ARe;
        tmpAIm = (AIm == CRe || AIm == CIm ) ? copy(AIm) : AIm;
        tmpBRe = (BRe == CRe || BRe == CIm ) ? copy(BRe) : BRe;
        tmpBIm = (BIm == CRe || BIm == CIm ) ? copy(BIm) : BIm;
        final int n = tmpBRe.length;
        final int numRows = CRe.length;
        final int numCols = CRe[0].length;
        for(int i=0; i<numRows; i++)
        {
            final double[] rowARe = tmpARe[i];
            final double[] rowAIm = tmpAIm[i];
            final double[] rowRe  = CRe[i];
            final double[] rowIm  = CIm[i];
            for(int j=0; j<numCols; j++)
            {
                double wRe = 0;
                double wIm = 0;
                for(int k=0; k<n; k++)
                {
                    wRe += rowARe[k]*tmpBRe[k][j] - rowAIm[k]*tmpBIm[k][j];
                    wIm += rowARe[k]*tmpBIm[k][j] + rowAIm[k]*tmpBRe[k][j];
                }
                rowRe[j] = wRe;
                rowIm[j] = wIm;
            }
        }
    }

    /** Performs matrix multiplication: A * B = C.
     *	If A and/or B equals C a temporary copy of C is assigned,
     *	which should be avoided if run times are critical. 
     * <p>
     * @param  A  First multiplicand.
     * @param  B  Second multiplicand.
     * @param  C  Product of both.
     */
   public final static void times( final int [][] A, final int [][] B, final int [][] C ) {

	int [][] tmpA, tmpB;
	
	if( A == C ) {
	    tmpA = copy( C );
	    tmpB = B != C ? B : tmpA;
	} else {
	    tmpA = A;
	    tmpB = B != C ? B : copy( C );
	}

        final int n       = tmpB.length;
	final int numRows = C   .length; 
	final int numCols = C[0].length; 

	for(int i=0; i<numRows; i++) {
	  final int [] rowA = tmpA[i];
	  final int [] rowC =    C[i];

	  for(int j=0; j<numCols; j++) {
	    int sum = 0;

	    for(int k=0; k<n; k++) {
	      sum += rowA[k]*tmpB[k][j];
	    }

	    rowC[j] = sum;
	  }
	}
    }
    
    /** Performs matrix multiplication: A * B = C.
     *	If A and/or B equals C a temporary copy of C is assigned,
     *	which should be avoided if run times are critical. 
     * <p>
     * @param  A  First multiplicand.
     * @param  B  Second multiplicand.
     * @param  C  Product of both.
     */
    public final static void times( final double [][] A, final double [][] B, final double [][] C ) {

	double [][] tmpA, tmpB;
	
	if( A == C ) {
	    tmpA = copy( C );
	    tmpB = B != C ? B : tmpA;
	} else {
	    tmpA = A;
	    tmpB = B != C ? B : copy( C );
	}

        final int n       = tmpB.length;
	final int numRows = C   .length; 
	final int numCols = C[0].length; 

	for(int i=0; i<numRows; i++) {
	  final double [] rowA = tmpA[i];
	  final double [] rowC =    C[i];

	  for(int j=0; j<numCols; j++) {
	    double sum = 0;

	    for(int k=0; k<n; k++) {
	      sum += rowA[k]*tmpB[k][j];
	    }

	    rowC[j] = sum;
	  }
	}
    }
      
    /** Performs matrix multiplication: A * B = C.
     *	If B equals C a temporary copy of C is assigned,
     *	which should be avoided if run times are critical. 
     * <p>
     * @param  A   First multiplicand.
     * @param  B   Second multiplicand.
     * @param  C   Product of both. 
     */
    public final static void times( final int [][] A, final double [][] B, final double [][] C ) {

	double [][] tmpB = (B != C ? B : copy( C ) );

        final int n       = tmpB.length;
	final int numRows = C   .length; 
	final int numCols = C[0].length; 

	for(int i=0; i<numRows; i++) {
	  final int    [] rowA =    A[i];
	  final double [] rowC =    C[i];

	  for(int j=0; j<numCols; j++) {
	    double sum = 0;

	    for(int k=0; k<n; k++) {
	      sum += rowA[k]*tmpB[k][j];
	    }

	    rowC[j] = sum;
	  }
	}
    }
  

    /** Performs matrix multiplication: A * B = C.
     *	If A equals C a temporary copy of C is assigned,
     *	which should be avoided if run times are critical. 
     * <p>
     * @param  A   First multiplicand.
     * @param  B   Second multiplicand.
     * @param  C   Product of both.   
     */
    public final static void times( final double [][] A, final int [][] B, final double [][] C ) {

	double [][] tmpA = (A != C ? A : copy( C ) );

        final int n       = B   .length;
	final int numRows = C   .length; 
	final int numCols = C[0].length; 

	for(int i=0; i<numRows; i++) {
	  final double [] rowA = tmpA[i];
	  final double [] rowC =    C[i];

	  for(int j=0; j<numCols; j++) {
	    double sum = 0;

	    for(int k=0; k<n; k++) {
	      sum += rowA[k]*B[k][j];
	    }

	    rowC[j] = sum;
	  }
	}
    }

   
    /** Performs matrix multiplication: A * B = C.
     *
     * @param  A   First multiplicand.
     * @param  B   Second multiplicand.
     * @param  C   Product of both. 
     */
    public final static void times( final int [][] A, final int [][] B, final double [][] C ) {

        final int n       = B   .length;
	final int numRows = C   .length; 
	final int numCols = C[0].length; 

	for(int i=0; i<numRows; i++) {
	  final int    [] rowA =    A[i];
	  final double [] rowC =    C[i];

	  for(int j=0; j<numCols; j++) {
	    double sum = 0;

	    for(int k=0; k<n; k++) {
	      sum += rowA[k]*B[k][j];
	    }

	    rowC[j] = sum;
	  }
	}
    }

    /** Performs A / b = C. 
     *
     * @param   ARe    Real part of dividend.
     * @param   AIm    Imag. part of dividend.
     * @param   ARe    Real part of divisor.
     * @param   AIm    Imag. part of divisor.
     * @param   CRe    Real part of product.
     * @param   CIm    Imag. part of product.
     */
    public final static void divide( final double [][] ARe, final double [][] AIm,
				     final double bRe, final double bIm,
				     final double [][] CRe, final double [][] CIm ) {
	
	final double nn = bRe*bRe + bIm*bIm;
	final int numRows = CRe.length;
	final int numCols = CRe[0].length;
	
	for(int i=0; i<numRows; i++) {
	    final double[] rowMRe = ARe[i];
	    final double[] rowRe  = CRe[i];
	    final double[] rowMIm = AIm[i];
	    final double[] rowIm  = CIm[i];
		
	    for(int j=0; j<numCols; j++) {
		final double rr=rowMRe[j];
		final double ii=rowMIm[j];
		
		rowRe[j]=(rr*bRe+ii*bIm)/nn;
		rowIm[j]=(ii*bRe-rr*bIm)/nn;
	    }
	}	
    }

    /** Performs A * b = b * A = C.  
     * @param   ARe  Real part of first multiplicand.
     * @param   AIm  Imag. part of first multiplicand.
     * @param   bRe  Real part of second multiplicand.
     * @param   bIm  Imag. part of second multiplicand.
     * @param   CRe  Real part of product.
     * @param   CIm  Imag. part of product. */
    public final static void times( final double [][] ARe, final double [][] AIm,
				    final double bRe, final double bIm,
				    final double [][] CRe, final double [][] CIm ) {
	final int numRows = ARe.length;
	final int numCols = ARe[0].length;

	for(int i=0; i<numRows; i++) {

	    final double[] rowRe  = CRe[i];
	    final double[] rowIm  = CIm[i];
	    final double[] rowReA = ARe[i];
	    final double[] rowImA = AIm[i];
		
	    for(int j=0; j<numCols; j++) {
		final double rr=rowReA[j];
		final double ii=rowImA[j];
		
		rowRe[j] = rr*bRe - ii*bIm;
		rowIm[j] = rr*bIm + ii*bRe;
	    }
	}
    }

    /** Performs A * b = b * A = C.   
     * @param   A    First multiplicand.
     * @param   b    Scalar.
     * @param   C    Product. 
     */
    public final static void times( final int [][] A, final int b, final int [][] C ) {

	final int numRows = C   .length; 
	final int numCols = C[0].length; 

	for(int i=0; i<numRows; i++) {
	  final int [] rowA = A[i];
	  final int [] rowC = C[i];

	  for(int j=0; j<numCols; j++)
	      rowC[j] = rowA[j] * b;
	}
    }

    /** Performs A / b = C. 
     * 
     * @param A  Dividend.
     * @param b  Divisor.
     * @param C  Product.     
     */
    public final static void divide( final int [][] A, final int b, final int [][] C ) {

	final int numRows = C   .length; 
	final int numCols = C[0].length; 

	for(int i=0; i<numRows; i++) {
	  final int [] rowA = A[i];
	  final int [] rowC = C[i];

	  for(int j=0; j<numCols; j++)
	      rowC[j] = rowA[j] / b;
	}
    }

    /** Performs A / b = C. 
     * 
     * @param A  Dividend.
     * @param b  Divisor.
     * @param C  Product.     
     */
    public final static void divide( final int [][] A, final double b, final double [][] C ) {

	final int numRows = C   .length; 
	final int numCols = C[0].length; 

	for(int i=0; i<numRows; i++) {
	  final int [] rowA = A[i];
	  final double [] rowC = C[i];

	  for(int j=0; j<numCols; j++)
	      rowC[j] = ((double)rowA[j]) / b;
	}
    }

    /** Performs b / A = C. 
     * 
     * @param A  Dividend.
     * @param b  Divisor.
     * @param C  Product.     
     */
    public final static void divide( final int b, final int [][] A, final int [][] C ) {

	final int numRows = C   .length; 
	final int numCols = C[0].length; 

	for(int i=0; i<numRows; i++) {
	  final int [] rowA = A[i];
	  final int [] rowC = C[i];

	  for(int j=0; j<numCols; j++)
	      rowC[j] =  b / rowA[j];
	}
    }

    /** Performs b / A = C. 
     * 
     * @param A  Dividend.
     * @param b  Divisor.
     * @param C  Product.     
     */
    public final static void divide( final double b, final double [][] A, final double [][] C ) {

	final int numRows = C   .length; 
	final int numCols = C[0].length; 

	for(int i=0; i<numRows; i++) {
	  final double [] rowA = A[i];
	  final double [] rowC = C[i];

	  for(int j=0; j<numCols; j++)
	      rowC[j] =  b / rowA[j];
	}
    }

    /** Performs b / (Are+i*Aim) = (Cre+i*Cim). 
     * 
     * @param A  Dividend.
     * @param b  Divisor.
     * @param C  Product.     
     */
    public final static void divide( final double b, final double [][] Are, final double [][]Aim,
				     final double [][] Cre, final double [][] Cim ) {

	final int numRows = Cre   .length; 
	final int numCols = Cre[0].length; 

	for(int i=0; i<numRows; i++) {
	  final double [] rowAre = Are[i];
	  final double [] rowCre = Cre[i];
	  final double [] rowAim = Aim[i];
	  final double [] rowCim = Cim[i];
	  
	  for(int j=0; j<numCols; j++){
	      final double Re = rowAre[j];
	      final double Im = rowAim[j];
	      final double nn = Re*Re+Im*Im;
	      if ( nn == 0. ){
		  rowCre[j] = b /nn;
		  rowCim[j] = 0.;
	      }
	      else
		  {
		      rowCre[j] =  b * Re /nn;
		      rowCim[j] = -b * Im /nn; 
		  }
	  }
	}
    }

     /** Performs b / A = C. 
     * 
     * @param A  Dividend.
     * @param b  Divisor.
     * @param C  Product.     
     */
    public final static void divide( final double b, final int [][] A, final double [][] C ) {

	final int numRows = C   .length; 
	final int numCols = C[0].length; 

	for(int i=0; i<numRows; i++) {
	  final int [] rowA = A[i];
	  final double [] rowC = C[i];

	  for(int j=0; j<numCols; j++)
	      rowC[j] =  b / rowA[j];
	}
    }

   /** Performs A / b = C. 
     * 
     * @param A  Dividend.
     * @param b  Divisor.
     * @param C  Product.     
      */
    public final static void divide( final double [][] A, final double b, final double [][] C ) {

	final int numRows = C   .length; 
	final int numCols = C[0].length; 

	for(int i=0; i<numRows; i++) {
	  final double [] rowA = A[i];
	  final double [] rowC = C[i];

	  for(int j=0; j<numCols; j++)
	      rowC[j] = rowA[j] / b;
	}
    }

    /** Performs A * b = b * A = C.    
     *  @param   A    First multiplicand.
     *  @param   b    Scalar.
     *  @param   C    Product. 
     */
    public final static void times( final double [][] A, final double b, final double [][] C ) {

	final int numRows = C   .length; 
	final int numCols = C[0].length; 

	for(int i=0; i<numRows; i++) {
	  final double [] rowA = A[i];
	  final double [] rowC = C[i];

	  for(int j=0; j<numCols; j++)
	      rowC[j] = rowA[j] * b;
	}
    }

    /** Performs A * b = b * A = C.    
     * @param   A    First multiplicand.
     * @param   b    Scalar.
     * @param   C    Product. 
     */
    public final static void times( final int [][] A, final double b, final double [][] C ) {

	final int numRows = C   .length; 
	final int numCols = C[0].length; 

	for(int i=0; i<numRows; i++) {
	  final int    [] rowA = A[i];
	  final double [] rowC = C[i];

	  for(int j=0; j<numCols; j++)
	      rowC[j] = rowA[j] * b;
	}
    }

    /** Performs matrix - vector multiplication: A * b = c.
     *	If b equals c a temporary clone of c is created,
     *	which should be avoided if run times are critical. 
     * @param MRe Real part of matrix.
     * @param MIm Imag. part of matrix.
     * @param VRe Real part of vector.
     * @param VIm Imag. part of vector.
     * @param WRe Real part of product.
     * @param WIm Imag. part of product.
     */
    public final static void times(final double[][] MRe, final double[][] MIm, final double[] VRe, final double[] VIm, final double[] WRe, final double[] WIm)
    {
        double[] tmpVRe, tmpVIm;
        tmpVRe = (VRe == WRe || VRe == WIm ) ? (double[]) VRe.clone() : VRe;
        tmpVIm = (VIm == WRe || VIm == WIm ) ? (double[]) VIm.clone() : VIm;
	final int numOfRows = MRe.length;
	final int numOfCols = MRe[0].length;
	for(int i = 0; i < numOfRows; i++)
        {
	    double re = 0.;
	    double im = 0.;
	    double[] rowRe = MRe[i];
	    double[] rowIm = MIm[i];
	    for(int j = 0; j < numOfCols; j++)
            {
                re += rowRe[j] * tmpVRe[j] - rowIm[j] * tmpVIm[j];
                im += rowRe[j] * tmpVIm[j] + rowIm[j] * tmpVRe[j];
            }
            WRe[i] = re;
            WIm[i] = im;
        }	
    }

    /** Performs matrix - vector multiplication: b * A = c.
     *	If b equals c a temporary clone of c is created,
     *	which should be avoided if run times are critical. 
     * @param MRe Real part of matrix.
     * @param MIm Imag. part of matrix.
     * @param VRe Real part of vector.
     * @param VIm Imag. part of vector.
     * @param WRe Real part of product.
     * @param WIm Imag. part of product.
     */
    public final static void times( final double [] VRe, final double [] VIm,
			     final double [][] MRe ,final double [][] MIm,
 			     final double [] WRe, final double [] WIm)
    {
	final int numOfRows = MRe.length;
	final int numOfCols = MRe[0].length;
	
	for(int i=0; i<numOfCols; i++) {
	    double thisRe=0;
	    double thisIm=0;

	    for(int j=0; j<numOfRows; j++)
		{
		    thisRe += MRe[j][i]*VRe[j] - MIm[j][i]*VIm[j];
		    thisIm += MRe[j][i]*VIm[j] + MIm[j][i]*VRe[j];
		}
	    WRe[i]=thisRe;
	    WIm[i]=thisIm;
	}
    }

    
    /** Performs matrix - vector multiplication: A * b = c.
     *	If b equals c a temporary clone of c is created,
     *	which should be avoided if run times are critical 
     * @param A Matrix.
     * @param b Vector.
     * @param c Product.
     */
    public final static void times( final int [][] A, final int [] b, final int [] c ) {

	int [] tmpB = b != c ? b : VectorOperations.copy( c );

        final int n       = b.length;
	final int numRows = A.length;

	for(int i=0; i<numRows; i++) {
	  final int [] rowA = A[i];

	  int sum = 0;

	  for(int j=0; j<n; j++) {
	      sum += rowA[j]*tmpB[j];
	  }
	  c[i] = sum;
	}
    }


    /** Performs matrix - vector multiplication: A * b = c. 
     * @param A Matrix.
     * @param b Vector.
     * @param c Product.
     */
    public final static void times( final int [][] A, final int [] b, final double [] c ) {

        final int n       = b.length;
	final int numRows = A.length;

	for(int i=0; i<numRows; i++) {
	  final int [] rowA = A[i];

	  int sum = 0;

	  for(int j=0; j<n; j++) {
	      sum += rowA[j]*b[j];
	  }

	  c[i] = sum;
	}
    }

    /** Performs matrix - vector multiplication: A * b = c.
     *	If b equals c a temporary clone of c is created,
     *	which should be avoided if run times are critical. 
     * @param A Matrix.
     * @param b Vector.
     * @param c Product. 
     */
    public final static void times( final double [][] A, final double [] b, final double [] c ) {

	double [] tmpB = b != c ? b : VectorOperations.copy( c );

        final int n       = b.length;
	final int numRows = A.length;

	for(int i=0; i<numRows; i++) {
	  final double [] rowA = A[i];

	  double sum = 0;

	  for(int j=0; j<n; j++) {
	      sum += rowA[j]*tmpB[j];
	  }

	  c[i] = sum;
	}
    }


    /** Performs matrix - vector multiplication: A * b = c.
     *	If b equals c a temporary clone of c is created,
     *	which should be avoided if run times are critical. 
     * @param A Matrix.
     * @param b Vector.
     * @param c Product. 
     */
    public final static void times( final int [][] A, final double [] b, final double [] c ) {

	double [] tmpB = b != c ? b : VectorOperations.copy( c );

        final int n       = b.length;
	final int numRows = A.length;

	for(int i=0; i<numRows; i++) {
	  final int [] rowA = A[i];

	  double sum = 0;

	  for(int j=0; j<n; j++) {
	      sum += rowA[j]*tmpB[j];
	  }

	  c[i] = sum;
	}
    }


    /** Performs matrix - vector multiplication: A * b = c. 
     * @param A Matrix.
     * @param b Vector.
     * @param c Product.
     */
    public final static void times( final double [][] A, final int [] b, final double [] c ) {

        final int n       = b.length;
	final int numRows = A.length;

	for(int i=0; i<numRows; i++) {
	  final double [] rowA = A[i];

	  double sum = 0;

	  for(int j=0; j<n; j++) {
	      sum += rowA[j]*b[j];
	  }

	  c[i] = sum;
	}
    }


    /** Performs vector - matrix multiplication: a^t * B = c^t.
     *	If a equals c a temporary clone of c is created,
     *	which should be avoided if run times are critical. 
     * @param B Matrix.
     * @param a Vector.
     * @param c Product.
     */
    public final static void times( final int [] a, final int [][] B, final int [] c ) {

	int [] tmpA = a != c ? a : VectorOperations.copy( c );

        final int n       = a.length;
	final int numCols = c.length;

	for(int i=0; i<numCols; i++) {

	  int sum = 0;

	  for(int j=0; j<n; j++) {
	      sum += tmpA[j]*B[j][i];
	  }

	  c[i] = sum;
	}
    }

    /** Performs vector - matrix multiplication: a^t * B = c^t. 
     * @param B Matrix.
     * @param a Vector.
     * @param c Product.
     */
    public final static void times( final int [] a, final int [][] B, final double [] c ) {

        final int n       = a.length;
	final int numCols = c.length;

	for(int i=0; i<numCols; i++) {

	  int sum = 0;

	  for(int j=0; j<n; j++) {
	      sum += a[j]*B[j][i];
	  }

	  c[i] = sum;
	}
    }

    /** Performs vector - matrix multiplication: a^t * B = c^t.
     *	If a equals c a temporary clone of c is created,
     *	which should be avoided if run times are critical. 
     * @param B Matrix.
     * @param a Vector.
     * @param c Product.
     */
    public final static void times( final double [] a, final double [][] B, final double [] c ) {

	double [] tmpA = a != c ? a : VectorOperations.copy( c );

        final int n       = a.length;
	final int numCols = c.length;

	for(int i=0; i<numCols; i++) {

	  double sum = 0;

	  for(int j=0; j<n; j++) {
	      sum += tmpA[j]*B[j][i];
	  }

	  c[i] = sum;
	}
    }


    /** Performs vector - matrix multiplication: a^t * B = c^t. 
     * @param B Matrix.
     * @param a Vector.
     * @param c Product.
     */
    public final static void times( final int [] a, final double [][] B, final double [] c ) {

        final int n       = a.length;
	final int numCols = c.length;

	for(int i=0; i<numCols; i++) {

	  double sum = 0;

	  for(int j=0; j<n; j++) {
	      sum += a[j]*B[j][i];
	  }

	  c[i] = sum;
	}
    }


    /** Performs vector - matrix multiplication: a^t * B = c^t.
     *	If a equals c a temporary clone of c is created,
     *	which should be avoided if run times are critical.
     * @param B Matrix.
     * @param a Vector.
     * @param c Product.
     */
    public final static void times( final double [] a, final int [][] B, final double [] c ) {

	double [] tmpA = a != c ? a : VectorOperations.copy( c );

        final int n       = a.length;
	final int numCols = c.length;

	for(int i=0; i<numCols; i++) {

	  double sum = 0;

	  for(int j=0; j<n; j++) {
	      sum += tmpA[j]*B[j][i];
	  }

	  c[i] = sum;
	}
    }


    /** Performs vector - vector multiplication: a * b^t = C.
     * @param  a  First Vector.
     * @param  b  Second Vector.
     * @param  C  The resulted matrix.
     */
    public final static void times( final int [] a, final int [] b, final int [][] C ) {

        final int n = a.length;

	for(int i=0; i<n; i++)
	    for(int j=0; j<n; j++) {
		C[i][j] = a[i] * b[j];
	    }
    }


    /** Performs vector - vector multiplication: a * b^t = C.
     * @param  a  First Vector.
     * @param  b  Second Vector.
     * @param  C  The resulted matrix.
     */
    public final static void times( final double [] a, final double [] b, final double [][] C ) {

        final int n = a.length;

	for(int i=0; i<n; i++)
	    for(int j=0; j<n; j++) {
		C[i][j] = a[i] * b[j];
	    }
    }


    /** Performs vector - vector multiplication: a * b^t = C.
     * @param  a  First Vector.
     * @param  b  Second Vector.
     * @param  C  The resulted matrix.
     */
    public final static void times( final int [] a, final double [] b, final double [][] C ) {

        final int n = a.length;

	for(int i=0; i<n; i++)
	    for(int j=0; j<n; j++) {
		C[i][j] = a[i] * b[j];
	    }
    }


    /** Performs vector - vector multiplication: a * b^t = C. 
     * @param  a  First Vector.
     * @param  b  Second Vector.
     * @param  C  The resulted matrix.
     */
    public final static void times( final double [] a, final int [] b, final double [][] C ) {

        final int n = a.length;

	for(int i=0; i<n; i++)
	    for(int j=0; j<n; j++) {
		C[i][j] = a[i] * b[j];
	    }
    }





    /* ************************ */
    /* ******** plus ********** */
    /* ************************ */

    /** Performs: A + B = C .
     * @param A  First summend. 
     * @param B  Second summend. 
     * @param C  Sum. 
     */
    public final static void plus( final int [][] A, final int [][] B, final int [][] C ) {

	final int numRows = C   .length; 
	final int numCols = C[0].length; 

	for(int i=0; i<numRows; i++) {
	  final int    [] rowA = A[i];
	  final int    [] rowB = B[i];
	  final int    [] rowC = C[i];

	  for(int j=0; j<numCols; j++)
	    rowC[j] = rowA[j] + rowB[j];

	}
    }
    
    /** Performs: A + B = C.
     * @param A  First summend. 
     * @param B  Second summend. 
     * @param C  Sum. 
     */
    public final static void plus( final double [][] A, final double [][] B, final double [][] C ) {

	final int numRows = C   .length; 
	final int numCols = C[0].length; 

	for(int i=0; i<numRows; i++) {
	  final double [] rowA = A[i];
	  final double [] rowB = B[i];
	  final double [] rowC = C[i];

	  for(int j=0; j<numCols; j++)
	    rowC[j] = rowA[j] + rowB[j];

	}
    }
    
    /** Performs: A + B = C .
     * @param A  First summend. 
     * @param B  Second summend. 
     * @param C  Sum. 
     */
    public final static void plus( final double [][] A, final int [][] B, final double [][] C ) {

	final int numRows = C   .length; 
	final int numCols = C[0].length; 

	for(int i=0; i<numRows; i++) {
	  final double [] rowA = A[i];
	  final int    [] rowB = B[i];
	  final double [] rowC = C[i];

	  for(int j=0; j<numCols; j++)
	    rowC[j] = rowA[j] + rowB[j];

	}
    }

    /** Performs: A + B = C .
     * @param A  First summend. 
     * @param B  Second summend. 
     * @param C  Sum. 
     */
    public final static void plus( final int [][] A, final double [][] B, final double [][] C ) {
	plus( B, A, C );
    }
    
    /** Performs: A + b = C .
     * @param A  First summend. 
     * @param b  Second summend. 
     * @param C  Sum. 
     */
    public final static void plus( final int [][] A, final int b, final int [][] C ) {

	final int numRows = C   .length; 
	final int numCols = C[0].length; 

	for(int i=0; i<numRows; i++) {
	  final int    [] rowA = A[i];
	  final int    [] rowC = C[i];

	  for(int j=0; j<numCols; j++)
	    rowC[j] = rowA[j] + b;

	}
    }
    
    /** Performs: A + b = C .
     * @param A  First summend. 
     * @param b  Second summend. 
     * @param C  Sum. 
     */
    public final static void plus( final double [][] A, final double b, final double [][] C ) {

	final int numRows = C   .length; 
	final int numCols = C[0].length; 

	for(int i=0; i<numRows; i++) {
	  final double [] rowA = A[i];
	  final double [] rowC = C[i];

	  for(int j=0; j<numCols; j++)
	    rowC[j] = rowA[j] + b;

	}
    }
    
    /** Performs: A + b = C .
     * @param A  First summend. 
     * @param b  Second summend. 
     * @param C  Sum. 
    */
    public final static void plus( final int [][] A, final double b, final double [][] C ) {

	final int numRows = C   .length; 
	final int numCols = C[0].length; 

	for(int i=0; i<numRows; i++) {
	  final int    [] rowA = A[i];
	  final double [] rowC = C[i];

	  for(int j=0; j<numCols; j++)
	    rowC[j] = rowA[j] + b;

	}
    }

    /* ************************ */
    /* ******** minus ********* */
    /* ************************ */

    /** Performs: A - B = C.
     *  @param  B   The subtrahend.
     *  @param  C   The Product.
     */
    public final static void minus( final int [][] A, final int [][] B, final int [][] C ) {

	final int numRows = C   .length; 
	final int numCols = C[0].length; 

	for(int i=0; i<numRows; i++) {
	  final int    [] rowA = A[i];
	  final int    [] rowB = B[i];
	  final int    [] rowC = C[i];

	  for(int j=0; j<numCols; j++)
	    rowC[j] = rowA[j] - rowB[j];

	}
    }
    
    /** Performs: A - B = C .
     *  @param  B   The subtrahend.
     *  @param  C   The Product.
     */
    public final static void minus( final double [][] A, final double [][] B, final double [][] C ) {

	final int numRows = C   .length; 
	final int numCols = C[0].length; 

	for(int i=0; i<numRows; i++) {
	  final double [] rowA = A[i];
	  final double [] rowB = B[i];
	  final double [] rowC = C[i];

	  for(int j=0; j<numCols; j++)
	    rowC[j] = rowA[j] - rowB[j];

	}
    }
    
    /** Performs: A - B = C.
     *  @param  B   The subtrahend.
     *  @param  C   The Product.
     */
    public final static void minus( final double [][] A, final int [][] B, final double [][] C ) {

	final int numRows = C   .length; 
	final int numCols = C[0].length; 

	for(int i=0; i<numRows; i++) {
	  final double [] rowA = A[i];
	  final int    [] rowB = B[i];
	  final double [] rowC = C[i];

	  for(int j=0; j<numCols; j++)
	    rowC[j] = rowA[j] - rowB[j];

	}
    }
    
    /** Performs: A - B = C .
     *  @param  B   The subtrahend.
     *  @param  C   The Product.
     */
    public final static void minus( final int [][] A, final double[][] B, final double [][] C ) {

	final int numRows = C   .length; 
	final int numCols = C[0].length; 

	for(int i=0; i<numRows; i++) {
	  final int    [] rowA = A[i];
	  final double [] rowB = B[i];
	  final double [] rowC = C[i];

	  for(int j=0; j<numCols; j++)
	    rowC[j] = rowA[j] - rowB[j];

	}
    }

    /** Performs: A - b = C .
     *  @param  b   The subtrahend.
     *  @param  C   The Product.
     */
    public final static void minus( final int [][] A, final int b, final int [][] C ) {

	final int numRows = C   .length; 
	final int numCols = C[0].length; 

	for(int i=0; i<numRows; i++) {
	  final int    [] rowA = A[i];
	  final int    [] rowC = C[i];

	  for(int j=0; j<numCols; j++)
	    rowC[j] = rowA[j] - b;

	}
    }

     /** Performs: a - B = C .
     *  @param  B   The subtrahend.
     *  @param  C   The Product.
     */
    public final static void minus( final int a, final int[][] B, final int [][] C ) {

	final int numRows = C   .length; 
	final int numCols = C[0].length; 

	for(int i=0; i<numRows; i++) {
	  final int    [] rowB = B[i];
	  final int    [] rowC = C[i];

	  for(int j=0; j<numCols; j++)
	    rowC[j] = a - rowB[j];

	}
    }
    
    /** Performs: A - b = C .
     *  @param  b   The subtrahend.
     *  @param  C   The Product.
     */
    public final static void minus( final double [][] A, final double b, final double [][] C ) {

	final int numRows = C   .length; 
	final int numCols = C[0].length; 

	for(int i=0; i<numRows; i++) {
	  final double [] rowA = A[i];
	  final double [] rowC = C[i];

	  for(int j=0; j<numCols; j++)
	    rowC[j] = rowA[j] - b;

	}
    }

     /** performs: a - B = C .
     *  @param  B   The subtrahend.
     *  @param  C   The Product.
     */
    public final static void minus( final double a, final double [][] B, final double [][] C ) {

	final int numRows = C   .length; 
	final int numCols = C[0].length; 

	for(int i=0; i<numRows; i++) {
	  final double [] rowB = B[i];
	  final double [] rowC = C[i];

	  for(int j=0; j<numCols; j++)
	    rowC[j] = a - rowB[j];

	}
    }
    
    /** Performs: A - b = C .
     *  @param  b   The subtrahend.
     *  @param  C   The Product.
     */
    public final static void minus( final int [][] A, final double b, final double [][] C ) {

	final int numRows = C   .length; 
	final int numCols = C[0].length; 

	for(int i=0; i<numRows; i++) {
	  final int    [] rowA = A[i];
	  final double [] rowC = C[i];

	  for(int j=0; j<numCols; j++)
	    rowC[j] = rowA[j] - b;

	}
    }

    /** Performs: A - B = C .
     *  @param  B   The subtrahend.
     *  @param  C   The Product.
     */
    public final static void minus( final int a, final double [][] B, final double [][] C ) {

	final int numRows = C   .length; 
	final int numCols = C[0].length; 

	for(int i=0; i<numRows; i++) {
	  final double [] rowB = B[i];
	  final double [] rowC = C[i];

	  for(int j=0; j<numCols; j++)
	    rowC[j] = a - rowB[j];

	}
    }

    /** Assign matrix  B to matrix  A. **/
    public final static void assign( final double [][] A, final double [][] B ) {
	final int numRows = A   .length; 
	final int numCols = A[0].length; 
	for(int i=0; i<numRows; i++)
	    System.arraycopy( A[i], 0, B[i], 0, numCols );
    }

    /** Assign matrix  B to matrix  A. **/
    public final static void assign( final int [][] A, final int [][] B ) {
	final int numRows = A   .length; 
	final int numCols = A[0].length; 
	for(int i=0; i<numRows; i++)
	    System.arraycopy( A[i], 0, B[i], 0, numCols );
    }

    /** Assign matrix  B to matrix  A. **/
    public final static void assign( final int [][] A, final double [][] B ) {
	final int numRows = A   .length; 
	final int numCols = A[0].length; 

	for(int i=0; i<numRows; i++) {
	    final int    [] rowA = A[i];
	    final double [] rowB = B[i];
	    
	    for(int j=0; j<numCols; j++)
		rowB[j] = rowA[j];
	}
    }
    
    /** Returns a copy of matrix  A. 
     * @param   A    The matrix to copy.
     * @return  The copy-matrix.
     */
    public final static int[][] copy( final int [][] A ) {
	int [][] cloneOfA = new int[ A.length ][ A[0].length ];
	assign( A, cloneOfA );
	return cloneOfA;
    }

    /** Returns a copy of matrix  A. 
     * @param   A    The matrix to copy.
     * @return  The copy-matrix.
     */
    public final static double[][] copy( final double [][] A ) {
	double [][] cloneOfA = new double[ A.length ][ A[0].length ];
	assign( A, cloneOfA );
	return cloneOfA;
    }

    /** Assigns all entries of matrix  A with  <code>value</code>.
     * @param   A        The matrix to assign.
     * @param   value    The assigned value.
     */
    public final static void assign( final double [][] A, final double value ) {

	final int numRows = A   .length;
	final int numCols = A[0].length; 

	final double []tmp = A[0];
	
	for(int j=0;j<numCols;j++) // erste Reihe erzeugen
	    tmp[j]=value;

	for(int i=1;i<numRows;i++) // copy them to the rest
	    System.arraycopy( tmp, 0, A[i], 0, numCols );
    }

    /** Assigns all entries of matrix  A with  <code>value</code>.
     * @param   A        The matrix to assign.
     * @param   value    The assigned value.
     */
    public final static void assign( final int [][] A, final int value ) {

	final int numRows = A   .length; 
	final int numCols = A[0].length; 

	final int []tmp = A[0];
	
	for(int j=0;j<numCols;j++) // erste Reihe erzeugen
	    tmp[j]=value;

	for(int i=1;i<numRows;i++) // copy them to the rest
	    System.arraycopy( tmp, 0, A[i], 0, numCols );
    }

    /** Makes matrix A to id.
     * @param   A        The matrix to assign.
     * @note the matrix must be quadratic.
     */
    public final static void assignId( final int [][] A ){
	final int length = A.length;
	
	int []tmp = new int[length];

	for( int i=0; i<length; i++)
	    for( int j=0; j<length; j++){
		tmp = A[i];
		if( i==j ) tmp[j] = 1;
		else tmp[j] = 0;
	    }
    }

    /** Makes matrix A to id.
     * @param   A        The matrix to assign.
     * @note the matrix must be quadratic.
     */
    public final static void assignId( final double [][] A ){
	final int length = A.length;
	
	double []tmp = new double[length];

	for( int i=0; i<length; i++)
	    for( int j=0; j<length; j++){
		tmp = A[i];
		if( i==j ) tmp[j] = 1.;
		else tmp[j] = 0.;
	    }
    }

    /** Assigns all entries of matrix  A with  <code>zero</code>.
     * @param   A        The matrix to assign.
     */
    public final static void assignZero( final int [][] A ){

	final int numRows = A   .length; 
	final int numCols = A[0].length; 

	final int []tmp = A[0];
	
	for(int j=0;j<numCols;j++) // erste Reihe erzeugen
	    tmp[j]=0;

	for(int i=1;i<numRows;i++) // copy them to the rest
	    System.arraycopy( tmp, 0, A[i], 0, numCols );
    }

    /** Assigns all entries of matrix  A with  <code>zero</code>.
     * @param   A        The matrix to assign.
     */
    public final static void assignZero( final double [][] A ){

	final int numRows = A   .length; 
	final int numCols = A[0].length; 

	final double []tmp = A[0];
	
	for(int j=0;j<numCols;j++) // erste Reihe erzeugen
	    tmp[j]=0;

	for(int i=1;i<numRows;i++) // copy them to the rest
	    System.arraycopy( tmp, 0, A[i], 0, numCols );
    }

    /** Stores the diagonal of matrix  A to <code>v</code>.
     * @param   A        The matrix.
     * @param   v        The vector.
     */
    public final static void getDiagonal( final int [][] A, final int []v ){
	final int length = v.length;

	for( int i=0; i<length; i++)
	     v[i] = A[i][i] ;
    }

    /** Assigns the diagonal of matrix  A with  <code>v</code>.
     * @param   A        The matrix to assign.
     * @param   v        The vector.
     */
    public final static void assignDiagonal( final int [][] A, final int []v ){
	final int length = v.length;

	for( int i=0; i<length; i++)
	    A[i][i] = v[i];
    }

    /** Assigns the diagonal of matrix  A with  <code>v</code>.
     * @param   A        The matrix to assign.
     * @param   v        The value.
     */
    public final static void assignDiagonal( final int [][] A, final int v ){
	final int length = Math.min(A.length, A[0].length );
	for( int i=0; i<length; i++)
	    A[i][i] = v;
    }

      /** Assigns the diagonal of matrix  A with  <code>v</code>.
     * @param   A        The matrix to assign.
     * @param   v        The value.
     */
    public final static void assignDiagonal( final double [][] A, final double v ){
	final int length = Math.min(A.length, A[0].length );
	for( int i=0; i<length; i++)
	    A[i][i] = v;
    }

    /** Assigns the diagonal of matrix  A with  <code>v</code>.
     * @param   A        The matrix to assign.
     * @param   v        The vector.
     */
    public final static void assignDiagonal( final double [][] A, final double []v ){
	final int length = v.length;

	for( int i=0; i<length; i++)
	    A[i][i] = v[i];
    }

     /** Stores the diagonal of matrix  A to <code>v</code>.
     * @param   A        The matrix.
     * @param   v        The vector.
     */
    public final static void getDiagonal( final double [][] A, final double []v ){
	final int length = v.length;

	for( int i=0; i<length; i++)
	    v[i] = A[i][i];
    }

   /** Assigns the diagonal of matrix  A with  <code>v</code>.
     * @param   A        The matrix to assign.
     * @param   v        The vector.
     */
    public final static void assignDiagonal( final double [][] A, final int []v ){
	final int length = v.length;

	for( int i=0; i<length; i++)
	    A[i][i] = (double)v[i];
    }

    /** Assigns the row of matrix A with vector v.
     * @param   A        The matrix to assign.
     * @param   v        The vector.
     * @param   rowNum   The number or row to assign.
     */ 
    public final static void assignRow( final int [][] A, final int []v, int rowNum ){
	System.arraycopy(v, 0, A[rowNum], 0, v.length);
    }

    /**Stores the row of matrix A in vector v.
     * @param   A        The matrix.
     * @param   v        The vector.
     * @param   rowNum   The number or row to assign.
     */ 
    public final static void getRow( final int [][] A, final int []v, int rowNum ){
	System.arraycopy(A[rowNum], 0, v, 0, v.length);
    }
    
     /**Stores the row of matrix A in vector v.
     * @param   A        The matrix.
     * @param   v        The vector.
     * @param   rowNum   The number or row to assign.
     */ 
    public final static void getRow( final double [][] A, final double []v, int rowNum ){
	System.arraycopy(A[rowNum], 0, v, 0, v.length);
    }

   /** Assigns the row of matrix A with vector v.
     * @param   A        The matrix to assign.
     * @param   v        The vector.
     * @param   rowNum   The number or row to assign.
     */  
    public final static void assignRow( final double [][] A, final double[]v, int rowNum ){
	 System.arraycopy(v, 0, A[rowNum], 0, v.length);
    }
    
    /** Assigns the row of matrix A with vector v.
     * @param   A        The matrix to assign.
     * @param   v        The vector.
     * @param   rowNum   The number or row to assign.
     */  
    public final static void assignRow( final double [][] A, final int[]v, int rowNum ){
	final int l = v.length;
	final double [] row = A[rowNum];
	for( int i=0; i<l; i++)
	    row[i] = (double) v[i];
    }

  /** Assigns the column of matrix A with vector v.
     * @param   A        The matrix to assign.
     * @param   v        The vector.
     * @param   rowNum   The number or column to assign.
     */ 
    public final static void assignCol( final int [][] A, final int []v, int colNum ){
	final int l = v.length;
	for( int i=0; i<l; i++ )
	    A[i][colNum] = v[i];
    }
    
  /** Assigns the column of matrix A in vector v.
     * @param   A        The matrix.
     * @param   v        The vector.
     * @param   rowNum   The number or column to assign.
     */ 
    public final static void getCol( final int [][] A, final int []v, int colNum ){
	final int l = v.length;
	for( int i=0; i<l; i++ )
	    v[i] = A[i][colNum];
    }

  /** Assigns the column of matrix A in vector v.
     * @param   A        The matrix.
     * @param   v        The vector.
     * @param   rowNum   The number or column to assign.
     */ 
    public final static void getCol( final double [][] A, final double []v, int colNum ){
	final int l = v.length;
	for( int i=0; i<l; i++ )
	    v[i] = A[i][colNum];
    }

   /** Assigns the column of matrix A with vector v.
     * @param   A        The matrix to assign.
     * @param   v        The vector.
     * @param   rowNum   The column or row to assign.
     */  
    public final static void assignCol( final double [][] A, final double[]v, int colNum ){
	final int l = v.length;
	for( int i=0; i<l; i++ )
	    A[i][colNum] = v[i];
    }
    
    /** Assigns the column of matrix A with vector v.
     * @param   A        The matrix to assign.
     * @param   v        The vector.
     * @param   rowNum   The number or column to assign.
     */  
    public final static void assignCol( final double [][] A, final int[]v, int colNum ){
	final int l = v.length;
	for( int i=0; i<l; i++ )
	    A[i][colNum] = (double)v[i];
    }
 
    /** Rounds matrix  A and stores result in  B. 
     ** A and B may coinside. 
     * @param  A  The matrix to round.
     * @param  B  The rounded matrix.
     */
    public final static void round( final double [][] A, final double [][] B ) {

	final int numRows = A   .length; 
	final int numCols = A[0].length; 

	for(int i=0; i<numRows; i++) {
	    final double    [] rowA =    A[i];
	    final double    [] rowB =    B[i];
	    
	    for(int j=0; j<numCols; j++)
		rowB[j] = Math.floor( rowA[j] + 0.5 );
	}
    }

    
    /** Rounds matrix  A and stores result in  B. 
     ** A and B may coinside. 
     * @param  A  The matrix to round.
     * @param  B  The rounded matrix.
     */
    public final static void round( final double [][] A, final int [][] B ) {

	final int numRows = A   .length; 
	final int numCols = A[0].length; 

	for(int i=0; i<numRows; i++) {
	    final double    [] rowA =    A[i];
	    final int       [] rowB =    B[i];
	    
	    for(int j=0; j<numCols; j++)
		rowB[j] = (int)Math.floor( rowA[j] + 0.5 );
	}
    }

    /** Floors matrix  A and stores result in  B. 
     ** A and B may coinside. 
     * @param  A  The matrix to round.
     * @param  B  The rounded matrix.
     */
    public final static void floor( final double [][] A, final double [][] B ) {

	final int numRows = A   .length; 
	final int numCols = A[0].length; 

	for(int i=0; i<numRows; i++) {
	    final double    [] rowA =    A[i];
	    final double    [] rowB =    B[i];
	    
	    for(int j=0; j<numCols; j++)
		rowB[j] = Math.floor( rowA[j] );
	}
    }

    
    /** Floors matrix  A and stores result in  B. 
     ** A and B may coinside. 
     * @param  A  The matrix to round.
     * @param  B  The rounded matrix.
     */
    public final static void floor( final double [][] A, final int [][] B ) {

	final int numRows = A   .length; 
	final int numCols = A[0].length; 

	for(int i=0; i<numRows; i++) {
	    final double    [] rowA =    A[i];
	    final int       [] rowB =    B[i];
	    
	    for(int j=0; j<numCols; j++)
		rowB[j] = (int)Math.floor( rowA[j] );
	}
    }

    /** Negates matrix A and stores result in B . 
     ** A and B may coinside. 
     * @param  A   The matrix to negate.
     * @param  B   The negated matrix.
     */
    public final static void neg( final int [][] A, final int [][] B ) {

	final int numRows = A   .length; 
	final int numCols = A[0].length; 

	for(int i=0; i<numRows; i++) {
	    final int    [] rowA =    A[i];
	    final int    [] rowB =    B[i];
	    
	    for(int j=0; j<numCols; j++)
		rowB[j] = -rowA[j];
	}
    }

    /** Negates matrix A and stores result in B . 
     ** A and B may coinside. 
     * @param  A   The matrix to negate.
     * @param  B   The negated matrix.
     */
    public final static void neg( final double [][] A, final double [][] B ) {

	final int numRows = A   .length; 
	final int numCols = A[0].length; 

	for(int i=0; i<numRows; i++) {
	    final double    [] rowA =    A[i];
	    final double    [] rowB =    B[i];
	    
	    for(int j=0; j<numCols; j++)
		rowB[j] = -rowA[j];
	}
    }

    /** Negates matrix A and stores result in B . 
     ** A and B may coinside. 
     * @param  A   The matrix to negate.
     * @param  B   The negated matrix.
     */
    public final static void neg( final int [][] A, final double [][] B ) {

	final int numRows = A   .length; 
	final int numCols = A[0].length; 

	for(int i=0; i<numRows; i++) {
	    final int       [] rowA =    A[i];
	    final double    [] rowB =    B[i];
	    
	    for(int j=0; j<numCols; j++)
		rowB[j] = -rowA[j];
	}
    }

    /** Performs A % b = C. */
    public final static void mod( final int [][] A, final int b, final int [][] C ) {

	final int numRows = C   .length; 
	final int numCols = C[0].length; 

	for(int i=0; i<numRows; i++) {
	  final int [] rowA = A[i];
	  final int [] rowC = C[i];

	  for(int j=0; j<numCols; j++)
	      rowC[j] = rowA[j] % b;
	}
    }

 /** Return the sum of the squares of all entries. */
    public static final double normSqr( final int[][] re) {
	int res = 0;
	int numRows = re.length;
	int numCols = re[0].length;
        for (int i = 0; i < numRows; i++) {
            final int[] rowRe = re[i];
            for (int j = 0; j < numCols; j++)
                res += rowRe[j] * rowRe[j];
        }
        return res;
      }

 /** Return the sum of the squares of all entries. */
    public  static final double normSqr( final double[][] re) {
	double res = 0;
	int numRows = re.length;
	int numCols = re[0].length;
        for (int i = 0; i < numRows; i++) {
            final double[] rowRe = re[i];
            for (int j = 0; j < numCols; j++)
                res += rowRe[j] * rowRe[j];
        }
        return res;
      }

    /**
     * Locates the entry with the greatest absolute value.
     * 
     * @param re real part of matrix
     * @param im imaginary part of matrix
     * @return int[] with indices of biggest entry (<code>{row, col}</code>)
     */
    public static final int[] maxAbs(final double[][] re)
    {
        final int r = re.length;
        final int c = re[0].length;
        int[] pos = new int[] {0, 0};
        double abs, max = 0.;
        for (int i = 0; i < r; i++)
            for (int j = 0; j < c; j++)
            {
                abs = re[i][j] * re[i][j];
                if (abs > max)
                {
                    max = abs;
                    pos[0] = i;
                    pos[1] = j;
                }
            }
        return pos;
    }

    /**
     * Locates the entry with the greatest absolute value.
     * 
     * @param re real part of matrix
     * @param im imaginary part of matrix
     * @return int[] with indices of biggest entry (<code>{row, col}</code>)
     */
    public static final int[] maxAbs(final double[][] re, final double[][] im)
    {
        final int r = re.length;
        final int c = re[0].length;
        int[] pos = new int[] {0, 0};
        double abs, max = 0.;
        for (int i = 0; i < r; i++)
            for (int j = 0; j < c; j++)
            {
                abs = re[i][j] * re[i][j] + im[i][j] * im[i][j];
                if (abs > max)
                {
                    max = abs;
                    pos[0] = i;
                    pos[1] = j;
                }
            }
        return pos;
    }

   
   /** Transpose <i>m</i> and store the result in <i>trans</i> matrix. 
    * @note make checking before use: the sizes must fit to each other 
    * @note if <i>m</i> is a square matrix, <i>m</i>  and <i>trans</i> may coinside
    */
    public static final void transpose( final int [][] m, final int [][] trans ) {

       final int numRows = trans.length;
       if( m == trans ) {
	   for (int i = 0; i < numRows; i++) {
	       for (int j = 0; j < i; j++) {
		   final int tmp = m[i][j];
		   m[i][j] = m[j][i];
		   m[j][i] = tmp;
	       }
	   }
       } else {
	final int numCols = trans[0].length;      
	   for (int i = 0; i < numRows; i++) {
	       int [] rowRe = trans[i];
	       for (int j = 0; j < numCols; j++) {
		   rowRe[j] = m[j] [i];
	       }
	   }
       }
    }
    
   /** Transpose <i>m</i> and store the result in <i>trans</i> matrix. 
    * @note make checking before use: the sizes must fit to each other 
    * @note if <i>m</i> is a square matrix, <i>m</i>  and <i>trans</i> may coinside
    */
    public static final void transpose( final double [][] m, final double [][] trans ) {

       final int numRows = trans.length;
       if( m == trans ) {
	   for (int i = 0; i < numRows; i++) {
	       for (int j = 0; j < i; j++) {
		   final double tmp = m[i][j];
		   m[i][j] = m[j][i];
		   m[j][i] = tmp;
	       }
	   }
       } else {
	   final int numCols = trans[0].length;      
	   for (int i = 0; i < numRows; i++) {
	       double[] rowRe = trans[i];
	       for (int j = 0; j < numCols; j++) {
		   rowRe[j] = m[j] [i];
	       }
	   }
       }
    }
    
    /** Assigns random all entries of m.*/ 
    public static final void random( double [][] m) {
	final int numRows = m.length;
	final int numCols = m[0].length;
	
        for (int i = 0; i < numRows; i++)
            for (int j = 0; j < numCols; j++)
                m[i][j] = 2. * Math.random() - 1.;
    } 

     /** Assigns random all entries of m.
      * all entries are from 0 to range. 
      */ 
    public static final void random( int [][] m, int range ) {
	final int numRows = m.length;
	final int numCols = m[0].length;
	
        for (int i = 0; i < numRows; i++)
            for (int j = 0; j < numCols; j++)
                m[i][j] = (int) Math.random() *range;
    }     

    /** Checks whether a two-dimensional int-array is rectangular.
     * A two-dimension array without entries is considered as
     * rectangular.
     * A two-dimensional int-array with a "null" row is allways 
     * considered to be not rectangular.
     * @param m two-dimension int-array
     * @return true iff the array is rectangular.
     */
    public static boolean isRectangular( final int [][] m ) {
    	if( m.length == 0)
    		return true;
    	final int numCols = m[0].length;
    	final int numRows = m.length;
    	
    	for( int i=1; i<numRows; i++)
    		if( m[i] == null || m[i].length != numCols )
    			return false;
    	return true;
    }
    
    /** Checks whether a two-dimensional double-array is rectangular.
     * A two-dimension array without entries is considered as
     * rectangular.
     * A two-dimensional array with a "null" row is allways 
     * considered to be not rectangular.
     * @param m two-dimension double-array
     * @return true iff the array is rectangular.
     */
    public static boolean isRectangular( final double [][] m ) {
    	if( m.length == 0)
    		return true;
    	final int numCols = m[0].length;
    	final int numRows = m.length;
    	
    	for( int i=1; i<numRows; i++)
    		if( m[i] == null || m[i].length != numCols )
    			return false;
    	return true;
    }
}//end of class
