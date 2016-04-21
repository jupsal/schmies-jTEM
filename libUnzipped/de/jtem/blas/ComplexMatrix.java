/**
This file is part of a jTEM project.
All jTEM projects are licensed under the FreeBSD license 
or 2-clause BSD license (see http://www.opensource.org/licenses/bsd-license.php). 

Copyright (c) 2002-2010, Technische Universität Berlin, jTEM
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

package de.jtem.blas;


import de.jtem.mfc.field.Complex;
import de.jtem.mfc.field.ComplexConstant;
import de.jtem.mfc.field.Field;
import de.jtem.numericalMethods.algebra.linear.Determinant;
import de.jtem.numericalMethods.algebra.linear.Inversion;
import de.jtem.numericalMethods.algebra.linear.MatrixOperations;


/**
 * This class represents general complex matrices.
 *
 * Following the number package's philosophy, operations on complex matrices
 * generally occur three-fold: For instances a,b,c of this class,
 * <table>
 * <tr><td valign="up"><tt>a.plus(b)</tt></td><td>creates and returns
 * a new instance containing the sum <i>a+b</i></td></tr>
 * <tr><td valign="up"><tt>a.assignPlus(b)</tt></td><td>stores the
 * sum <i>a+b</i> in <i>a</i></td></tr>
 * <tr><td valign="up"><tt>c.assignPlus(a,b)</tt></td><td>stores the
 * sum <i>a+b</i> in c</td></tr>
 * </table>
 * <p>
 * Non-commutative operations are <em>right operations</em>:
 * <table>
 * <tr><td valign="up"></td><td><tt>a.assignTimes(b)</tt></td><td>
 * assigns <i>a</i> to <i>a*b</i>,</td></tr>
 * <tr><td valign="up">use</td><td><tt>a.assignTimes(b,a)</tt></td><td>
 * to assign <i>a</i> to <i>b*a</i></td></tr>
 * </table>
 * <p>
 * There are assign/get-methods to interact with instances of the class
 * {@link RealMatrix}.
 * <p>
 * Internally, real and imaginary parts of the matrix are separated
 * and the corresponding arrays can be accessed directly using
 * <tt>re(),im()</tt> to get and <tt>re(double[][])</tt>,
 * <tt>im(double[][])</tt> to assign <em>references (!)</em>.
 *
 * @author	marina
 */
public final class ComplexMatrix extends AbstractMatrix
{
    private static final long serialVersionUID = 1L;

    public double[][] re;
    public double[][] im;

    /** create a zero matrix */
    public ComplexMatrix()    {
        this(0,0);
    }

    /** create a matrix with all entries equal to zero */
    public ComplexMatrix( final int newNumRows,final int newNumCols){
        numRows = newNumRows;
	numCols = newNumCols;
	re=new double[numRows][numCols];
	im=new double[numRows][numCols];
    }

    /** create a matrix and initialize all entries with <i>initValue</i> */
    public ComplexMatrix(final int numrows,final int numcols, final Field.Complex initValue){
	this(numrows,numcols);
	assign( initValue );
    }

    /** create a matrix and initialize all entries with the <i>initRe,initIm</i> */
    public ComplexMatrix(final int numrows,final int numcols, final double initRe,final double initIm){
	this(numrows,numcols);
	assign( initRe, initIm );
    }

    /** create copy of matrix M */
    public ComplexMatrix(final ComplexMatrix M){
	this(M.numRows,M.numCols);
	double[][] MRe = M.re;
	double[][] MIm = M.im;

	for(int i=0; i<numRows; i++) {
	    System.arraycopy( MRe[i], 0, re[i], 0, numCols );
	    System.arraycopy( MIm[i], 0, im[i], 0, numCols );
	}
    }

    /** create copy of matrix M */
    public ComplexMatrix(final RealMatrix M){
	this(M.numRows,M.numCols);

	MatrixOperations.assign(M.re,re);
	MatrixOperations.assign(im,0.);
    }

     /** create copy of matrix M */
    public ComplexMatrix( final IntegerMatrix M ){
	this(M.numRows,M.numCols);
	MatrixOperations.assign( M.re, re );
	MatrixOperations.assign(  im, 0. );
    }

  /** create a quadratic matrix and initialize all entries with zero */
    public ComplexMatrix( final int numrows ){
	this(numrows,numrows);
    }

    /** create a quadratic matrix and initialize all entries with <i>initValue</i> */
    public ComplexMatrix(final int numrows, final Field.Complex initValue){
	this(numrows,numrows,initValue);
    }

    /** create a quadratic matrix and initialize all entries with <i>initRe,initIm</i> */
    public ComplexMatrix( final int numrows, final double initRe, final double initIm ){
	this(numrows,numrows,initRe,initIm);
    }

    /** create a Matrix consisting only of one row or column:
        If <i>covariant</i> is <tt>true</tt>, the matrix has one row
        equal to <i>V</i>, otherwise it has one column equal to <i>V</i>  */
    public ComplexMatrix(final ComplexVector V, final boolean covariant ){
	this( covariant ? 1 : V.size(), covariant ? V.size(): 1 );

	if( covariant ) {
	    System.arraycopy( V.re, 0, re[0], 0, numCols );
	    System.arraycopy( V.im, 0, im[0], 0, numCols );

	} else {
	    double[] VRe = V.re;
	    double[] VIm = V.im;
	    for(int i=0; i<numRows; i++) {
		re[i][0] = VRe[i];
		im[i][0] = VIm[i];
	      }
	}
    }

    /** create a matrix with initialized with <i>ca</i> if possible */
    public ComplexMatrix( final Complex[][] ca ){
      if (ca.length==0) newSize(0,0);
      else {
	checkArray(ca);

	newSize(ca.length,ca[0].length);

	for (int i=0;i<numRows;i++) {
	  double[] Rre = re[i];
	  double[] Rim = im[i];
	  Complex[] Rca = ca[i];
	  for (int j=0;j<numCols;j++) {
	    Rre[j]=Rca[j].re;
	    Rim[j]=Rca[j].im;
	  }
	}
      }
    }

      /** create a matrix with initialized with <i>ca</i> if possible */
    public ComplexMatrix( final Field.Complex[][] ca ){
      if (ca.length==0) newSize(0,0);
      else {
checkArray(ca);

	newSize(ca.length,ca[0].length);

	for (int i=0;i<numRows;i++) {
	  double[] Rre = re[i];
	  double[] Rim = im[i];
	  Field.Complex[] Rca = ca[i];
	  for (int j=0;j<numCols;j++) {
	    Rre[j]=Rca[j].getRe();
	    Rim[j]=Rca[j].getIm();
	  }
	}
      }
    }

    /** create a matrix with initialized with <i>ca</i> if possible */
    public ComplexMatrix( final ComplexConstant[][] ca ){
      if (ca.length==0) newSize(0,0);
      else {
	checkArray(ca);

	newSize(ca.length,ca[0].length);

	for (int i=0;i<numRows;i++) {
	  double[] Rre = re[i];
	  double[] Rim = im[i];
	  ComplexConstant[] Rca = ca[i];
	  for (int j=0;j<numCols;j++) {
	    Rre[j]=Rca[j].getRe();
	    Rim[j]=Rca[j].getIm();
	  }
	}
      }
    }


   /** create a matrix with initialized with <i>ca</i> if possible */
    public ComplexMatrix(final double[][] rm){
      if (rm.length==0) newSize(0,0);
      else {
	checkArray(rm);

	newSize(rm.length,rm[0].length);

	MatrixOperations.assign( rm, re );
	MatrixOperations.assign( im, 0. );
      }
    }

     /** create a matrix with initialized with <i>ca</i> if possible */
    public ComplexMatrix( final int[][] rm ){
      if (rm.length==0) newSize(0,0);
      else {
	checkArray(rm);

	newSize(rm.length,rm[0].length);

	MatrixOperations.assign( rm, re );
	MatrixOperations.assign( im, 0. );
      }
    }

   /** create a matrix with initialized with <i>Re</i> and <i>Im</i> if possible */
    ComplexMatrix( final double[][] Re, final double[][] Im, final boolean performArrayCheck ){
      if (Re.length==0 && Im.length==0) newSize(0,0);
      else {
	  if( performArrayCheck ) {
	      checkArray(Re);
	      checkArray(Im,Re.length,Re[0].length);
	  }

	newSize(Re.length,Re[0].length);

	MatrixOperations.assign( Re, re );
	MatrixOperations.assign( Im, im );
      }
    }

    /** create a matrix with initialized with <i>Re</i> and <i>Im</i> if possible */
    ComplexMatrix( final int[][] Re, final int[][] Im, final boolean performArrayCheck ){
      if (Re.length==0 && Im.length==0) newSize(0,0);
      else {
	  if( performArrayCheck ) {
	      checkArray(Re);
	      checkArray(Im,Re.length,Re[0].length);
	  }

	newSize(Re.length,Re[0].length);

	MatrixOperations.assign( Re, re );
	MatrixOperations.assign( Im, im );
      }
    }

   /** create a matrix with initialized with <i>Re</i> and <i>Im</i> if possible */
    public ComplexMatrix( final double[][] Re, final double[][] Im ){
	this( Re, Im, true );
    }

    /** create a matrix with initialized with <i>Re</i> and <i>Im</i> if possible */
    public ComplexMatrix( final int[][] Re, final int[][] Im ){
	this( Re, Im, true );
    }

    /** create a quadratic matrix with its diagonal initialized to the vector <i>diag</i> */
    public ComplexMatrix( final ComplexVector diag )
    {
	this( diag.size() );

	double[] diagRe = diag.re;
	double[] diagIm = diag.im;

	for(int i=0; i<numRows; i++)
	    {
		re[i][i]=diagRe[i];
		im[i][i]=diagIm[i];
	    }
    }

  void newSize( final int newNumRows, final int newNumCols, final boolean copyValues ) {
      double[][] newRe = null;
      double[][] newIm = null;
      if( newNumCols != numCols ) {
	  newRe = new double[newNumRows][newNumCols];
	  newIm = new double[newNumRows][newNumCols];
	  if( copyValues ) {
	      int minNumCols = newNumCols < numCols ? newNumCols : numCols;
	      for( int i=0; i<newNumRows && i<numRows; i++ ) {
		  System.arraycopy( re[i], 0, newRe[i], 0, minNumCols );
		  System.arraycopy( im[i], 0, newIm[i], 0, minNumCols );
	      }
	  }
      } else if( newNumRows != numRows ) {
	  newRe = new double[newNumRows][];
	  newIm = new double[newNumRows][];
	  int i;
	  for( i=0; i<newNumRows && i<numRows; i++ ) {
	      newRe[i] = re[i];
	      newIm[i] = im[i];
	  }
	  for( ; i<newNumRows; i++ ) {
	      newRe[i] = new double[newNumCols];
	      newIm[i] = new double[newNumCols];
	  }
      } else return;
      re = newRe;
      im = newIm;
      numCols = newNumCols;
      numRows = newNumRows;
    }

    /** create and return new <b>Field.Complex</b> containing entry at <i>rowIndex,colIndex</i> */
    public Complex get( final int rowIndex, final int colIndex ){
	return new Complex(re[rowIndex][colIndex],im[rowIndex][colIndex]);
    }

    /** store entry at <i>rowIndex,colIndex</i> in <i>c</i> */
    public void get(final int rowIndex, final int colIndex, final Complex c){
	c.re = re[rowIndex][colIndex];
	c.im = im[rowIndex][colIndex];
    }

    /** store <i>c</i> in entry <i>rowIndex,colIndex</i> */
    public void set(final int rowIndex, final int colIndex, final Field.Complex c){
	re[rowIndex][colIndex]=c.getRe(); im[rowIndex][colIndex]=c.getIm();

    }

    /** store <i>(aRe,aIm)</i> in entry <i>rowIndex,colIndex</i> */
    public void set(final int rowIndex, final int colIndex, final double aRe, final double aIm ){
	re[rowIndex][colIndex]=aRe; im[rowIndex][colIndex]=aIm;

    }

    /** assign all entries to <i>b</i> */
    public void assign( final Field.Complex b ) {
	assign( b.getRe(), b.getIm() );
    }

    /** assign all entries to <i>b</i> */
    public void assign( final double b ) {
	assign( b, 0 );
    }

    /** assign all entries to <i>(aRe,aIm)</i> */
    public void assign( final double aRe,final double aIm ) {
	MatrixOperations.assign( re, aRe );
	MatrixOperations.assign( im, aIm );
    }

    /** assign matrix to <i>M</i> */
    public void assign( final ComplexMatrix M ) {
	if (this == M) return;
	newSize( M.numRows, M.numCols );
	MatrixOperations.assign(M.re,re);
	MatrixOperations.assign(M.im,im);
    }

    /** assign matrix to <i>M</i> */
    public void assign( final Complex[][] M ) {
        checkArray(M,M.length, M[0].length);
	newSize( M.length, M[0].length );
	for(int i=0; i<numRows; i++)
	  for(int j=0; j<numCols; j++) {
	    re[i][j] = M[i][j].re;
	    im[i][j] = M[i][j].im;
	  }
    }

    /** assign matrix to <i>M</i> */
    public void assign( final ComplexConstant[][] M ) {
        checkArray(M,M.length, M[0].length);
	newSize( M.length, M[0].length );
	for(int i=0; i<numRows; i++)
	  for(int j=0; j<numCols; j++) {
	    re[i][j] = M[i][j].getRe();
	    im[i][j] = M[i][j].getIm();
	  }
    }

    /** assign matrix to <i>M</i> */
    public void assign( final Field.Complex[][] M ) {
        checkArray(M,M.length, M[0].length);
	newSize( M.length, M[0].length );
	for(int i=0; i<numRows; i++)
	  for(int j=0; j<numCols; j++) {
	    re[i][j] = M[i][j].getRe();
	    im[i][j] = M[i][j].getIm();
	  }
    }

    /** assign matrix to <i>Re</i> and <i>Im</i> if possible */
    public void assign(final double[][] Re,final double[][] Im){
      if (Re.length==0 && Im.length==0) newSize(0,0);
      else {
	checkArray(Re);
	checkArray(Im,Re.length,Re[0].length);

	newSize(Re.length,Re[0].length);
	MatrixOperations.assign(Re,re);
	MatrixOperations.assign(Im,im);
      }
    }

    /** assign matrix to <i>M</i> */
    public void assign( final RealMatrix M ) {
	newSize( M.numRows, M.numCols );

	MatrixOperations.assign(M.re,re);
	MatrixOperations.assign(im,0.0);
    }

    /** assign matrix to <i>M</i> */
    public void assign( final double[][] M ) {
	newSize( M.length, M[0].length );

	MatrixOperations.assign(M,re);
	MatrixOperations.assign(im,0.0);
    }

    /** assign matrix to <i>M</i> */
    public void assign( final IntegerMatrix M ) {
	newSize( M.numRows, M.numCols );

	MatrixOperations.assign(M.re,re);
	MatrixOperations.assign(im,0.0);
    }

   /** assign matrix to <i>re</i> and <i>im</i> if possible  */
    public void assign( final RealMatrix re,final RealMatrix im ) {
      assign( re.re, im.re );
    }

    /** assign matrix to <i>re</i> and <i>im</i> if possible  */
    public void assign( final IntegerMatrix re,final IntegerMatrix im ) {
      assign( re.re, im.re );
    }


    /** assign matrix to <i>re</i> and <i>im</i> if possible  */
    public void assign( final int[][] rre,final int[][] iim ) {
	MatrixOperations.assign ( rre, re );
	MatrixOperations.assign ( iim, im );
    }

    /** The following four methods return/assign REFERENCES only ! */

    /** return a REFERENCE to <i>this</i>' real part */
    public double[][] re() {
      return re;
    }

    /** return a REFERENCE to <i>this</i>' imaginary part */
    public double[][] im() {
      return im;
    }

    /** assign <i>this</i>' real part by REFERENCING <i>M</i> */
    public void re( final double[][] M ) {
      checkShape(M);
      re = M;
    }

    /** assign <i>this</i>' imaginary part by REFERENCING <i>M</i> */
    public void im( final double[][] M ) {
      checkShape(M);
      im = M;
    }

    /** assign <i>this</i>' real part to <i>M</i>, leave imaginary part unchanged */
    public void setRe( final double[][] M ) {
      checkShape(M);
      MatrixOperations.assign(M,re);
    }

    /** assign <i>this</i>' imaginary part to <i>M</i>, leave real part unchanged */
    public void setIm( final double[][] M ) {
	checkShape(M);
	MatrixOperations.assign(M,im);
    }

     /** assign <i>this</i>' real part to <i>M</i>, leave imaginary part unchanged */
    public void setRe( final int[][] M ) {
      checkShape(M);
      MatrixOperations.assign(M,re);
    }

    /** assign <i>this</i>' imaginary part to <i>M</i>, leave real part unchanged */
    public void setIm( final int[][] M ) {
	checkShape(M);
	MatrixOperations.assign(M,im);
    }

   /** assign <i>M</i> to <i>this</i>' real part  */
    public void getRe( double[][] M ) {
      checkShape(M);
      MatrixOperations.assign(re,M);
    }

    /** assign <i>M</i> to <i>this</i>' imaginary part  */
    public void getIm( double[][] M ) {
      checkShape(M);
      MatrixOperations.assign(im,M);
    }

    /** create and return new <b>de.jtem.blas.RealMatrix</b> containing <i>this</i>' real part */
    public RealMatrix getRe() {
	return new RealMatrix (re, false );
    }

    /** create and return new <b>de.jtem.blas.RealMatrix</b> containing <i>this</i>' imaginary part */
    public RealMatrix getIm() {
	return new RealMatrix (im, false );
    }

    /** store <i>this</i>' real part in <i>M</i> if possible */
    public void getRe(final RealMatrix M ) {
	M.assign(re);
    }

    /** store <i>this</i>' imaginary part in <i>M</i> if possible */
    public void getIm(final RealMatrix M) {
	M.assign(im);
    }

    /** assign <i>this</i>' real part to <i>M</i> if it has the right shape,
	leave imaginary part unchanged */
    public void setRe( final RealMatrix M ) {
	checkShape( M.re );
	setRe(M.re);
    }

    /** assign <i>this</i>' imaginary part to <i>M</i>, leave real part unchanged */
    public void setIm( final RealMatrix M ) {
	checkShape( M.re );
	setIm(M.re);
    }

    /** create and return new <b>ComplexVector</b> containing row <i>aRowNum</i> */
    public ComplexVector getRow( final int aRowNum ) {
	ComplexVector V = new ComplexVector( numCols );
	getRow( aRowNum, V );
	return V;
    }

    /** store row <i>aRowNum</i> in <i>V</i> */
    public void getRow(final  int aRowNum,final ComplexVector V ) {
	checkRowIndex( aRowNum );
	V.newSize( numCols );
	MatrixOperations.getRow( re, V.re, aRowNum );
	MatrixOperations.getRow( im, V.im, aRowNum );
    }

    /** store <i>V</i> in row <i>aRowNum</i> */
    public void setRow( final int aRowNum, final ComplexVector V ) {
	checkRowIndex( aRowNum );
	if( V.size() != numCols )
	    throw new IllegalArgumentException("size of vector does not match with number of colums");
	MatrixOperations.assignRow( re, V.re, aRowNum );
	MatrixOperations.assignRow( im, V.im, aRowNum );
    }

    /** create and return new <b>ComplexVector</b> containing column <i>aColNum</i> */
    public ComplexVector getCol( final int aColNum ) {
	ComplexVector V = new ComplexVector( numRows );

	getCol( aColNum, V );
	return V;
    }

    /** store column <i>aColNum</i> in <i>V</i> */
    public void getCol( final int aColNum, final ComplexVector V ) {
	checkColIndex( aColNum );
	V.newSize( numRows );

	MatrixOperations.getCol( re, V.re, aColNum );
	MatrixOperations.getCol( im, V.im, aColNum );
    }


    /** store <i>V</i> in column <i>aColNum</i> */
    public void setCol( final int aColNum, final ComplexVector V ) {
	checkColIndex( aColNum );
	if( V.size() != numRows )
	    throw new IllegalArgumentException("size of vector does not match with number of rows");

	MatrixOperations.assignCol( re, V.re, aColNum );
	MatrixOperations.assignCol( im, V.im, aColNum );
    }

    /** create and return a block from <i>(rowIndex,colIndex)</i> to the lower right corner */
    public ComplexMatrix getBlock( final int rowIndex, final int colIndex ){
	return getBlock( rowIndex, colIndex, numRows-rowIndex, numCols-colIndex );
    }

    /** create and return a block of size <i>(numOfRows,numOfCols)</i> beginning at <i>(rowIndex,colIndex)</i> */
    public ComplexMatrix getBlock( final int rowIndex, final int colIndex,final int numOfRows,final int numOfCols ) {
	ComplexMatrix M = new ComplexMatrix( numOfRows, numOfCols );
	getBlock( rowIndex, colIndex, numOfRows, numOfCols, M );
	return M;
    }

    /** store a block from <i>(rowIndex,colIndex)</i> to the lower right corner in <i>M</i> */
    public ComplexMatrix getBlock( final int rowIndex, final int colIndex, final ComplexMatrix M ){
	return getBlock( rowIndex, colIndex, numRows-rowIndex, numCols-colIndex, M );
    }

    /** store a block of size <i>(numOfRows,numOfCols)</i> beginning at <i>(rowIndex,colIndex)</i> in <i>M</i> */
    public ComplexMatrix getBlock( final int rowIndex, final int colIndex, final int numOfRows, final int numOfCols, ComplexMatrix M ) {

	checkRowIndex( rowIndex );
	checkRowIndex( rowIndex + numOfRows - 1 );
	checkColIndex( colIndex );
	checkColIndex( colIndex + numOfCols - 1 );

	M.newSize( numOfRows, numOfCols );

	double [][] MRe = M.re;
	double [][] MIm = M.im;

	for( int i=rowIndex, I=0; I<numOfRows; i++, I++ ) {

	    System.arraycopy( re[i], colIndex, MRe[I], 0, numOfCols );
	    System.arraycopy( im[i], colIndex, MIm[I], 0, numOfCols );
	}

	return M;
    }

    /** store <i>M</i> as a block in <i>this</i> matrix, beginning at <i>(rowIndex,colIndex)</i> */
    public void setBlock( final int rowIndex, final int colIndex, final ComplexMatrix M ) {
	setBlock( rowIndex, colIndex, M, M.numRows, M.numCols );
    }

    /** store a block from M of size <i>(numOfRows,numOfCols)</i> in <i>this</i> matrix, beginning at <i>(rowIndex,colIndex)</i> */
    public void setBlock( final ComplexMatrix M, final int numOfRows, final int numOfCols ) {
	setBlock( 0, 0, M, numOfRows, numOfCols );
    }

    /** store a block of size <i>(numOfRows,numOfCols)</i> beginning at <i>M</i>'s upper left corner
        in <i>this</i> matrix, beginning at <i>(rowIndex,colIndex)</i> */
    public void setBlock( final int rowIndex, final int colIndex, final ComplexMatrix M, final int numOfRows, final int numOfCols ) {

	checkRowIndex( rowIndex );
	checkRowIndex( rowIndex + numOfRows - 1 );
	checkColIndex( colIndex );
	checkColIndex( colIndex + numOfCols - 1 );

	M.checkRowIndex( numOfRows - 1 );
	M.checkColIndex( numOfCols - 1 );

	double [][] MRe = M.re;
	double [][] MIm = M.im;

	for( int i=rowIndex, I=0; I<numOfRows; i++, I++ ) {

	    System.arraycopy( MRe[I], 0, re[i], colIndex, numOfCols );
	    System.arraycopy( MIm[I], 0, im[i], colIndex, numOfCols );
	}
    }


    /** assign this matrix to the adjoint of M */
    public void assignAdjoint( final ComplexMatrix M, final int ithRow, final int jthCol ) {
	M.checkRowIndex( ithRow );
	M.checkColIndex( jthCol );

	newSize( M.numRows-1, M.numCols-1 );

	for( int i=0; i<ithRow; i++ ) {
	    System.arraycopy( M.re[i], 0,          re[i], 0,                  jthCol     );
	    System.arraycopy( M.re[i], jthCol + 1, re[i], jthCol, M.numCols - jthCol - 1 );

	    System.arraycopy( M.im[i], 0,          im[i], 0,                  jthCol     );
	    System.arraycopy( M.im[i], jthCol + 1, im[i], jthCol, M.numCols - jthCol - 1 );
	}

	for( int i=ithRow; i<numRows; i++ ) {
	    System.arraycopy( M.re[i+1], 0,          re[i], 0,                  jthCol     );
	    System.arraycopy( M.re[i+1], jthCol + 1, re[i], jthCol, M.numCols - jthCol - 1 );

	    System.arraycopy( M.im[i+1], 0,          im[i], 0,                  jthCol     );
	    System.arraycopy( M.im[i+1], jthCol + 1, im[i], jthCol, M.numCols - jthCol - 1 );
	}

	if( (ithRow + jthCol) % 2 == 1 )
	    assignTimes( -1 );
    }

    public ComplexMatrix getAdjoint(  final int ithRow, final int jthCol ) {
	ComplexMatrix adjoint = new ComplexMatrix( numRows-1, numCols-1 );
	adjoint.assignAdjoint( this, ithRow, jthCol );
	return adjoint;
    }

    /** stores <i>this</i>' diagonal in <i>diag</i> */
    public void getDiagonal( final ComplexVector diag ) {
      int numDiag = Math.min(numRows,numCols);

      if (diag.size() != numDiag) throw new IllegalArgumentException("vector has wrong length");

      MatrixOperations.getDiagonal( re,diag.re);
      MatrixOperations.getDiagonal( im,diag.im);
    }

    /** create and return a <b>ComplexVector</b> containing <i>this</i>' diagonal */
    public ComplexVector getDiagonal() {
      ComplexVector diag = new ComplexVector(Math.min(numRows,numCols));
      getDiagonal(diag);
      return diag;
    }

    /** stores <i>diag</i> in <i>this</i>' diagonal while other entries remain unchanged */
    public void setDiagonal( final ComplexVector diag ) {
      int numDiag = Math.min(numRows,numCols);

      if (diag.size() != numDiag) throw new IllegalArgumentException("vector has wrong length");

      MatrixOperations.assignDiagonal( re,diag.re);
      MatrixOperations.assignDiagonal( im,diag.im);
    }


    /** stores <i>diag</i> in <i>this</i>' diagonal while other entries remain unchanged */
    public void setDiagonal( final Complex[] diag ) {
      int numDiag = Math.min(numRows,numCols);

      if (diag.length != numDiag) throw new IllegalArgumentException("vector has wrong length");

      for(int i=0;i<numDiag;i++){
	  re[i][i]=diag[i].re;
	  im[i][i]=diag[i].im;
      }
    }

     /** stores <i>diag</i> in <i>this</i>' diagonal while other entries remain unchanged */
    public void setDiagonal( final ComplexConstant[] diag ) {
      int numDiag = Math.min(numRows,numCols);

      if (diag.length != numDiag) throw new IllegalArgumentException("vector has wrong length");

      for(int i=0;i<numDiag;i++){
	  re[i][i]=diag[i].getRe();
	  im[i][i]=diag[i].getIm();
      }
    }

     /** stores <i>diag</i> in <i>this</i>' diagonal while other entries remain unchanged */
    public void setDiagonal( final Field.Complex[] diag ) {
      int numDiag = Math.min(numRows,numCols);

      if (diag.length != numDiag) throw new IllegalArgumentException("vector has wrong length");

      for(int i=0;i<numDiag;i++){
	  re[i][i]=diag[i].getRe();
	  im[i][i]=diag[i].getIm();
      }
    }

   /** assign diagonal to @param diagRe + i @param diagIm while other entries remain unchanged */
    public void setDiagonal( final double diagRe, final double diagIm  ) {
	int numDiag = Math.min(numRows,numCols);
	MatrixOperations.assignDiagonal(  re, diagRe);
	MatrixOperations.assignDiagonal(  im, diagIm);
    }

  /** @deprecated  */
    public void assignDiagonal( final double diagRe, final double diagIm  ) {
	int numDiag = Math.min(numRows,numCols);
	MatrixOperations.assignDiagonal(  re, diagRe);
	MatrixOperations.assignDiagonal(  im, diagIm);
    }

   /** assign diagonal to @param diag while other entries remain unchanged */
    public void setDiagonal( final Field.Complex diag ) {
	setDiagonal( diag.getRe(), diag.getIm() );
    }

    /** create and return an array <b>Field.Complex[][]</b> containing <i>this</i>' entries */
    public Complex[][] toArray(){
	Complex[][] res=new Complex[numRows][numCols];
	for(int i=0; i<numRows; i++)
	    for( int j=0; j<numCols; j++ )
		{
		    res[i][j]=new Complex(re[i][j],im[i][j]);
		}
	return res;
    }

    /** store <i>this</i>' entries in <i>M</i>*/
    public void toArray(Complex[][] M){
        checkShape(M);
	for (int i=0; i<numRows; i++) {
	  double[] Rre = re[i];
	  double[] Rim = im[i];
	  for (int j=0; j<numCols; j++) M[i][j].assign(Rre[j],Rim[j]);
	}
    }

    /** create and return a new matrix containing <i>this</i>' entries */
    public ComplexMatrix copy(){
	return new ComplexMatrix(this);
    }

    /** return a clone */
    public Object clone() { return copy(); }


    /** create and return an identity matrix of size <i>numrows*numrows</i> */
    public static ComplexMatrix id(int numrows){
	ComplexMatrix ID = new ComplexMatrix(numrows);
	for(int i=0; i<numrows; i++) ID.re[i][i]=1;
	return ID;
    }

    /** assign <i>this</i> to an identity matrix of size <i>numrows*numrows</i> */
    public void assignId(final int numrows){
        newSize(numrows,numrows);
	assignZero();
	for(int i=0; i<numrows; i++) re[i][i]=1;
    }

    /** assign <i>this</i> to an identity matrix of size <i>numRows*numRows</i> */
    public void assignId(){
	assignId(numRows);
    }

    /** add <i>b</i> and <i>this</i> and create and return a matrix containing the result */
    public ComplexMatrix plus(final ComplexMatrix b){
        checkShape(b);

	double[][] sumRe=new double[numRows][numCols];
	double[][] sumIm=new double[numRows][numCols];
	MatrixOperations.plus(re,b.re,sumRe);
	MatrixOperations.plus(im,b.im,sumIm);

	return new ComplexMatrix(sumRe,sumIm);
    }

    /** subtract <i>b</i> from <i>this</i> and create and return a matrix containing the result */
    public ComplexMatrix minus(final ComplexMatrix b){
        checkShape(b);

	double[][] sumRe=new double[numRows][numCols];
	double[][] sumIm=new double[numRows][numCols];
	MatrixOperations.minus(re,b.re,sumRe);
	MatrixOperations.minus(im,b.im,sumIm);

	return new ComplexMatrix(sumRe,sumIm);
    }

    /** add <i>b</i> and <i>this</i> and create and return a matrix containing the result */
    public ComplexMatrix plus(final RealMatrix b){
        checkShape(b);

	double[][] sumRe=new double[numRows][numCols];
	MatrixOperations.plus(re,b.re,sumRe);
	return new ComplexMatrix(sumRe,im);
    }

    /** subtract <i>b</i> from <i>this</i> and create and return a matrix containing the result */
    public ComplexMatrix minus(final RealMatrix b){
        checkShape(b);
	double[][] sumRe=new double[numRows][numCols];
	MatrixOperations.minus(re,b.re,sumRe);
	return new ComplexMatrix(sumRe,im);
    }

    /** add <i>b</i> and <i>this</i> and create and return a matrix containing the result */
    public ComplexMatrix plus(final IntegerMatrix b){
        checkShape(b);

	double[][] sumRe=new double[numRows][numCols];
	MatrixOperations.plus(re,b.re,sumRe);
	return new ComplexMatrix(sumRe,im);
    }

    /** subtract <i>b</i> from <i>this</i> and create and return a matrix containing the result */
    public ComplexMatrix minus(final IntegerMatrix b){
        checkShape(b);
	double[][] sumRe=new double[numRows][numCols];
	MatrixOperations.minus(re,b.re,sumRe);
	return new ComplexMatrix(sumRe,im);
    }

    /** add <i>a</i> and <i>b</i> and store result in <i>this</i> matrix */
    public void assignPlus(final ComplexMatrix a, final ComplexMatrix b){
        a.checkShape(b);

	resize( a.numRows, a.numCols );
	MatrixOperations.plus(a.re,b.re,re);
	MatrixOperations.plus(a.im,b.im,im);
    }

     /** add <i>a</i> and <i>b</i> and store result in <i>this</i> matrix */
    public void assignPlus(final RealMatrix a, final RealMatrix b){
        a.checkShape(b);
	resize( a.numRows, a.numCols );
	MatrixOperations.plus(a.re,b.re,re);
	MatrixOperations.assign( im, 0. );
    }

    /** add <i>a</i> and <i>b</i> and store result in <i>this</i> matrix */
    public void assignPlus(final IntegerMatrix a, final IntegerMatrix b){
        a.checkShape(b);
	resize( a.numRows, a.numCols );
	int[][] t = new int [numRows][numCols];
	MatrixOperations.plus(a.re,b.re,t);
	MatrixOperations.assign( t, re );
	MatrixOperations.assign( im, 0. );
    }

   /** add <i>a</i> and <i>b</i> and store result in <i>this</i> matrix */
    public void assignPlus(final ComplexMatrix a, final RealMatrix b){
        a.checkShape(b);
	resize( a.numRows, a.numCols );
	MatrixOperations.plus(a.re,b.re,re);
	MatrixOperations.assign( a.im, im );
    }

    /** add <i>a</i> and <i>b</i> and store result in <i>this</i> matrix */
    public void assignPlus(final IntegerMatrix a, final RealMatrix b){
        a.checkShape(b);
	resize( a.numRows, a.numCols );
	MatrixOperations.plus(a.re,b.re,re);
	MatrixOperations.assign( im, 0. );
    }

     /** add <i>a</i> and <i>b</i> and store result in <i>this</i> matrix */
    public void assignPlus( final RealMatrix a, final IntegerMatrix b){
	assignPlus(b,a);
    }


  /** add <i>a</i> and <i>b</i> and store result in <i>this</i> matrix */
    public void assignPlus( final RealMatrix a, final ComplexMatrix b ){
	assignPlus(b, a);
    }

    /** add <i>a</i> and <i>b</i> and store result in <i>this</i> matrix */
    public void assignPlus( final ComplexMatrix a, final IntegerMatrix b ){
        a.checkShape(b);
	resize( a.numRows, a.numCols );
	MatrixOperations.plus(a.re,b.re,re);
	MatrixOperations.assign( a.im, im );
    }

    /** add <i>a</i> and <i>b</i> and store result in <i>this</i> matrix */
    public void assignPlus( final IntegerMatrix a, final ComplexMatrix b ){
	assignPlus(b, a);
    }

    /** subtract <i>b</i> from <i>a</i> and store result in <i>this</i> matrix */
    public void assignMinus( final ComplexMatrix a, final ComplexMatrix b ){
        a.checkShape(b);

	resize( a.numRows, a.numCols );
	MatrixOperations.minus(a.re,b.re,re);
	MatrixOperations.minus(a.im,b.im,im);
    }

    /** subtract <i>b</i> from <i>a</i> and store result in <i>this</i> matrix */
    public void assignMinus( final RealMatrix a, final RealMatrix b ){
        a.checkShape(b);
	resize( a.numRows, a.numCols );
	MatrixOperations.minus(a.re,b.re,re);
	MatrixOperations.assign( im, 0. );
    }

    /** subtract <i>b</i> from <i>a</i> and store result in <i>this</i> matrix */
    public void assignMinus(final IntegerMatrix a, final IntegerMatrix b){
        a.checkShape(b);
	resize( a.numRows, a.numCols );
	int[][] t = new int [numRows][numCols];
	MatrixOperations.minus(a.re,b.re,t);
	MatrixOperations.assign( t, re );
	MatrixOperations.assign( im, 0. );
    }

   /** subtract <i>b</i> from <i>a</i> and store result in <i>this</i> matrix */
    public void assignMinus(final ComplexMatrix a, final RealMatrix b){
        a.checkShape(b);

	resize( a.numRows, a.numCols );
	MatrixOperations.minus(a.re,b.re,re);
	MatrixOperations.assign( a.im, im );

    }

    /** subtract <i>b</i> from <i>a</i> and store result in <i>this</i> matrix */
    public void assignMinus(final IntegerMatrix a, final RealMatrix b){
        a.checkShape(b);

	resize( a.numRows, a.numCols );
	MatrixOperations.minus(a.re,b.re,re);
	MatrixOperations.assign( im, 0. );
    }

    /** subtract <i>b</i> from <i>a</i> and store result in <i>this</i> matrix */
    public void assignMinus(final RealMatrix a, final IntegerMatrix b){
        a.checkShape(b);

	resize( a.numRows, a.numCols );
	MatrixOperations.minus(a.re,b.re,re);
	MatrixOperations.assign( im, 0. );
    }

   /** subtract <i>b</i> from <i>a</i> and store result in <i>this</i> matrix */
    public void assignMinus(final ComplexMatrix a, final IntegerMatrix b){
        a.checkShape(b);

	resize( a.numRows, a.numCols );
	MatrixOperations.minus(a.re,b.re,re);
	MatrixOperations.assign( a.im, im );
    }

    /** subtract <i>b</i> from <i>a</i> and store result in <i>this</i> matrix */
    public void assignMinus(final IntegerMatrix a, final ComplexMatrix b){
        a.checkShape(b);

	resize( a.numRows, a.numCols );
	MatrixOperations.minus(a.re,b.re,re);
	MatrixOperations.neg( b.im, im );
    }

    /** subtract <i>b</i> from <i>a</i> and store result in <i>this</i> matrix */
    public void assignMinus(final RealMatrix a, final ComplexMatrix b){
        a.checkShape(b);

	resize( a.numRows, a.numCols );
	MatrixOperations.minus(a.re,b.re,re);
	MatrixOperations.neg( b.im, im );
    }

    /** add <i>b</i> and <i>this</i> and store the result in <i>this</i> matrix */
    public void assignPlus( final IntegerMatrix b){
        checkShape(b);
	MatrixOperations.plus(re,b.re,re);
    }

     /** add <i>b</i> and <i>this</i> and store the result in <i>this</i> matrix */
    public void assignPlus( final RealMatrix b){
        checkShape(b);
	MatrixOperations.plus(re,b.re,re);
    }

   /** add <i>b</i> and <i>this</i> and store the result in <i>this</i> matrix */
    public void assignPlus(final ComplexMatrix b){
        checkShape(b);
	MatrixOperations.plus(re,b.re,re);
	MatrixOperations.plus(im,b.im,im);
    }

    /** subtract <i>b</i> from <i>this</i> and store the result in <i>this</i> matrix */
    public void assignMinus(final ComplexMatrix b){
        checkShape(b);

	MatrixOperations.minus(re,b.re,re);
	MatrixOperations.minus(im,b.im,im);
    }

    /** subtract <i>b</i> from <i>this</i> and store the result in <i>this</i> matrix */
    public void assignMinus(final RealMatrix b){
        checkShape(b);
	MatrixOperations.minus(re,b.re,re);
    }

    /** subtract <i>b</i> from <i>this</i> and store the result in <i>this</i> matrix */
    public void assignMinus(final IntegerMatrix b){
        checkShape(b);
	MatrixOperations.minus(re,b.re,re);
    }

    /* **************************************************************************************** */
    /* **************       matrix-scalar multiplication      ********************************* */
    /* **************************************************************************************** */

    /** multiply <i>this</i> and <i>a</i> and create and return a matrix containing the result */
    public ComplexMatrix times( final Field.Complex a ){
	ComplexMatrix result = new ComplexMatrix( numRows, numCols );
	result.assignTimes( this, a );
	return result;
    }

    /** multiply <i>this</i> and <i>a</i> and create and return a matrix containing the result */
    public ComplexMatrix times( final double a ){
	ComplexMatrix result = new ComplexMatrix( numRows, numCols );
	result.assignTimes( this, a );
	return result;
    }

    /** multiply <i>this</i> and <i>a</i> and store the result in <i>this</i> matrix */
    public void assignTimes( final Field.Complex a ){
	assignTimes( this, a );
    }

    /** multiply <i>this</i> and <i>a</i> and store the result in <i>this</i> matrix */
    public void assignTimes( final double a ){
	assignTimes( this, a );
    }

    /** multiply <i>A</i> and <i>b</i> and store the result in <i>this</i> matrix */
    public void assignTimes( final ComplexMatrix A,final Field.Complex b  ){
	newSize( A.numRows, A.numCols );
	MatrixOperations.times( A.re, A.im, b.getRe(), b.getIm(), re, im );
    }
    /** multiply <i>A</i> and <i>b</i> and store the result in <i>this</i> matrix */
    public void assignTimes( final Field.Complex b , final ComplexMatrix A){
	assignTimes( A,b  );
    }

    /** multiply <i>A</i> and <i>b</i> and store the result in <i>this</i> matrix */
    public void assignTimes( final ComplexMatrix A,  final double b ) {
	newSize( A.numRows, A.numCols );
	MatrixOperations.times( A.re, b, re );
	MatrixOperations.times( A.im, b, im );
    }

    /** multiply <i>A</i> and <i>b</i> and store the result in <i>this</i> matrix */
    public void assignTimes( final double b ,final ComplexMatrix A) {
	assignTimes( A, b );
    }

    /** multiply <i>A</i> and <i>b</i> and store the result in <i>this</i> matrix */
    public void assignTimes( final RealMatrix A, final Field.Complex b ) {

	newSize( A.numRows, A.numCols );
	MatrixOperations.times( A.re, b.getRe(), re );
	MatrixOperations.times( A.re, b.getIm(), im );
    }

    /** multiply <i>A</i> and <i>b</i> and store the result in <i>this</i> matrix */
    public void assignTimes( final IntegerMatrix A, final Field.Complex b ) {

	newSize( A.numRows, A.numCols );

	MatrixOperations.times( A.re, b.getRe(), re );
	MatrixOperations.times( A.re, b.getIm(), im );
    }

    /** multiply <i>A</i> and <i>b</i> and store the result in <i>this</i> matrix */
    public void assignTimes( final RealMatrix A, final double b ) {

	newSize( A.numRows, A.numCols );

	MatrixOperations.times( A.re, b, re );
	MatrixOperations.assign( im, 0 );
    }

    /** multiply <i>A</i> and <i>b</i> and store the result in <i>this</i> matrix */
    public void assignTimes( final IntegerMatrix A, final  double b ) {

	newSize( A.numRows, A.numCols );

	MatrixOperations.times( A.re, b, re );
	MatrixOperations.assign( im, 0 );
    }

    /* **************************************************************************************** */
    /* **************       matrix-matrix multiplication      ********************************* */
    /* **************************************************************************************** */

    /** multiply <i>this</i> and <i>A</i> and create and return a matrix containing the result */
    public ComplexMatrix times(final ComplexMatrix A){
	ComplexMatrix result = new ComplexMatrix( numRows, A.numCols );
	result.assignTimes( this, A );
	return result;
    }

    /** multiply <i>this</i> and <i>A</i> and create and return a matrix containing the result */
    public ComplexMatrix times(final RealMatrix A){
	ComplexMatrix result = new ComplexMatrix( numRows, A.numCols );
	result.assignTimes( this, A );
	return result;
    }

    /** multiply <i>this</i> and <i>A</i> and create and return a matrix containing the result */
    public ComplexMatrix times(final IntegerMatrix A){
	ComplexMatrix result = new ComplexMatrix( numRows, A.numCols );
	result.assignTimes( this, A );
	return result;
    }

    /** multiply <i>this</i> and <i>A</i> and store the result in <i>this</i> matrix */
    public void assignTimes(final ComplexMatrix A){
	assignTimes( this, A );
    }

    /** multiply <i>this</i> and <i>A</i> and store the result in <i>this</i> matrix */
    public void assignTimes(final RealMatrix A){
	assignTimes( this, A );
    }

    /** multiply <i>this</i> and <i>A</i> and store the result in <i>this</i> matrix */
    public void assignTimes(final IntegerMatrix A) {
	assignTimes( this, A );
    }

    /** multiply <i>A</i> and <i>B</i> and store the result in <i>this</i> matrix */
    public void assignTimes( final ComplexMatrix A, final ComplexMatrix B) {
	if(A.numCols!=B.numRows)
	    throw new IllegalArgumentException("matrices have wrong sizes");

	final ComplexMatrix tmpA = ( A == this )? new ComplexMatrix( this ) : A;
	final ComplexMatrix tmpB = ( B == this )? new ComplexMatrix( this ) : B;

	newSize( tmpA.numRows, tmpB.numCols );

	MatrixOperations.times( tmpA.re, tmpA.im, tmpB.re, tmpB.im, re, im);
    }

      /** multiply <i>A</i> and <i>B</i> and store the result in <i>this</i> matrix */
    public void assignTimes( final RealMatrix A, final RealMatrix B) {
	if(A.numCols!=B.numRows)
	    throw new IllegalArgumentException("matrices have wrong sizes");
	newSize( A.numRows, B.numCols );

	MatrixOperations.times( A.re, B.re, re );
	MatrixOperations.assign( im , 0. );
    }

       /** multiply <i>A</i> and <i>B</i> and store the result in <i>this</i> matrix */
    public void assignTimes( final IntegerMatrix A, final IntegerMatrix B) {
	if(A.numCols!=B.numRows)
	    throw new IllegalArgumentException("matrices have wrong sizes");
	newSize( A.numRows, B.numCols );

	MatrixOperations.times( A.re, B.re, re );
	MatrixOperations.assign( im , 0. );
    }

 /** multiply <i>A</i> and <i>B</i> and store the result in <i>this</i> matrix */
    public void assignTimes(final ComplexMatrix A, final RealMatrix B) {

	if(A.numCols!=B.numRows)
	    throw new IllegalArgumentException("matrices have wrong sizes");
	newSize( A.numRows, B.numCols );

	MatrixOperations.times( A.re, B.re, re );
	MatrixOperations.times( A.im, B.re, im );
    }

    /** multiply <i>A</i> and <i>B</i> and store the result in <i>this</i> matrix */
    public void assignTimes(final ComplexMatrix A, final IntegerMatrix B) {

	if(A.numCols!=B.numRows)
	    throw new IllegalArgumentException("matrices have wrong sizes");
	newSize( A.numRows, B.numCols );

	MatrixOperations.times( A.re, B.re, re );
	MatrixOperations.times( A.im, B.re, im );
    }

    /** multiply <i>A</i> and <i>B</i> and store the result in <i>this</i> matrix */
    public void assignTimes( final RealMatrix A, final ComplexMatrix B) {

	if(A.numCols!=B.numRows)
	    throw new IllegalArgumentException("matrices have wrong sizes");
	newSize( A.numRows, B.numCols );

	MatrixOperations.times( A.re, B.re, re );
	MatrixOperations.times( A.re, B.im, im );
    }

    /** multiply <i>A</i> and <i>B</i> and store the result in <i>this</i> matrix */
    public void assignTimes( final IntegerMatrix A, final ComplexMatrix B) {

	if(A.numCols!=B.numRows)
	    throw new IllegalArgumentException("matrices have wrong sizes");
	newSize( A.numRows, B.numCols );

	MatrixOperations.times( A.re, B.re, re );
	MatrixOperations.times( A.re, B.im, im );
    }

    /** multiply <i>A</i> and <i>B</i> and store the result in <i>this</i> matrix */
    public void assignTimes(final RealMatrix A, final IntegerMatrix B) {

	if(A.numCols!=B.numRows)
	    throw new IllegalArgumentException("matrices have wrong sizes");
	newSize( A.numRows, B.numCols );

	MatrixOperations.times( A.re, B.re, re );
	MatrixOperations.assign( im, 0 );
    }

    /** multiply <i>A</i> and <i>B</i> and store the result in <i>this</i> matrix */
    public void assignTimes(final IntegerMatrix A, final RealMatrix B) {

	if(A.numCols!=B.numRows)
	    throw new IllegalArgumentException("matrices have wrong sizes");
	newSize( A.numRows, B.numCols );

	MatrixOperations.times( A.re, B.re, re );
	MatrixOperations.assign( im, 0 );
    }

    /* ***** end of matrix-matrix multiplication ***** */


    /** multiply <i>this</i> and <i>V</i> and create and return a vector containing the result */
    public ComplexVector times( final ComplexVector V ){
	ComplexVector res=new ComplexVector( numRows );
	res.assignTimes( this, V );
	return res;
    }

     public ComplexVector times( final RealVector V ){
	ComplexVector res=new ComplexVector( numRows );
	res.assignTimes( this, V );
	return res;
    }

      public ComplexVector times( final IntegerVector V ){
	ComplexVector res=new ComplexVector( numRows );
	res.assignTimes( this, V );
	return res;
    }

   /** divide all entries by <i>b</i> and create and return a matrix containing the result */
    public ComplexMatrix divide( final Field.Complex b ){
	ComplexMatrix result = new ComplexMatrix();
	result.assignDivide( this, b );
	return result;
    }

    /** divide all entries by <i>r</i> and create and return a matrix containing the result */
    public ComplexMatrix divide(final double r){
	ComplexMatrix result = new ComplexMatrix( numRows, numCols );
	result.assignDivide( this, r );
	return result;
    }

    /** divide all entries of <i>a</i> by <i>b</i> and store the result in <i>this</i> matrix */
    public void assignDivide( final ComplexMatrix M, final Field.Complex b){
	newSize( M.numRows, M.numCols );
	MatrixOperations.divide( M.re, M.im, b.getRe(), b.getIm(), re, im );
    }

    /** divide all entries of <i>a</i> by <i>r</i> and store the result in <i>this</i> matrix */
    public void assignDivide(final ComplexMatrix a, final double r ){
	newSize( a.numRows, a.numCols );
	MatrixOperations.divide( a.re, r, re );
	MatrixOperations.divide( a.im, r, im );
    }

    /** divide all entries of <i>a</i> by <i>r</i> and store the result in <i>this</i> matrix */
    public void assignDivide(final RealMatrix a, final double r ){
	newSize( a.numRows, a.numCols );

	MatrixOperations.divide( a.re , r, re );
	MatrixOperations.assign (im, 0. );
    }

    /** divide all entries of <i>a</i> by <i>b</i> and store the result in <i>this</i> matrix */
    public void assignDivide( final RealMatrix M, final Field.Complex b){
	newSize( M.numRows, M.numCols );
	MatrixOperations.divide( M.re, M.re, b.getRe(), b.getIm(), re, im );
    }

    /** divide all entries of <i>a</i> by <i>r</i> and store the result in <i>this</i> matrix */
    public void assignDivide(final IntegerMatrix a, final double r ){
	newSize( a.numRows, a.numCols );

	MatrixOperations.divide( a.re , r, re );
	MatrixOperations.assign (im, 0. );
    }

    /** divide all entries of <i>a</i> by <i>b</i> and store the result in <i>this</i> matrix */
    public void assignDivide( final IntegerMatrix M, final Field.Complex b){
	newSize( M.numRows, M.numCols );
	MatrixOperations.divide( M.re, b.getRe(), re );
	MatrixOperations.divide( M.re, b.getIm(), im );
    }

    /** divide all entries of <i>a</i> by <i>b</i> and store the result in <i>this</i> matrix */
    public void assignDivide( final Field.Complex b, final ComplexMatrix M ){
	newSize( M.numRows, M.numCols );
	MatrixOperations.divide(  b.getRe(), M.re, re );
	MatrixOperations.divide(  b.getIm(), M.im, im );
    }

    /** divide all entries of <i>a</i> by <i>r</i> and store the result in <i>this</i> matrix */
    public void assignDivide( final double r, final ComplexMatrix a ){
	newSize( a.numRows, a.numCols );
	MatrixOperations.divide(  r, a.re, a.im, re, im );
    }

    /** divide all entries of <i>a</i> by <i>r</i> and store the result in <i>this</i> matrix */
    public void assignDivide(final double r,final RealMatrix a  ){
	newSize( a.numRows, a.numCols );
	MatrixOperations.divide(  r, a.re , re );
	MatrixOperations.assign (im, 0. );
    }

    /** divide all entries of <i>a</i> by <i>b</i> and store the result in <i>this</i> matrix */
    public void assignDivide( final Field.Complex b,final RealMatrix M ){
	newSize( M.numRows, M.numCols );
	MatrixOperations.divide(  b.getRe(), M.re, re );
	MatrixOperations.divide(  b.getIm(), M.re, im );
    }

    /** divide all entries of <i>a</i> by <i>r</i> and store the result in <i>this</i> matrix */
    public void assignDivide( final double r, final IntegerMatrix a ){
	newSize( a.numRows, a.numCols );

	MatrixOperations.divide(  r,a.re , re );
	MatrixOperations.assign (im, 0. );
    }

    /** divide all entries of <i>a</i> by <i>b</i> and store the result in <i>this</i> matrix */
    public void assignDivide(  final Field.Complex b, final IntegerMatrix M){
	newSize( M.numRows, M.numCols );
	MatrixOperations.divide(  b.getRe(), M.re, re );
	MatrixOperations.divide(  b.getIm(), M.re, im );
    }

    /** divide all entries of <i>a</i> by <i>r</i> */
    public void assignDivide( final double r ){
	MatrixOperations.divide( re , r, re );
	MatrixOperations.divide( im , r, im );
    }

      /** divide all entries of <i>a</i> by <i>b</i> */
    public void assignDivide( final Field.Complex b ){
	assignDivide(this,b);
    }

    /** conjugate <i>M</i> and create and return a matrix containing the result */
    public static ComplexMatrix conjugate( final ComplexMatrix M ) {
	ComplexMatrix res = new ComplexMatrix( M.numRows, M.numCols );
	res.assignConjugate( M );
	return res;
    }

    /** conjugate <i>this</i> matrix and create and return a matrix containing the result */
    public ComplexMatrix conjugate() {
	ComplexMatrix res = new ComplexMatrix( numRows, numCols );
	res.assignConjugate( this );
	return res;
    }

    /** conjugate <i>M</i> and store the result in <i>this</i> matrix */
    public void assignConjugate( final ComplexMatrix M ) {
	resize( M.numRows, M.numCols );

	MatrixOperations.assign(M.re, re);
	MatrixOperations.neg(   M.im, im);
    }

    /** conjugate <i>this</i> matrix */
   public void assignConjugate() {
       MatrixOperations.neg(im,im);
   }

    /** return <tt>true</tt> if <i>this</i> is symmetric */
    public boolean isSymmetric() {
	if( !isSquared() )
	    return false;
	for( int i=0; i<numRows; i++ ) {
	    double[] Rre = re[i];
	    double[] Rim = im[i];
	    for( int j=1; j<i; j++ ) {
	        double dummy1 = Rre[j]-re[j][i];
	        double dummy2 = Rim[j]-im[j][i];
		if( dummy1*dummy1+dummy2*dummy2 > EPSILON )
		    return false;
	    }
	}
	return true;
    }

    /** return the sum of the squares of all entries */
    public double normSqr(){
	return (MatrixOperations.normSqr(re)+MatrixOperations.normSqr(im));
    }

    /** return the sum of the squares of all entries of <i>M</i> */
    public static double normSqr( final ComplexMatrix M){
	return M.normSqr();
    }

   /** return <tt>true</tt> if all entries are zero */
    public static boolean isZero(final ComplexMatrix M) {
      return M.isZero();
    }

    /** return <tt>true</tt> if <i>this</i> equals <i>m</i> */
    public boolean equals(final ComplexMatrix m){
	return equals( m, EPSILON );
    }

     /** return <tt>true</tt> if <i>this</i> equals <i>m</i> */
    public boolean equals(final ComplexMatrix m, final double eps ){
	if( m== null || m.numRows!=numRows || m.numCols!=numCols )
          return false;
	for(int i=0; i<numRows; i++) {
	    double[] Rre1 = re[i];
	    double[] Rim1 = im[i];
	    double[] Rre2 = m.re[i];
	    double[] Rim2 = m.im[i];
	    for(int j=0; j<numCols; j++) {
              if( Math.abs( Rre1[j]-Rre2[j] ) > eps ||
                  Math.abs( Rim1[j]-Rim2[j] ) > eps )
                  return false;
              }
	}
	return true;
    }

   /** return <tt>true</tt> if <i>this</i> equals <i>o</i> */
    public boolean equals(final Object o ){
      if (o == null)
        return false;
      try {
        return equals( (ComplexMatrix) o);
      }
      catch (ClassCastException ex) {
        return false;
      }
    }

    /** creates a <b>String</b> containing this' matrix' entries */
    public String toString(){
	StringBuffer sb=new StringBuffer(300);
	sb.append("(");
	for(int i=0; i<numRows; i++)
	    {
		sb.append('(');
		if(numCols>0)
		    {
			sb.append(re[i][0]);
			if(im[i][0]>=0.) sb.append('+');
			sb.append(im[i][0]);
			sb.append('i');
			for(int j=1; j<numCols; j++)
			    {
				sb.append(", ");
				sb.append(re[i][j]);
				if(im[i][j]>=0.) sb.append('+');
				sb.append(im[i][j]);
				sb.append('i');
			    }
		    }
		sb.append(')');
		if(i<numRows-1) sb.append(",\n");
	    }
	sb.append(')');
	return sb.toString();
    }

    /** transpose <i>M</i> and create and return a matrix containing the result */
    public static ComplexMatrix transpose( final ComplexMatrix M ) {
        ComplexMatrix res = new ComplexMatrix( M.numCols, M.numRows );
	res.assignTranspose( M );
	return res;
    }

    /** transpose <i>this</i> matrix create and return a matrix containing the result */
    public ComplexMatrix transpose() {
        ComplexMatrix res = new ComplexMatrix( numCols, numRows );
	res.assignTranspose( this );
	return res;
    }

    /** transpose <i>M</i> and store the result in <i>this</i> matrix */
    public void assignTranspose( final ComplexMatrix M ) {

	ComplexMatrix tmpM = ( M.isSquared() || M != this ) ? M : (new ComplexMatrix( this ));

	newSize( tmpM.numCols, tmpM.numRows );

	MatrixOperations.transpose( tmpM.re, re );
	MatrixOperations.transpose( tmpM.im, im );
    }

    /** transpose <i>this</i> matrix */
    public void assignTranspose() {
	assignTranspose(this);
    }


    /**
     * Computes the determinant.
     * <p />
     * This method temporary uses one additional matrix and two vectors.
     * <p />
     * @return	The determinant.
     * @see de.jtem.numericalMethods.algebra.linear.Determinant#compute(double[][],double[][]) Determinant
     */
    public Complex determinant(){
        if (numRows != numCols)
            throw new IllegalArgumentException("matrix must be square");
        double[] det = Determinant.compute(re, im);
        return new Complex(det[0], det[1]);
    }

    /**
     * Computes and returns the inverse of this matrix.
     * <p />
     * This method temporary uses two additional matrices and one vector.
     * <p />
     * @return	The inverse.
     * @see de.jtem.numericalMethods.algebra.linear.Inversion#compute(double[][],double[][],double[][],double[][]) Inversion
     */
    public ComplexMatrix invert() {
        ComplexMatrix inv = new ComplexMatrix();
        inv.assignInvert(this);
        return inv;
    }

    /**
     * Computes and returns the inverse of a matrix.
     * <p />
     * This method temporary uses two additional matrices and one vector.
     * <p />
     * @param	M   The matrix to invert.
     * @return	The inverse.
     * @see de.jtem.numericalMethods.algebra.linear.Inversion#compute(double[][],double[][],double[][],double[][]) Inversion
     */
    public static ComplexMatrix invert( ComplexMatrix M ) {
        ComplexMatrix inv = new ComplexMatrix();
        inv.assignInvert(M);
        return inv;
    }

    /**
     * Computes the inverse of a matrix and stores the result in this matrix.
     * <p />
     * This method temporary uses one additional matrix and one vector.
     * <p />
     * @param	M   The matrix to invert.
     * @see de.jtem.numericalMethods.algebra.linear.Inversion#compute(double[][],double[][],double[][],double[][]) Inversion
     */
    public void assignInvert(ComplexMatrix M){
        if (M.numRows != M.numCols)
            throw new IllegalArgumentException("matrix must be square");
        resize(M.numRows);
        if (!Inversion.compute(M.re, M.im, re, im))
            throw new IllegalArgumentException("matrix was singular");
    }

    /**
     * Inverts this matrix.
     * <p />
     * This method temporary uses one additional matrix and one vector.
     * <p />
     * @see de.jtem.numericalMethods.algebra.linear.Inversion#compute(double[][],double[][],double[][],double[][]) Inversion
     */
    public void assignInvert() {
        ComplexMatrix inv = new ComplexMatrix();
        inv.assignInvert(this);
        assign(inv);
    }

    /** assign all entries to zero */
    public void assignZero() {
	MatrixOperations.assign( re, 0 );
	MatrixOperations.assign( im, 0 );
    }

    /** assign all entries randomly */
    public void assignRandom(){
	MatrixOperations.random( re );
	MatrixOperations.random( im );
    }

    /** returns round of <code>A</code> */
    public static ComplexMatrix round( final ComplexMatrix A ) {
	ComplexMatrix result = new ComplexMatrix( A.numRows, A.numCols );
	MatrixOperations.round( A.re, result.re );
	MatrixOperations.round( A.im, result.im );
	return result;
    }

   /** returns round of<em>this</em> */
    public ComplexMatrix round() {
	ComplexMatrix result = new ComplexMatrix( numRows, numCols );
	MatrixOperations.round( re, result.re );
	MatrixOperations.round( im, result.im );
	return result;
    }

    /** assign <em>this</em> to round of <code>A</code> */
    public void assignRound( final ComplexMatrix A ) {
	newSize( A.numRows, A.numCols );
	MatrixOperations.round( A.re, re );
	MatrixOperations.round( A.im, im );
    }

     /** assign <em>this</em> to round of <code>A</code> */
    public void assignRound( final RealMatrix A ) {
	newSize( A.numRows, A.numCols );
	MatrixOperations.round( A.re, re );
	MatrixOperations.assign(   im, 0. );
    }

    /**
     * Rounds <i>this</i>.
     */
    public void assignRound() {
	MatrixOperations.round( re, re );
	MatrixOperations.round( im, im );
    }

     /** returns round of <code>A</code> */
    public static ComplexMatrix floor( final ComplexMatrix A ) {
	ComplexMatrix result = new ComplexMatrix( A.numRows, A.numCols );
	MatrixOperations.floor( A.re, result.re );
	MatrixOperations.floor( A.im, result.im );
	return result;
    }

   /** returns round of<em>this</em> */
    public ComplexMatrix floor() {
	ComplexMatrix result = new ComplexMatrix( numRows, numCols );
	MatrixOperations.floor( re, result.re );
	MatrixOperations.floor( im, result.im );
	return result;
    }

    /** assign <em>this</em> to round of <code>A</code> */
    public void assignFloor( final ComplexMatrix A ) {
	newSize( A.numRows, A.numCols );
	MatrixOperations.floor( A.re, re );
	MatrixOperations.floor( A.im, im );
    }

     /** assign <em>this</em> to round of <code>A</code> */
    public void assignFloor( final RealMatrix A ) {
	newSize( A.numRows, A.numCols );
	MatrixOperations.floor( A.re, re );
	MatrixOperations.assign(   im, 0. );
    }

    /**
     * Rounds <i>this</i>.
     */
    public void assignFloor() {
	MatrixOperations.floor( re, re );
	MatrixOperations.floor( im, im );
    }


   /**
     * Returns negative of <code>A</code>.
     */
    public static ComplexMatrix neg( final ComplexMatrix A ) {
	ComplexMatrix result = new ComplexMatrix( A.numRows, A.numCols );
	MatrixOperations.neg( A.re, result.re );
	MatrixOperations.neg( A.im, result.im );
	return result;
    }

    /**
     * Returns negative of <code>this</code>.
     */
    public ComplexMatrix neg() {
	ComplexMatrix result = new ComplexMatrix( numRows, numCols );
	MatrixOperations.neg( re, result.re );
	MatrixOperations.neg( im, result.im );
	return result;
    }

    /**
     * Assigns <code>this</code> with negative of <code>A</code>.
     */
    public void assignNeg( final ComplexMatrix A ) {
	newSize( A.numRows, A.numCols );
	MatrixOperations.neg( A.re, re );
	MatrixOperations.neg( A.im, im );

    }


    /** assign neg of <code>v</code> */
    public final void assignNeg( final RealMatrix v ) {
	newSize( v.numRows, v.numCols );
	MatrixOperations.neg( v.re, re );
	MatrixOperations.assign( im, 0.0 );
    }

     /** assign neg of <code>v</code> */
    public final void assignNeg( final IntegerMatrix v ) {
	newSize( v.numRows ,v. numCols );
	MatrixOperations.neg( v.re, re );
	MatrixOperations.assign( im , 0.0 );
    }

    /**
     * Negates <code>this</code>.
     */
    public void assignNeg() {
	MatrixOperations.neg( re, re );
	MatrixOperations.neg( im, im );
    }



}



