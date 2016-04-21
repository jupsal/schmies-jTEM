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

import de.jtem.mfc.field.Field;
import de.jtem.numericalMethods.algebra.linear.Determinant;
import de.jtem.numericalMethods.algebra.linear.Inversion;
import de.jtem.numericalMethods.algebra.linear.MatrixOperations;



/**
 * This class represents general double matrices.
 * 
 * Following the number package's philosophy, operations on double matrices 
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
 * The array containing the matrix's content can be accessed directly using
 * <tt>re()</tt> to get and <tt>re(double[][])</tt>, to set <em>references (!)</em>.
 *
 * @author	marina
 */
public final class RealMatrix extends AbstractMatrix
{
    private static final long serialVersionUID = 1L;

    public double[][] re;

    /** create a zero matrix */
    public RealMatrix()    {
       this(0,0);
    }

    /** create a matrix with all entries equal to zero */
    public RealMatrix(final int newNumRows, final int newNumCols){
        numRows = newNumRows;
	numCols = newNumCols;
	re=new double[numRows][numCols];
    }

    /** create a matrix and initialize all entries with <i>initValue</i> */
    public RealMatrix(final int newNumRows, final int newNumCols, final double initValue){
	this(newNumRows, newNumCols);
	assign( initValue );	
    }  

    /** create copy of matrix A */
    public RealMatrix(final RealMatrix A){
	this(A.numRows,A.numCols);
	double[][] ARe = A.re;
	
	for(int i=0; i<numRows; i++) {
	    System.arraycopy( ARe[i], 0, re[i], 0, numCols );
	}
    }
    
    /** create copy of matrix A */
    public RealMatrix(final IntegerMatrix A){
	this( A.re, true );
    }

    /** create a quadratic matrix with all entries equal to zero */
    public RealMatrix( final int numrows){
	this(numrows,numrows);
    }

    /** create a quadratic matrix with all entries equal to the initValue */
    public RealMatrix(final int numrows, final double initValue){
	this(numrows,numrows,initValue);   
    }  

    /** create a Matrix consisting only of one row or column:
        If <i>covariant</i> is <tt>true</tt>, the matrix has one row 
        equal to <i>V</i>, otherwise it has one column equal to <i>V</i>  */
    public RealMatrix( final RealVector V, final boolean covariant ){
	this( covariant ? 1 : V.size(), covariant ? V.size(): 1 );

	if( covariant ) {
	    System.arraycopy( V.re, 0, re[0], 0, numCols );
	    
	} else {
	    double[] VRe = V.re;

	    for(int i=0; i<numRows; i++)
		re[i][0] = VRe[i];
	}
    }

    /** create a quadratic matrix with its diagonal initialized to the vector <i>diag</i> */
    public RealMatrix( final RealVector diag ){
	this( diag.size() );
	for(int i=0; i<numRows; i++) re[i][i]=diag.re[i];
    }

    /** create a matrix with initialized with <i>m</i> if possible */
    RealMatrix ( final double[][] m, boolean performArrayCheck ) {
	assign( m, performArrayCheck );
    }

    /** create a matrix with initialized with <i>m</i> if possible */
    RealMatrix ( final int[][] m, boolean performArrayCheck ) {
	assign( m, performArrayCheck );
    }

    public RealMatrix ( double[][] Re ){
	this( Re, true );
    }

    public RealMatrix ( int[][] In ){
	this( In, true );
    }

    void newSize( final int newNumRows, final int newNumCols,final boolean copyValues ) {

      double[][] newRe = null;	    

      if( newNumCols != numCols ) {
	  newRe = new double[newNumRows][newNumCols];
	  if( copyValues ) {
	      int minNumCols = newNumCols < numCols ? newNumCols : numCols;
	      for( int i=0; i<newNumRows && i<numRows; i++ ) {
		  System.arraycopy( re[i], 0, newRe[i], 0, minNumCols );
	      }
	  }
      } else if( newNumRows != numRows ) {
	  newRe = new double[newNumRows][];
	  int i;
	  for( i=0; i<newNumRows && i<numRows; i++ ) {
	      newRe[i] = re[i];
	  }
	  for( ; i<newNumRows; i++ ) {
	      newRe[i] = new double[newNumCols];
	  }
      } else return;
      re = newRe;
      numCols = newNumCols;
      numRows = newNumRows;
    }

    /** return entry at <i>rowIndex,colIndex</i> */
    public double get(final int rowIndex, final int colIndex)    {
	return re[rowIndex][colIndex];
    }

    /** store <i>c</i> in entry <i>rowIndex,colIndex</i> */
    public void set(final int rowIndex,final  int colIndex,final  double c){
	re[rowIndex][colIndex]=c;
    }

    /** assign all entries to <i>c</i> */
    public void assign(final double c ) {
	MatrixOperations.assign( re, c );
    }

    void assign(final double [][] m, final boolean performArrayCheck ) {
	if (m.length==0) {
	    newSize(0,0);
	} else {
	    if( performArrayCheck )
		checkArray(m);
	    newSize(m.length,m[0].length);
	    MatrixOperations.assign(m,re);
	}
    }

    void assign( final int [][] m, final boolean performArrayCheck ) {
	if (m.length==0) {
	    newSize(0,0);
	} else {
	    if( performArrayCheck )
		checkArray(m);
	    newSize(m.length,m[0].length);
	    MatrixOperations.assign(m,re);
	}
    }

    public void assign( final double[][] m )  {
	assign( m, true );
    }

    /** assign matrix to <i>M</i> */
    public void assign(final RealMatrix M ) {
	assign( M.re, false );
    }

    public void assign( final int [][] m )  {
	assign( m, true );
    }

    /** assign matrix to <i>M</i> */
    public void assign( final IntegerMatrix M ) {
	assign( M.re, false );
    }

    /** assign <i>M</i> to <i>this</i> matrix */
    public void get(final ComplexMatrix M )     {
	M.setRe(this);
    }

    /** The following two methods return/assign REFERENCES only ! */
 
    /** return the REFERENCE to <i>this</i>' array */
    public double[][] re()     {
      return re;
    }

    /** assign <i>this</i>' array by REFERENCING <i>M</i> */
    public void re(final double[][] M) {
      checkShape(M);
      re = M;
    }

    /** create and return new <b>de.jtem.blas.RealVector</b> containing row <i>aRowNum</i> */
    public RealVector getRow( final int aRowNum ) {
	return getRow( aRowNum, null );
    }

    /** store row <i>aRowNum</i> in <i>V</i> */
    public RealVector getRow( final int aRowNum, RealVector V ) {
	checkRowIndex( aRowNum );
	if( V==null)                             // why do you do that ?!?
	    V = new RealVector ( numCols );
	else
	    V.newSize(numCols);
	System.arraycopy( re[aRowNum], 0, V.re, 0, numCols );
	return V;
    }

    /** store <i>V</i> in row <i>aRowNum</i> */
    public void setRow( final int aRowNum, RealVector V ) {
	checkRowIndex( aRowNum );
	if( V.size() != numCols )
	    throw new IllegalArgumentException("size of vector does not match with number of colums");
	System.arraycopy( V.re, 0, re[aRowNum], 0, numCols );
    }

    /** store <i>V</i> in row <i>aRowNum</i> */
    public void setRow( final int aRowNum, double[] V ) {
	checkRowIndex( aRowNum );
	if( V.length != numCols )
	    throw new IllegalArgumentException("size of vector does not match with number of colums");
	System.arraycopy( V, 0, re[aRowNum], 0, numCols );
    }

 
    /** store <i>V</i> in row <i>aRowNum</i> */
    public void setRow( final int aRowNum, IntegerVector V ) {
	checkRowIndex( aRowNum );
	if( V.size() != numCols )
	    throw new IllegalArgumentException("size of vector does not match with number of colums");
	MatrixOperations.assignRow( re, V.re , aRowNum );
    }

    /** store <i>V</i> in row <i>aRowNum</i> */
    public void setRow( final int aRowNum, int[] V ) {
	checkRowIndex( aRowNum );
	if( V.length != numCols )
	    throw new IllegalArgumentException("size of vector does not match with number of colums");
	MatrixOperations.assignRow( re, V , aRowNum );
    }

   /** create and return new <b>de.jtem.blas.RealVector</b> containing column <i>aColNum</i> */
    public RealVector getCol( final int aColNum ) {
	return getCol( aColNum, null );
    }

    /** store column <i>aColNum</i> in <i>V</i> */
    public RealVector getCol( final int aColNum, RealVector V ) {
	checkColIndex( aColNum );
	if( V==null)
	    V = new RealVector ( numRows );
	else
	    V.newSize(numRows);
	double [] VRe = V.re;
	for( int i=0; i<numRows; i++ ) {
	    VRe[i] = re[i][aColNum];
	}
	return V;
    }

    /** store <i>V</i> in column <i>aColNum</i> */
    public void setCol( final int aColNum,final RealVector V ) {
	checkColIndex( aColNum );
	if( V.size() != numRows )
	    throw new IllegalArgumentException("size of vector does not match with number of rows");
	MatrixOperations.assignCol( this.re, V.re, aColNum );
    }

     /** store <i>V</i> in column <i>aColNum</i> */
    public void setCol( final int aColNum, final double[] V ) {
	checkColIndex( aColNum );
	if( V.length != numRows )
	    throw new IllegalArgumentException("size of vector does not match with number of rows");
	MatrixOperations.assignCol( this.re, V, aColNum );
    }
    
    /** store <i>V</i> in column <i>aColNum</i> */
    public void setCol( final int aColNum,final IntegerVector V ) {
	checkColIndex( aColNum );
	if( V.size() != numRows )
	    throw new IllegalArgumentException("size of vector does not match with number of rows");
	MatrixOperations.assignCol( this.re, V.re, aColNum );
    }

     /** store <i>V</i> in column <i>aColNum</i> */
    public void setCol( final int aColNum, final int[] V ) {
	checkColIndex( aColNum );
	if( V.length != numRows )
	    throw new IllegalArgumentException("size of vector does not match with number of rows");
	MatrixOperations.assignCol( this.re, V, aColNum );
    }

   /** create and return a block from <i>(rowIndex,colIndex)</i> to the lower right corner */
    public RealMatrix getBlock( final int rowIndex,final int colIndex ) {
	RealMatrix M = new RealMatrix ( numRows-rowIndex, numCols-colIndex );
	return getBlock(rowIndex,colIndex,numRows-rowIndex,numCols-colIndex,M);
    }

    /** create and return a block of size <i>(numOfRows,numOfCols)</i> beginning at <i>(rowIndex,colIndex)</i> */
    public RealMatrix getBlock( final int rowIndex,final int colIndex,final int numOfRows,final int numOfCols ) {
	RealMatrix M = new RealMatrix ( numOfRows, numOfCols );
	return getBlock(rowIndex,colIndex,numOfRows,numOfCols,M);
    }

    /** store a block of size <i>(numOfRows,numOfCols)</i> beginning at <i>(rowIndex,colIndex)</i> in <i>M</i> */
    public RealMatrix getBlock(final int rowIndex,final int colIndex, RealMatrix M ) {
        return getBlock(rowIndex,colIndex,numRows-rowIndex,numCols-colIndex,M);
    }

    /** store a block of size <i>(numOfRows,numOfCols)</i> beginning at <i>(rowIndex,colIndex)</i> in <i>M</i> */
    public RealMatrix getBlock( final int rowIndex, final int colIndex,final int numOfRows,final int numOfCols, RealMatrix M ) {

	checkRowIndex( rowIndex );
	checkRowIndex( rowIndex + numOfRows - 1 );
	checkColIndex( colIndex );
	checkColIndex( colIndex + numOfCols - 1 );

	M.newSize( numOfRows, numOfCols );

	double [][] MRe = M.re;

	for( int i=rowIndex, j=0; j<numOfRows; i++, j++ ) 
	    System.arraycopy( re[i], colIndex, MRe[j], 0, numOfCols );
	
	return M;
    }
 
    /** store <i>M</i> as a block in </i>this</i> matrix beginning at <i>(numOfRows,numOfCols)</i>  */
    public void setBlock( final int rowIndex,final  int colIndex,final  RealMatrix M ) {
	setBlock( rowIndex, colIndex, M, M.numRows, M.numCols );
    }

    /** store a block of size <i>(numOfRows,numOfCols)</i> from<i>M</i> as a block in </i>this</i> matrix beginning at <i>(0,0)</i>  */
    public void setBlock( final RealMatrix M, final int numRows, final int numCols ) {
	setBlock( 0, 0, M, numRows, numCols );
    }

    /** store a block of size <i>(numOfRows,numOfCols)</i> from <i>M</i> in </i>this</i> matrix beginning at <i>(rowIndex,colIndex)</i>  */
    public void setBlock( final int rowIndex, final int colIndex, final RealMatrix M, final int numOfRows, final int numOfCols ) {
	
	checkRowIndex( rowIndex );
	checkRowIndex( rowIndex + numOfRows - 1 );
	checkColIndex( colIndex );
	checkColIndex( colIndex + numOfCols - 1 );

	M.checkRowIndex( numOfRows - 1 );
	M.checkColIndex( numOfCols - 1 );
	
	double [][] MRe = M.re;

	for( int i=0, j=rowIndex; i<numOfRows; i++, j++ ) 
	    System.arraycopy( MRe[i], 0, re[j], colIndex, numOfCols );

	
    }

    /** assign this matrix to the adjoint of M */
    public void assignAdjoint( final RealMatrix M, final int ithRow, final int jthCol ) {
	M.checkRowIndex( ithRow );
	M.checkColIndex( jthCol );
	newSize( M.numRows-1, M.numCols-1 );
	for( int i=0; i<ithRow; i++ ) {
	    System.arraycopy( M.re[i], 0,          re[i], 0,                  jthCol     );
	    System.arraycopy( M.re[i], jthCol + 1, re[i], jthCol, M.numCols - jthCol - 1 );
	}
	for( int i=ithRow; i<numRows; i++ ) {
	    System.arraycopy( M.re[i+1], 0,          re[i], 0,                  jthCol     );
	    System.arraycopy( M.re[i+1], jthCol + 1, re[i], jthCol, M.numCols - jthCol - 1 );
	}
	if( (ithRow + jthCol) % 2 == 1 )
	    assignTimes( -1 );
    }

    public RealMatrix getAdjoint( final int ithRow, final int jthCol ) {
	RealMatrix adjoint = new RealMatrix ( numRows-1, numCols-1 );
	adjoint.assignAdjoint( this, ithRow, jthCol );
	return adjoint;
    }

    /** stores <i>this</i>' diagonal in <i>diag</i> */
    public void getDiagonal( final RealVector diag ) {
      int numDiag = Math.min(numRows,numCols);
      if (diag.size() != numDiag) throw new IllegalArgumentException("vector has wrong length");
      double[] diagRe = diag.re;
      MatrixOperations.getDiagonal( re, diag.re );
    }

    /** create and return a <b>de.jtem.blas.RealVector</b> containing <i>this</i>' diagonal */
    public RealVector getDiagonal() {
      RealVector diag = new RealVector (Math.min(numRows,numCols));
      getDiagonal(diag);
      return diag;
    }

    /** stores <i>diag</i> in <i>this</i>' diagonal while other entries remain unchanged */
    public void setDiagonal( final RealVector diag ) {
      int numDiag = Math.min(numRows,numCols);
      if (diag.size() != numDiag) throw new IllegalArgumentException("vector has wrong length");
      MatrixOperations.assignDiagonal( re, diag.re );
    }

    /** stores <i>diag</i> in <i>this</i>' diagonal while other entries remain unchanged */
    public void setDiagonal( final double[] diag ) {
      int numDiag = Math.min(numRows,numCols);
      if (diag.length != numDiag) throw new IllegalArgumentException("vector has wrong length");
      MatrixOperations.assignDiagonal( re, diag );
    }

     /** stores <i>diag</i> in <i>this</i>' diagonal while other entries remain unchanged */
    public void setDiagonal( final IntegerVector diag ) {
      int numDiag = Math.min(numRows,numCols);
      if (diag.size() != numDiag) throw new IllegalArgumentException("vector has wrong length");
      MatrixOperations.assignDiagonal( re, diag.re );
    }

    /** stores <i>diag</i> in <i>this</i>' diagonal while other entries remain unchanged */
    public void setDiagonal( final int[] diag ) {
      int numDiag = Math.min(numRows,numCols);
      if (diag.length != numDiag) throw new IllegalArgumentException("vector has wrong length");
      MatrixOperations.assignDiagonal( re, diag );
    }

   /** assign diagonal to @param diag while other entries remain unchanged */
    public void setDiagonal( final double diag ) {
	int numDiag = Math.min(numRows,numCols);
	MatrixOperations.assignDiagonal(re, diag );	
    }

    /** create and return an array <b>Field.Complex[][]</b> containing <i>this</i>' entries */
    public double[][] toArray()    {
	double[][] res=new double[numRows][numCols];
	for(int i=0; i<numRows; i++)
	    System.arraycopy(re[i],0,res[i],0,numCols);
	return res;
    }

    /** store <i>this</i>' entries in <i>M</i>*/
    public void toArray(final double[][] M){
        checkShape(M);
	for(int i=0; i<numRows; i++)
	  System.arraycopy(re[i], 0, M[i], 0, numCols);
    }
    
    /** create and return a new matrix containing <i>this</i>' entries */
    public RealMatrix copy(){
	return new RealMatrix (this);
    }

    /** return a clone */
    public Object clone() { return copy(); }


    /** create and return an identity matrix of size <i>numrows*numrows</i> */
    public static RealMatrix id(final int numrows){
	RealMatrix ID = new RealMatrix (numrows);
	ID.setDiagonal( 1.0 );
	return ID;
    }

    /** assign <i>this</i> to an identity matrix of size <i>numrows*numrows</i> */
    public void assignId(final int numrows){
        newSize(numrows,numrows);
	assignZero();
	setDiagonal( 1 );
    }

    /** assign <i>this</i> to an identity matrix of size <i>numRows*numRows</i> */
    public void assignId(){
	assignId(numRows);
    }



    /** add <i>b</i> and <i>this</i> and create and return a matrix containing the result */
    public RealMatrix plus( final RealMatrix b ){
	checkShape(b);
        RealMatrix result = new RealMatrix( numRows, numCols );	
	MatrixOperations.plus( re, b.re, result.re );
	return result;
    }

    /** subtract <i>b</i> from <i>this</i> and create and return a matrix containing the result */
    public RealMatrix minus( final RealMatrix b ){
	checkShape(b);
	RealMatrix result = new RealMatrix( numRows, numCols );
	MatrixOperations.minus( re, b.re, result.re );
	return result;
    }

    /** add <i>b</i> and <i>this</i> and create and return a matrix containing the result */
    public RealMatrix plus( final IntegerMatrix b ){
	checkShape(b);
        RealMatrix result = new RealMatrix( numRows, numCols );	
	MatrixOperations.plus( re, b.re, result.re );
	return result;
    }

    /** subtract <i>b</i> from <i>this</i> and create and return a matrix containing the result */
    public RealMatrix minus( final IntegerMatrix b ){
	checkShape(b);
	RealMatrix result = new RealMatrix( numRows, numCols );
	MatrixOperations.minus( re, b.re, result.re );
	return result;
    }

    /** add <i>b</i> and <i>this</i> and create and return a matrix containing the result */
    public ComplexMatrix plus( final ComplexMatrix b ){
	checkShape(b);
        ComplexMatrix result = new ComplexMatrix( numRows, numCols );	
	MatrixOperations.plus( re, b.re, result.re );
	return result;
    }

    /** subtract <i>b</i> from <i>this</i> and create and return a matrix containing the result */
    public ComplexMatrix minus( final ComplexMatrix b ){
	checkShape(b);
	ComplexMatrix result = new ComplexMatrix( numRows, numCols );
	MatrixOperations.minus( re, b.re, result.re );
	return result;
    }

    /** add <i>a</i> and <i>b</i> and store result in <i>this</i> matrix */
    public void assignPlus(final RealMatrix a, final RealMatrix b){
        a.checkShape(b);
	newSize(a);
	MatrixOperations.plus(a.re,b.re,re);
    }

    /** add <i>a</i> and <i>b</i> and store result in <i>this</i> matrix */
    public void assignPlus(final IntegerMatrix a, final IntegerMatrix b){
        a.checkShape(b);
	newSize(a);
	int [][] t = new int[numRows][numCols];
	MatrixOperations.plus(a.re,b.re,t);
	assign( t );
    }

    /** add <i>a</i> and <i>b</i> and store result in <i>this</i> matrix */
    public void assignPlus(final RealMatrix a, final IntegerMatrix b){
        a.checkShape(b);
	newSize(a);
	MatrixOperations.plus(a.re,b.re,re);
    }

    /** add <i>a</i> and <i>b</i> and store result in <i>this</i> matrix */
    public void assignPlus(final IntegerMatrix a, final RealMatrix b){
	assignPlus( b, a);
    }

    /** subtract <i>b</i> from <i>a</i> and store result in <i>this</i> matrix */
    public void assignMinus( final RealMatrix a, final RealMatrix b ){
      a.checkShape(b);
      newSize(a);
      MatrixOperations.minus( a.re, b.re, re);
    }

   /** subtract <i>b</i> from <i>a</i> and store result in <i>this</i> matrix */
    public void assignMinus(final IntegerMatrix a, final IntegerMatrix b){
        a.checkShape(b);
	newSize(a);
	int [][] t = new int[numRows][numCols];
	MatrixOperations.minus(a.re,b.re,t);
	assign( t );
    }

   /** subtract <i>b</i> from <i>a</i> and store result in <i>this</i> matrix */
    public void assignMinus( final RealMatrix a,final IntegerMatrix b ){
	a.checkShape(b);
	newSize(a);
	MatrixOperations.minus( a.re, b.re, re);
    }

    /** subtract <i>b</i> from <i>a</i> and store result in <i>this</i> matrix */
    public void assignMinus( final IntegerMatrix a,final RealMatrix b ){
        a.checkShape(b);
	newSize(a);
	MatrixOperations.minus(a.re,b.re,re);
    }

    /** add <i>b</i> and <i>this</i> and store the result in <i>this</i> matrix */
    public void assignPlus( final RealMatrix b ){
        checkShape(b);
	MatrixOperations.plus(re,b.re,re);
    }

    /** subtract <i>b</i> from <i>this</i> and store the result in <i>this</i> matrix */
    public void assignMinus( final RealMatrix b ){
	checkShape(b);
	MatrixOperations.minus(re,b.re,re);
    }

    /** add <i>b</i> and <i>this</i> and store the result in <i>this</i> matrix */
    public void assignPlus( final IntegerMatrix b ){
        checkShape(b);
	MatrixOperations.plus(re,b.re,re);
    }

    /** subtract <i>b</i> from <i>this</i> and store the result in <i>this</i> matrix */
    public void assignMinus( final IntegerMatrix b ){
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
    public RealMatrix times( final double a ){
	RealMatrix result = new RealMatrix ( numRows, numCols );
	result.assignTimes( this, a );
	return result;
    }

   /** multiply <i>this</i> and <i>a</i> and store the result in <i>this</i> matrix */
    public void assignTimes( final double a ){
	MatrixOperations.times( re, a, re );
    }
    
    /** multiply <i>A</i> and <i>b</i> and store the result in <i>this</i> matrix */
    public void assignTimes(final RealMatrix A, final double b){
	newSize( A.numRows, A.numCols );
	MatrixOperations.times( A.re, b, re );
    }

    /** multiply <i>A</i> and <i>b</i> and store the result in <i>this</i> matrix */
    public void assignTimes(final double b,final RealMatrix A ){
	assignTimes( A, b);
    }

    /** multiply <i>A</i> and <i>b</i> and store the result in <i>this</i> matrix */
    public void assignTimes(final double b,final IntegerMatrix A ){
	assignTimes( A, b);
    }

    /** multiply <i>A</i> and <i>b</i> and store the result in <i>this</i> matrix */
    public void assignTimes(final IntegerMatrix A, final double b){
	newSize( A.numRows, A.numCols );
	MatrixOperations.times( A.re, b, re );
    }

    /* **************************************************************************************** */
    /* **************       matrix-matrix multiplication      ********************************* */
    /* **************************************************************************************** */


    /** multiply <i>this</i> and <i>A</i> and create and return a matrix containing the result */
    public ComplexMatrix times( final ComplexMatrix A ){
	ComplexMatrix result = new ComplexMatrix( numRows, A.numCols );
	result.assignTimes( this, A );
	return result;
    }

    /** multiply <i>this</i> and <i>A</i> and create and return a matrix containing the result */
    public RealMatrix times( final RealMatrix A ){
	RealMatrix result = new RealMatrix ( numRows, A.numCols );
	result.assignTimes( this, A );
	return result;
    }


    /** multiply <i>this</i> and <i>A</i> and create and return a matrix containing the result */
    public RealMatrix times( final IntegerMatrix A ){
	RealMatrix result = new RealMatrix ( numRows, A.numCols );
	result.assignTimes( this, A );
	return result;
    }

    /** multiply <i>this</i> and <i>A</i> and store the result in <i>this</i> matrix */
    public void assignTimes(final RealMatrix A){
      assignTimes( this, A );
    }
     
    /** multiply <i>this</i> and <i>A</i> and store the result in <i>this</i> matrix */
    public void assignTimes(final IntegerMatrix A){ 
      assignTimes( this, A );
    }
    
    /** multiply <i>A</i> and <i>B</i> and store the result in <i>this</i> matrix */
    public void assignTimes(final RealMatrix A, final RealMatrix B){
	if(A.numCols!=B.numRows)
	    throw new IllegalArgumentException("matrices have wrong sizes");
	double [][] tmpA = A != this ? A.re :                      MatrixOperations.copy(re)  ;
	double [][] tmpB = B != this ? B.re : ( A == this ? tmpA : MatrixOperations.copy(re) );
	newSize( A.numRows, B.numCols );
	MatrixOperations.times(tmpA,tmpB, re);
    }

     /** multiply <i>A</i> and <i>B</i> and store the result in <i>this</i> matrix */
    public void assignTimes(final IntegerMatrix A, final IntegerMatrix B){
	if(A.numCols!=B.numRows)
	    throw new IllegalArgumentException("matrices have wrong sizes");
	newSize( A.numRows, B.numCols );
	MatrixOperations.times(A.re,B.re, re);
    }

   /** multiply <i>A</i> and <i>B</i> and store the result in <i>this</i> matrix */
    public void assignTimes(final IntegerMatrix A, final RealMatrix B){
	if(A.numCols!=B.numRows)
	    throw new IllegalArgumentException("matrices have wrong sizes");
	double [][] tmpB = B != this ? B.re : MatrixOperations.copy(re);
	newSize( A.numRows, B.numCols );
	MatrixOperations.times( A.re, tmpB, re );
    }

    /** multiply <i>A</i> and <i>B</i> and store the result in <i>this</i> matrix */
    public void assignTimes(final RealMatrix A, final IntegerMatrix B){
	if(A.numCols!=B.numRows)
	    throw new IllegalArgumentException("matrices have wrong sizes");
	double [][] tmpA = A != this ? A.re : MatrixOperations.copy(re);
	newSize( A.numRows, B.numCols );
	MatrixOperations.times( tmpA, B.re, re );
    }

    /* ***** end of matrix-matrix multiplication ***** */

    /** multiply <i>this</i> and <i>V</i> and create and return a vector containing the result */
    public RealVector times(final RealVector V){
	RealVector res=new RealVector ( numRows );
	res.assignTimes( this, V );
	return res;
    }

    /** multiply <i>this</i> and <i>V</i> and create and return a vector containing the result */
    public RealVector times(final IntegerVector V){
	RealVector res=new RealVector ( numRows );
	res.assignTimes( this, V );
	return res;
    }

    /** multiply <i>this</i> and <i>V</i> and create and return a vector containing the result */
    public ComplexVector times(final ComplexVector V){
	ComplexVector res=new ComplexVector ( numRows );
	res.assignTimes( this, V );
	return res;
    }

   /** divide all entries by <i>b</i> and create and return a matrix containing the result */
    public RealMatrix divide( final double b ){
	RealMatrix result = new RealMatrix (this);
	result.assignDivide(b);
	return result;
    }

    /** divide all entries by <i>b</i> and create and return a matrix containing the result */
    public ComplexMatrix divide( final Field.Complex b ){
	ComplexMatrix result = new ComplexMatrix (this);
	result.assignDivide(b);
	return result;
    }

   /** divide all entries by <i>b</i>  */
    public void assignDivide(final double b){
	MatrixOperations.divide(re,b,re);
    }


    /** divide all entries of <i>a</i> by <i>b</i> and store the result in <i>this</i> matrix */
    public void assignDivide(final RealMatrix a, final double b){	
	newSize( a );
	MatrixOperations.divide( a.re, b, re );
    }

     /** divide all entries of <i>a</i> by <i>b</i> and store the result in <i>this</i> matrix */
    public void assignDivide( final double b,final RealMatrix a){	
	newSize( a );
	MatrixOperations.divide(  b, a.re, re );
    }

   /** divide all entries of <i>a</i> by <i>b</i> and store the result in <i>this</i> matrix */
    public void assignDivide(final IntegerMatrix a, final double b){	
	newSize( a );
	MatrixOperations.divide( a.re, b, re );
    }

    /** divide all entries of <i>a</i> by <i>b</i> and store the result in <i>this</i> matrix */
    public void assignDivide( final double b, final IntegerMatrix a ){	
	newSize( a );
	MatrixOperations.divide( b, a.re, re );
    }

  /** return <tt>true</tt> if <i>this</i> is symmetric */
    public boolean isSymmetric() {
	if( !isSquared() )
	    return false;	
	for( int i=0; i<numRows; i++ ) {
	    double[] Rre = re[i];
	    for( int j=0; j<i; j++ ) {
	        double dummy = Rre[j]-re[j][i];
		if( dummy*dummy > EPSILON )
		    return false;
	    }
	}
	return true;
    }

    /** return the sum of the squares of all entries */
    public double normSqr(){
	return MatrixOperations.normSqr(re);
    }

    /** return <tt>normSqr</tt> */
    public static double normSqr( final RealMatrix M ){
	return MatrixOperations.normSqr(M.re);
    }

     /** return the sum of the squares of all entries */
    public double norm(){
	return Math.sqrt(MatrixOperations.normSqr(re));
    }

    /** return <tt>norm</tt> */
    public static double norm( final RealMatrix M ){
	return Math.sqrt(MatrixOperations.normSqr(M.re));
    }

   /** return <tt>true</tt> if all entries are zero */
    public static boolean isZero( final RealMatrix M ) {
      return M.isZero();
    }

    /** return <tt>true</tt> if <i>this</i> equals <i>m</i> */
    public boolean equals( final RealMatrix m ){
	if(m.numRows!=numRows||m.numCols!=numCols) return false;
	for(int i=0; i<numRows; i++)
	    for(int j=0; j<numCols; j++)
		if ((m.re[i][j]-re[i][j])*(m.re[i][j]-re[i][j]) > EPSILON) return false;
	return true;
    }    

    /** return <tt>true</tt> if <i>this</i> equals <i>m</i> */
    public boolean equals( final RealMatrix m , final double diff){
	if(m.numRows!=numRows||m.numCols!=numCols) return false;
	for(int i=0; i<numRows; i++)
	    for(int j=0; j<numCols; j++)
		if ((m.re[i][j]-re[i][j])*(m.re[i][j]-re[i][j]) > diff) return false;
	return true;
    } 

    /** return <tt>true</tt> if <i>this</i> equals <i>o</i> */
    public boolean equals( final Object o ){
	try{
	    return equals( (RealMatrix)o ); 
	}
	catch(ClassCastException ex){
	    return false;
	}
    }

    /** creates a <b>String</b> containing this' matrix' entries */
    public String toString(){
	StringBuffer sb=new StringBuffer(300);
	sb.append("(");
	for(int i=0; i<numRows; i++){
	    sb.append('(');
	    if(numCols>0){
		sb.append(re[i][0]);
		for(int j=1; j<numCols; j++){
		    sb.append(", ");
		    sb.append(re[i][j]);
		}
	    }
	    sb.append(')');
	    if(i<numRows-1) sb.append(",\n");
	}
	sb.append(')');
	return sb.toString();
    }
    
    /** transpose <i>M</i> and create and return a matrix containing the result */
    public static RealMatrix transpose( final RealMatrix M ) {
        RealMatrix res = new RealMatrix ( M.numCols, M.numRows );	
	MatrixOperations.transpose( M.re, res.re );
	return res;
    }

    /** transpose <i>this</i> matrix create and return a matrix containing the result */
    public RealMatrix transpose() {
	RealMatrix res = new RealMatrix ( numCols, numRows );
	MatrixOperations.transpose( re, res.re );
	return res;
    }

    /** transpose <i>M</i> and store the result in <i>this</i> matrix */
    public void assignTranspose( final RealMatrix M ){
	RealMatrix tmpM = (  M.isSquared() || M != this ) ? M : new RealMatrix ( this );
	newSize( tmpM.numCols, tmpM.numRows );

	MatrixOperations.transpose( tmpM.re, re );
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
     * @see de.jtem.numericalMethods.algebra.linear.Determinant#compute(double[][]) Determinant
     */
    public double determinant(){
        if (numRows != numCols)
            throw new IllegalArgumentException("matrix must be square");
        return Determinant.compute(re);
    }

    /**
     * Computes and returns the inverse of this matrix.
     * <p />
     * This method temporary uses two additional matrices and one vector.
     * <p />
     * @return	The inverse.
     * @see de.jtem.numericalMethods.algebra.linear.Inversion#compute(double[][],double[][]) Inversion
     */
    public RealMatrix invert() {
        RealMatrix inv = new RealMatrix();
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
     * @see de.jtem.numericalMethods.algebra.linear.Inversion#compute(double[][],double[][]) Inversion
     */
    public static RealMatrix invert( final RealMatrix M )  {
        RealMatrix inv = new RealMatrix();
        inv.assignInvert(M);
        return inv;
    }

    /**
     * Computes the inverse of a matrix and stores the result in this matrix.
     * <p />
     * This method temporary uses one additional matrix and one vector.
     * <p />
     * @param	M   The matrix to invert.
     * @see de.jtem.numericalMethods.algebra.linear.Inversion#compute(double[][],double[][]) Inversion
     */
    public void assignInvert(final RealMatrix M){
        if (M.numRows != M.numCols)
            throw new IllegalArgumentException("matrix must be square");
        resize(M.numRows);
        if (!Inversion.compute(M.re, re))
            throw new IllegalArgumentException("matrix was singular");
    }

    /**
     * Inverts this matrix.
     * <p />
     * This method temporary uses one additional matrix and one vector.
     * <p />
     * @see de.jtem.numericalMethods.algebra.linear.Inversion#compute(double[][],double[][]) Inversion
     */
    public void assignInvert() {
        RealMatrix inv = new RealMatrix();
        inv.assignInvert(this);
        assign(inv);
    }

    /** assign all entries to zero */
    public void assignZero() {
	MatrixOperations.assign( re, 0 );	
    }
    
    /** assign all entries randomly */
    public void assignRandom(){
	MatrixOperations.random( re );
    }

    /** returns round of @param A */
    public static RealMatrix round( final RealMatrix A ) {
	RealMatrix result = new RealMatrix ( A.numRows, A.numCols );
	MatrixOperations.round( A.re, result.re );
	return result;
    }

   /** returns round of<em>this</em> */
    public RealMatrix round() {
	RealMatrix result = new RealMatrix ( numRows, numCols );
	MatrixOperations.round( re, result.re );
	return result;
    }

    /** assign <em>this</em> to round of @param A */
    public void assignRound( final RealMatrix A ) {
	newSize( A.numRows, A.numCols );
	MatrixOperations.round( A.re, re );
    }

     /** assign <em>this</em> to round of @this */
    public void assignRound() {
	MatrixOperations.round( re, re );
    }

    /** returns floor of @param A */
    public static RealMatrix floor( final RealMatrix A ) {
	RealMatrix result = new RealMatrix ( A.numRows, A.numCols );
	MatrixOperations.floor( A.re, result.re );
	return result;
    }

   /** returns floor of<em>this</em> */
    public RealMatrix floor() {
	RealMatrix result = new RealMatrix ( numRows, numCols );
	MatrixOperations.floor( re, result.re );
	return result;
    }

    /** floor <em>this</em> to round of @param A */
    public void assignFloor( final RealMatrix A ) {
	newSize( A.numRows, A.numCols );
	MatrixOperations.floor( A.re, re );
    }

    /** floor <em>this</em> to round of this */
    public void assignFloor() {
	MatrixOperations.floor( re, re );
    }

   /** creates real matrix  and assigns it with neg of <code>v</code> */
    public final static RealMatrix neg( final RealMatrix v ) {
	RealMatrix result = new RealMatrix ( v.numRows, v.numCols );
	MatrixOperations.neg( v.re, result.re );
	return result;
    }

    /** assign neg of <code>v</code> */
    public final void assignNeg( final RealMatrix v ) {
	newSize( v.numRows, v.numCols ,false );
	MatrixOperations.neg( v.re, re );
    }

   /** assign neg of <code>v</code> */
    public final void assignNeg( final IntegerMatrix v ) {
	newSize( v.numRows, v.numCols ,false );
	MatrixOperations.neg( v.re, re );
    }

    /** assign this to neg of <code>this</code> */
    public final void assignNeg() {
	MatrixOperations.neg( re, re );
    }
	
}
