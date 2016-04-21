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

import de.jtem.numericalMethods.algebra.linear.MatrixOperations;
import de.jtem.numericalMethods.algebra.linear.VectorOperations;

/**
 * This class represents general integerVectors.
 * 
 * Following the number package's philosophy, operations on int vectors 
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
 *
 * @author	marina,vitali
 */
public class IntegerVector extends AbstractVector
{
    private static final long serialVersionUID = 1L;
    
    public int[] re;

    /** create a zero vector */
    public IntegerVector(){
	re=new int[0];
	length=0;
    }

    /** create a vector with all entries equal to zero */
    public IntegerVector(int size){
	re=new int[size];
	length=size;
    }

    /** create a vector and initialize all entries with <i>initValue</i> */
    public IntegerVector(int size,int initValue){
	this( size );
	for(int ix=0; ix<size; ix++){
		re[ix]=initValue;
	    }
    }

    /** create copy of vector V */
    public IntegerVector( final IntegerVector V){
	re=(int[])V.re.clone();
	length=V.length;
    }

    /** create a vector with initialized with <i>r</i> */
    public IntegerVector( final int[] r){
	re=(int[])r.clone();
	length=r.length;
    }

    void newSize( int newSize, boolean copyValues ) {
	if( newSize != length ) {
	    int [] newRe = new int[newSize];
       
	    int minSize = newSize < length ? newSize : length;

	    if( copyValues ) {
		System.arraycopy( re, 0, newRe, 0, minSize );
	    }
	    length = newSize;
	    re = newRe;
	}
    }

    /** return entry <i>index</i> */
    public int get( int index ) {
	return re[index];
    }
    
     /** set entry <i>index</i> to <i>r</i> */
    public void set(int index, int r) {
	re[index]=r;
    }

    /** assign <i>this<i/> vector to <i>V</i> */
    public void assign( final IntegerVector V ) {
	if( this != V ) {                      
	    newSize( V.length );               
	    System.arraycopy( V.re, 0, re, 0, length );
	}
    }

    /** assign <i>this<i/> vector to <i>V</i> */
    public void assign( final int[] V ) {
	newSize(V.length);
	System.arraycopy( V, 0, re, 0, length );
    }
   
        
    /** assign all entries to <i>r</i> */
    public void assign( final int r ) {
	VectorOperations.assign( re, r );
    }
    
    /** stores <i>this</i> in <i>diag</i>' diagonal while other entries remain unchanged */
    public void assignDiagonal( IntegerMatrix diag ){
	diag.setDiagonal( this );
    }

    /** stores <i>this</i> in <i>diag</i>' diagonal while other entries remain unchanged */
    public void assignDiagonal(int[][] diag ){
        int numDiag = Math.min(diag.length, diag[0].length);
        if (size() != numDiag) throw new IllegalArgumentException("vector has wrong length");
	MatrixOperations.assignDiagonal( diag, re );
    }

    /** stores <i>this</i> in i. row of <i>row</i> while other entries remain unchanged.
     * i goes from 0 to length-1  
     */
    public void assignRow( IntegerMatrix row, int i ){
	row.setRow( i , this );
    }
   
    /** stores <i>this</i> in i. row of <i>row</i>  while other entries remain unchanged 
     * i goes from 0 to length-1  
     */
    public void assignRow( int[][] row, int i ){
	if( row.length != length ) throw new IllegalArgumentException("sizes does not match");
       	if( row[0].length <= i ) throw new IllegalArgumentException("wrong second parameter");
	MatrixOperations.assignRow( row, re, i );
    }

    /** stores <i>this</i> in i. column of <i>col</i>  while other entries remain unchanged.
     * i goes from 0 to length-1  
     */
    public void assignCol( IntegerMatrix col, int i ){
	col.setCol( i, this );
    }
   
    /** stores <i>this</i> in i. column of <i>col</i> while other entries remain unchanged
     * i goes from 0 to length-1  
     */
    public void assignCol( int[][] col, int i ){
       	if( col.length <= i ) throw new IllegalArgumentException("wrong second parameter");
       	if( col[0].length != length ) throw new IllegalArgumentException("sizes does not match");
	MatrixOperations.assignCol( col, re, i );
    }

    /** return a REFERENCE to <i>this</i> array */
    public int[] re() {
	return re;
    }

    /** set <i>this'</i> real part by REFERENCING <i>V</i>  */
    public void re(int[] V) {
	checkLength( V.length );
	re = V;
    }

    /** create and return a block from <i>anIndex</i> to the bottom */
    public IntegerVector getBlock( int anIndex ) {
	return getBlock( anIndex, length-anIndex );
    }

    /** store a block from an <i>anIndex</i> in <i>V</i> */
    public void getBlock( int anIndex, IntegerVector V ) {
	getBlock( anIndex, length-anIndex, V);
    }

    /** create and return a block of size <i>aSize</i> beginning at <i>anIndex</i> */
    public IntegerVector getBlock( int anIndex, int aSize ) {
	IntegerVector V = new IntegerVector (aSize);
	getBlock(anIndex,aSize,V);
	return V;
    }

    /** store a block of size <i>aSize</i> beginning at <i>anIndex</i> in <i>V</i> */
    public IntegerVector getBlock( int anIndex, int aSize, IntegerVector V ) {

	checkIndex( anIndex );
	checkIndex( anIndex + aSize - 1 );
	V.newSize( aSize );
	System.arraycopy( re, anIndex, V.re, 0, aSize );
	
	return V;
    }
 
    /** store <i>V</i> as a block in <i>this</i> vector, beginning at <i>anIndex</i> */
    public void setBlock( int anIndex, IntegerVector V ) {
	setBlock( anIndex, V, V.length );
    }

    /** store a block of size <i>aSize</i> beginning at <i>V</i>'s top
        in <i>this</i> vector, start at <i>anIndex</i> */
    public void setBlock( int anIndex, final IntegerVector V, int aSize ) {
	checkIndex( anIndex );
	checkIndex( anIndex + aSize - 1 );
	V.checkIndex( aSize -1 );
	System.arraycopy( V.re, 0, re, anIndex, aSize );
    }	
       			
    /** create and return an array <b>int[]</b> containing <i>this</i>' entries */
    public int[] toArray()    {
	return (int[])re.clone();
    }

    /** store <i>this</i>' entries in <i>res</i>*/
    public void toArray(final int[] res)   {
      checkLength(res.length);
      System.arraycopy(re, 0, res, 0, length);
    }

    /** create and return a new vector containing <i>this</i>' entries */
    public IntegerVector copy()    {
	return new IntegerVector (this);
    }

    /** return a clone */
    public Object clone() { return copy(); }
    

    /** add <i>V</i> and <i>this</i> and create and return a vector containing the result */
    public IntegerVector plus( final IntegerVector V )    {
	checkShape( V );
	IntegerVector res=new IntegerVector ( V.length );
	VectorOperations.plus( re, V.re, res.re);
	return res;
    }

    /** subtract <i>V</i> from <i>this</i> and create and return a vector containing the result */
    public IntegerVector minus( final IntegerVector V )  {
	checkShape( V );
	IntegerVector res=new IntegerVector ( V.length );
	VectorOperations.minus( re, V.re, res.re);
	return res;
    }

     /** add <i>V</i> and <i>this</i> and create and return a vector containing the result */
    public RealVector plus( final RealVector V )    {
	checkShape( V );
	RealVector res=new RealVector ( V.length );
	VectorOperations.plus( re, V.re, res.re);
	return res;
    }

    /** subtract <i>V</i> from <i>this</i> and create and return a vector containing the result */
    public RealVector minus( final RealVector V )  {
	checkShape( V );
	RealVector res=new RealVector ( V.length );
	VectorOperations.minus( re, V.re, res.re);
	return res;
    }

    /** add <i>V</i> and <i>this</i> and create and return a vector containing the result */
    public ComplexVector plus( final ComplexVector V )    {
	checkShape( V );
	ComplexVector res=new ComplexVector ( V.length );
	VectorOperations.plus( re, V.re, res.re);
	return res;
    }

    /** subtract <i>V</i> from <i>this</i> and create and return a vector containing the result */
    public ComplexVector minus( final ComplexVector V )  {
	checkShape( V );
	ComplexVector res=new ComplexVector ( V.length );
	VectorOperations.minus( re, V.re, res.re);
	return res;
    }


   /** add <i>V</i> and <i>this</i> and store result in <i>this</i> vector */
    public void assignPlus( final IntegerVector V)  {
	checkShape(V);
	VectorOperations.plus( re, V.re, re );
    }

    /** subtract <i>V</i> from <i>this</i> and store result in <i>this</i> vector */
    public void assignMinus( final IntegerVector V)  {
	checkShape(V);
	VectorOperations.minus( re, V.re, re );
    }
	
    /**  this = V + W. */
    public void assignPlus( final IntegerVector V, final IntegerVector W) {	
	V.checkShape(W);
	newSize( V.length );
	VectorOperations.plus( V.re, W.re, re );
    }

    /** subtract <i>W</i> from <i>V</i> and store result in <i>this</i> vector */
    public void assignMinus( final IntegerVector V, final IntegerVector W) {	
	V.checkShape(W);
	newSize( V.length );
	VectorOperations.minus( V.re, W.re, re );
    }

    /** return the dot product of <i>this</i> and <i>V</i> */
    public int dot( final IntegerVector V ) {
	checkShape(V);
	return VectorOperations.times(re,V.re);
    }

    /** return the dot product of <i>this</i> and <i>V</i> */
    public static final double dot(  final IntegerVector W, final IntegerVector V ){
	if( V.length != W.length )
	    throw new IllegalArgumentException("different vec sizes");
	return VectorOperations.times( V.re, W.re );
    }


    /** multiply all entries with r and create and return a vector containing the result */
    public IntegerVector times(final int r){
	IntegerVector res=new IntegerVector ( length );
	VectorOperations.times( re, r, res.re );
	return res;
    }

    /** multiply all entries with r */
    public void assignTimes(final int r ){
	VectorOperations.times(re,r,re);
    }

    /** multiply all entries of <i>V</i> with <i>r</i> and store the result in <i>this</i> vector */
    public void assignTimes(final IntegerVector V,final int r) {
	newSize( V.length );
	VectorOperations.times(V.re,r,re);
    }

    /** multiply all entries of <i>V</i> with <i>r</i> and store the result in <i>this</i> vector */
    public void assignTimes(final int r, final IntegerVector V) {
	assignTimes(V,r) ;
    }
    
    /** creates a new IntegerVector W = this / r, (divide all entries of this by r). */
    public IntegerVector divide(final int r )    {
	IntegerVector res=new IntegerVector ( length );
	VectorOperations.divide( re, r, res.re );
	return res;
    }

    /**  this /= r, (divide all entries of this by r). */
    public void assignDivide(final int r) {
	VectorOperations.divide(re,r,re);
    }

    /**  this = V / r, (divide all entries of V by r). */
    public void assignDivide(final IntegerVector V, final int r) {
	newSize( V.length );
	VectorOperations.divide( V.re, r, re ) ;
    }
    /**  this = r / V */
    public void assignDivide( final int r, final IntegerVector V) {
	newSize( V.length );
	VectorOperations.divide(r,V.re,re);
    }
  
    /** Creates and returns new IntegerVector that equals <code>this * M</code>. */  
    public IntegerVector times( final IntegerMatrix M ) {
	IntegerVector result = new IntegerVector ( M.getNumCols() );
	result.assignTimes( this, M );
	return result;
    }    

    /** this = M * V. */
    public void assignTimes( final IntegerMatrix M, final IntegerVector V ) {
	V.checkLength(M.numCols);
	int[] vRe = (this != V) ? V.re : (int[])V.re.clone() ;

	newSize(M.numRows);
	MatrixOperations.times(M.re,vRe,re);
    }

    /** this = V * M. */
    public void assignTimes( final IntegerVector V, final IntegerMatrix M ) {
	V.checkLength(M.numRows);
	int[] vRe = (this != V) ? V.re : (int[])V.re.clone() ;
	newSize(M.numCols);
	MatrixOperations.times(vRe,M.re,re);
    }

    /* return norm of <i>V</i> */
    public static double norm( final IntegerVector V ) {
	return Math.sqrt(VectorOperations.normSqr(V.re));
    }

    /* return norm of <i>this</i> vector */
    public double norm() {
	return Math.sqrt(VectorOperations.normSqr(re));
    }

    /* return square of norm of <i>V</i> */
    public static double normSqr( final IntegerVector V ) {
	return VectorOperations.normSqr(V.re);
    }

    /* return square of norm of <i>this</i> vector */
    public double normSqr() {
	return VectorOperations.normSqr(re);
    }

    /** return <tt>true</tt> if all entries are zero */
    public static boolean isZero( final IntegerVector V ) {
      return V.isZero();
    }
    

    /** return <tt>true</tt> if <i>this</i> vector is equal to <i>v</i> */
    public boolean equals( final IntegerVector v ) {
	if(v.length!=length) return false;
	for(int i=0; i<length; i++)
	    if( (v.re[i]!=re[i]) )
		return false;
	return true;	
    }


    /** return <tt>true</tt> if <i>this</i> vector is equal to <i>o</i> */
    public boolean equals( Object o ) {
	try{
	     return equals( (IntegerVector)o );
	    }
	catch(ClassCastException ex)
	    {
		return false;
	    }
    }

    /** creates a <b>String</b> containing this' matrix' entries */
    public String toString()  {
	StringBuffer sb=new StringBuffer(300);
	sb.append("(");
	if(length>0)
	    {
		sb.append(re[0]);
		for(int j=1; j<length; j++)
		    {
			sb.append(", ");
			sb.append(re[j]);
		    }
	    }
	sb.append(')');
	return sb.toString();
    }

    /** assign all entries to zero */
    public void assignZero() {
	VectorOperations.assign( re, 0 );
    }

    /** creates int vector and assigns it with round of <code>v</code> */
    public final static IntegerVector round( final RealVector v ) {
	IntegerVector result = new IntegerVector ( v.length );
	VectorOperations.round( v.re, result.re );
	return result;
    }

    /** assign round of <code>v</code> */
    public final void assignRound( final RealVector v ) {
	newSize( v.length );
	VectorOperations.round( v.re, re );
    }
  
    /** creates int vector and assigns it with round of <code>v</code> */
    public final static IntegerVector floor( final RealVector v ) {
	IntegerVector result = new IntegerVector ( v.length );
	VectorOperations.floor( v.re, result.re );
	return result;
    }

    /** assign round of <code>v</code> */
    public final void assignFloor( final RealVector v ) {
	newSize( v.length );
	VectorOperations.floor( v.re, re );
    }	
    public static java.util.Random iRandom = new java.util.Random();

    /** assign all entries randomly */
    public void assignRandom() {
      for( int i = 0; i<length; i++ ) re[i] = iRandom.nextInt(1024);
    }

    /** creates int vector and assigns it with neg of <code>v</code> */
    public final static IntegerVector neg( final IntegerVector v ) {
	IntegerVector result = new IntegerVector ( v.length );
	VectorOperations.neg( v.re, result.re );
	return result;
    }

    /** assign neg of <code>v</code> */
    public final void assignNeg( final IntegerVector v ) {
	newSize( v.length );
	VectorOperations.neg( v.re, re );
    }

    /** assign this to neg of <code>this</code> */
    public final void assignNeg() {
	VectorOperations.neg( re, re );
    }  
    

}

   









