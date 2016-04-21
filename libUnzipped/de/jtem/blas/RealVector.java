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
 * This class represents general double vectors.
 * 
 * Following the number package's philosophy, operations on double vectors 
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
public class RealVector extends AbstractVector
{
    private static final long serialVersionUID = 1L;
    
    public double[] re;

    /** create a zero vector */
    public RealVector() {
	re=new double[0];
	length=0;
    }

    /** create a vector with all entries equal to zero */
    public RealVector(int size){
	re=new double[size];
	length=size;
    }

    /** create a vector and initialize all entries with <i>initValue</i> */
    public RealVector(int size,double initValue){
	this( size );
	for(int ix=0; ix<size; ix++){
		re[ix]=initValue;
	    }
    }

    /** create copy of vector V */
    public RealVector( final RealVector V){
	re=(double[])V.re.clone();
	length=V.length;
    }  

    /** creates real version of vector V */
    public RealVector( final IntegerVector V){
	this( V.length );
	assign( V );
    }

    /** create a vector with initialized with <i>r</i> */
    public RealVector( final int[] r){
	this( r.length );
	VectorOperations.assign( r, re );
    }

    /** create a vector with initialized with <i>r</i> */
    public RealVector( final double[] r){
	re=(double[])r.clone();
	length = r.length;
    }

    void newSize( int newSize, boolean copyValues ) {

	if( newSize != length ) {
	    double [] newRe = new double[newSize];
       
	    int minSize = newSize < length ? newSize : length;

	    if( copyValues ) {
		System.arraycopy( re, 0, newRe, 0, minSize );
	    }

	    length = newSize;

	    re = newRe;
	}
    }

    /** return entry <i>index</i> */
    public double get( int index ) {
	return re[index];
    }
    
     /** set entry <i>index</i> to <i>r</i> */
    public void set(int index, double r) {
	re[index]=r;
    }
 
    /** assign <i>this<i/> vector to <i>V</i> */
    public void assign( final RealVector V ) {
    newSize( V.length ); 
	System.arraycopy( V.re, 0, re, 0, length );
    }
  
    /** assign <i>this<i/> vector to <i>V</i> */
    public void assign( final IntegerVector V ) {
	newSize( V.length );               
	VectorOperations.assign(V.re,re);
    }

    /** assign <i>this<i/> vector to <i>V</i> */
    public void assign( final double[] V ) {
	newSize(V.length);
	System.arraycopy( V, 0, re, 0, length );
    }
   
    /** assign <i>this<i/> vector to <i>V</i> */
    public void assign( final int[] V ) {
	newSize(V.length);
	VectorOperations.assign(V,re);
    }

   /** assign all entries to <i>r</i> */
    public void assign( final double r ) {
	VectorOperations.assign(re,r);
    }
  
    /** assign all entries to <i>r</i> */
    public void assign( final int r ) {
	VectorOperations.assign(re,r);
    }

    /** stores <i>this</i> in <i>diag</i>' diagonal while other entries remain unchanged */
    public void assignDiagonal( final RealMatrix diag ){
	diag.setDiagonal( this );
    }

    /** stores <i>this</i> in <i>diag</i>' diagonal while other entries remain unchanged */
    public void assignDiagonal( double[][] diag ){
        int numDiag = Math.min(diag.length, diag[0].length);
        if (size() != numDiag) throw new IllegalArgumentException("vector has wrong length");
	MatrixOperations.assignDiagonal( diag, re );
    }

    /** stores <i>this</i> in i. row of <i>row</i> while other entries remain unchanged.
     * i goes from 0 to length-1  
     */
    public void assignRow( RealMatrix row, int i ){
       	row.setRow( i, this);
    }
   
    /** stores <i>this</i> in i. row of <i>row</i>  while other entries remain unchanged 
     * i goes from 0 to length-1  
     */
    public void assignRow( double[][] row, int i ){
	if( row.length != length ) throw new IllegalArgumentException("sizes does not match");
       	if( row[0].length <= i ) throw new IllegalArgumentException("wrong second parameter");
	MatrixOperations.assignRow( row, re, i );
    }

    /** stores <i>this</i> in i. column of <i>col</i>  while other entries remain unchanged.
     * i goes from 0 to length-1  
     */
    public void assignCol( RealMatrix col, int i ){
       	col.setCol( i, this);
    }
   
    /** stores <i>this</i> in i. column of <i>col</i> while other entries remain unchanged
     * i goes from 0 to length-1  
     */
    public void assignCol( double[][] col, int i ){
       	if( col.length <= i ) throw new IllegalArgumentException("wrong second parameter");
       	if( col[0].length != length ) throw new IllegalArgumentException("sizes does not match");
	MatrixOperations.assignCol( col, re, i );
    }

    /** The following two methods return/set REFERENCES only ! */

    /** return a REFERENCE to <i>this</i> array */
    public double[] re() {
	return re;
    }

    /** set <i>this'</i> real part by REFERENCING <i>V</i>  */
    public void re( double[] V ) {
	checkLength( V.length );
	re = V;
    }

    /** create and return a block from <i>anIndex</i> to the bottom */
    public RealVector getBlock( int anIndex ) {
	return getBlock( anIndex, length-anIndex );
    }

    /** store a block from an <i>anIndex</i> in <i>V</i> */
    public void getBlock( int anIndex, RealVector V ) {
	getBlock( anIndex, length-anIndex, V);
    }

    /** create and return a block of size <i>aSize</i> beginning at <i>anIndex</i> */
    public RealVector getBlock( int anIndex, int aSize ) {
	RealVector V = new RealVector (aSize);
	getBlock(anIndex,aSize,V);
	return V;
    }

    /** store a block of size <i>aSize</i> beginning at <i>anIndex</i> in <i>V</i> */
    public RealVector getBlock( int anIndex, int aSize, RealVector V ) {

	checkIndex( anIndex );
	checkIndex( anIndex + aSize - 1 );

	V.newSize( aSize );

	System.arraycopy( re, anIndex, V.re, 0, aSize );
	
	return V;
    }
 
    /** store <i>V</i> as a block in <i>this</i> vector, beginning at <i>anIndex</i> */
    public void setBlock( int anIndex, RealVector V ) {
	setBlock( anIndex, V, V.length );
    }

    /** store a block of size <i>aSize</i> beginning at <i>V</i>'s top
        in <i>this</i> vector, start at <i>anIndex</i> */
    public void setBlock( int anIndex, RealVector V, int aSize ) {
	checkIndex( anIndex );
	checkIndex( anIndex + aSize - 1 );
	V.checkIndex( aSize -1 );
	System.arraycopy( V.re, 0, re, anIndex, aSize );
    }
		       			
       			
    /** create and return an array <b>double[]</b> containing <i>this</i>' entries */
    public double[] toArray() {
	return (double[])re.clone();
    }

    /** store <i>this</i>' entries in <i>res</i>*/
    public void toArray(double[] res){
      checkLength(res.length);
      System.arraycopy(re, 0, res, 0, length);
    }

    /** create and return a new vector containing <i>this</i>' entries */
    public RealVector copy(){
	return new RealVector (this);
    }

    /** return a clone */
    public Object clone() { return copy(); }
   
    /** add <i>V</i> and <i>this</i> and create and return a vector containing the result */
    public RealVector plus( final RealVector V ){
	checkShape( V ); 
	RealVector result = new RealVector (V.length);
	VectorOperations.plus( re, V.re, result.re );
	return result;
    }

    /** add <i>V</i> and <i>this</i> and create and return a vector containing the result */
    public RealVector plus( final IntegerVector V ){
	checkShape( V ); 
	RealVector result = new RealVector (V.length);
	VectorOperations.plus( re, V.re, result.re );
	return result;
    }

    /** subtract <i>V</i> from <i>this</i> and create and return a vector containing the result */
    public RealVector minus( final RealVector V ){
	checkShape( V ); 
	RealVector result = new RealVector (V.length);
	VectorOperations.minus( re, V.re, result.re );
	return result;
    }

      /** subtract <i>V</i> from <i>this</i> and create and return a vector containing the result */
    public RealVector minus( final IntegerVector V ){
	checkShape( V ); 
	RealVector result = new RealVector (V.length);
	VectorOperations.minus( re, V.re, result.re );
	return result;
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
    public void assignPlus( final RealVector V){
	checkShape(V);
	VectorOperations.plus(re,V.re,re);
    }

  /** add <i>V</i> and <i>this</i> and store result in <i>this</i> vector */
    public void assignPlus( final IntegerVector V){
	checkShape(V);
	VectorOperations.plus(re,V.re,re);
    }

    /** subtract <i>V</i> from <i>this</i> and store result in <i>this</i> vector */
    public void assignMinus( final RealVector V){
	checkShape(V);
	VectorOperations.minus(re,V.re,re);
    }
	
     /** subtract <i>V</i> from <i>this</i> and store result in <i>this</i> vector */
    public void assignMinus( final IntegerVector V){
	checkShape(V);
	VectorOperations.minus(re,V.re,re);
    }

   /**  this = V + W. */
    public void assignPlus( final RealVector V, final RealVector W){
	V.checkShape( W );
	newSize( V.length );
	VectorOperations.plus( V.re, W.re, re );
    }

    /**  this = V + W. */
    public void assignPlus( final IntegerVector V, final IntegerVector W){
	V.checkShape( W );
	newSize( V.length );
	int[] res = new int[V.length];
	VectorOperations.plus( V.re, W.re, res );
	VectorOperations.assign( res, re );
    }

     /**  this = V + W. */
    public void assignPlus( final RealVector V, final IntegerVector W){
	V.checkShape( W );
	newSize( V.length );
	VectorOperations.plus( V.re, W.re, re );
    }

     /**  this = V + W. */
    public void assignPlus( final IntegerVector V, final RealVector W){
	V.checkShape( W );
	newSize( V.length );
	VectorOperations.plus( V.re, W.re, re );
    }

    /** subtract <i>W</i> from <i>V</i> and store result in <i>this</i> vector */
    public void assignMinus( final RealVector V, final RealVector W){	
	V.checkShape( W );
	newSize( V.length );
	VectorOperations.minus( V.re, W.re, re );
    }

     /** subtract <i>W</i> from <i>V</i> and store result in <i>this</i> vector */
    public void assignMinus( final IntegerVector V, final RealVector W){	
	V.checkShape( W );
	newSize( V.length );
	VectorOperations.minus( V.re, W.re, re );
    }

    /** subtract <i>W</i> from <i>V</i> and store result in <i>this</i> vector */
    public void assignMinus( final RealVector V, final IntegerVector W){	
	V.checkShape( W );
	newSize( V.length );
	VectorOperations.minus( V.re, W.re, re );
    }

   /** subtract <i>W</i> from <i>V</i> and store result in <i>this</i> vector */
    public void assignMinus( final IntegerVector V, final IntegerVector W){	
	V.checkShape( W );
	newSize( V.length );
	int[] res = new int[V.length];	
	VectorOperations.minus( V.re, W.re, res );
	VectorOperations.assign( res, re );

    }

    /** return the dot product of <i>this</i> and <i>V</i> */
    public double dot( final IntegerVector V ){
	checkShape(V);
	return VectorOperations.times( re, V.re );
    }

    /** return the dot product of <i>this</i> and <i>V</i> */
    public double dot( final RealVector V){
	checkShape(V);
	return VectorOperations.times( re, V.re );
    }

    /** return the dot product of <i>this</i> and <i>V</i> */
    public static final double dot( final RealVector V, final IntegerVector W ){
	if( V.length != W.length )
	    throw new IllegalArgumentException("different vec sizes");
	return VectorOperations.times( V.re, W.re );
    }

    /** return the dot product of <i>this</i> and <i>V</i> */
    public static final double dot(  final IntegerVector W, final RealVector V ){
	if( V.length != W.length )
	    throw new IllegalArgumentException("different vec sizes");
	return VectorOperations.times( V.re, W.re );
    }

    /** multiply all entries with r and create and return a vector containing the result */
    public RealVector times( final double r)    {
	RealVector result = new RealVector ( length );
	VectorOperations.times( re, r, result.re );
	return result;
    }

     /** multiply all entries with r and create and return a vector containing the result */
    public RealVector times( final int r)    {
	RealVector result = new RealVector ( length );
	VectorOperations.times( re, r, result.re );
	return result;
    }

   /** multiply all entries with r */
    public void assignTimes( final double r){
	VectorOperations.times(re,r,re);
    }

   /** multiply all entries with r */
    public void assignTimes( final int r){
	VectorOperations.times(re,r,re);
    }

    /** multiply all entries of <i>V</i> with <i>r</i> and store the result in <i>this</i> vector */
    public void assignTimes(final RealVector V, final double r) {
	newSize( V.length );
	VectorOperations.times( V.re, r, re );
    }

    /** multiply all entries of <i>V</i> with <i>r</i> and store the result in <i>this</i> vector */
    public void assignTimes(final double r, final RealVector V) {
	assignTimes( V, r); 
    }  

     /** multiply all entries of <i>V</i> with <i>r</i> and store the result in <i>this</i> vector */
    public void assignTimes(final RealVector V, final int r) {
	newSize( V.length );
	VectorOperations.times( V.re, r, re );
    }

    /** multiply all entries of <i>V</i> with <i>r</i> and store the result in <i>this</i> vector */
    public void assignTimes( final int r, final RealVector V) {
	assignTimes( V, r); 
    }

    /** multiply all entries of <i>V</i> with <i>r</i> and store the result in <i>this</i> vector */
    public void assignTimes( final IntegerVector V, final double r ) {
	newSize( V.length );
	VectorOperations.times( V.re, r, re );
    }

    /** multiply all entries of <i>V</i> with <i>r</i> and store the result in <i>this</i> vector */
    public void assignTimes( final double r, final IntegerVector V ) {
	assignTimes( V, r); 
    } 


   /** creates a new RealVector W = this / r, (divide all entries of this by r). */
    public RealVector divide( final double r)    {
	RealVector result = new RealVector ( length );
	VectorOperations.divide( re, r, result.re);
	return result;
    }

    /**  this /= r, (divide all entries of this by r). */
    public void assignDivide( final double r) {
	VectorOperations.divide(re,r,re);
    }

    /**  this = V / r, (divide all entries of V by r). */
    public void assignDivide( final RealVector V, final double r) {
	newSize( V.length );
	VectorOperations.divide( V.re, r, re );
    }
    /**  this =  r / V */
    public void assignDivide( final double r, final RealVector V ){
	newSize( V.length );
	VectorOperations.divide( r, V.re, re ) ;
    }
    
    /**  this = V / r, (divide all entries of V by r). */
    public void assignDivide( final IntegerVector V, final double r) {
	newSize( V.length );
	VectorOperations.divide( V.re, r, re );
    }
    /**  this =  r / V */
    public void assignDivide( final double r, final IntegerVector V ){
	newSize( V.length );
	VectorOperations.divide( r, V.re, re ) ;
    }

    /** Creates and returns new RealVector that equals <code>this * M</code>. */  
    public RealVector times( final RealMatrix M ) {
	RealVector result = new RealVector ( M.getNumCols() );
	result.assignTimes( this, M );
	return result;
    }    

     /** Creates and returns new RealVector that equals <code>this * M</code>. */  
    public RealVector times( final IntegerMatrix M ) {
	RealVector result = new RealVector ( M.getNumCols() );
	result.assignTimes( this, M );
	return result;
    }

   /** this = M * V. */
    public void assignTimes( final RealMatrix M, final RealVector V ) {
	V.checkLength(M.numCols);
	double[] vRe = (this != V ) ? V.re : (double[])V.re.clone() ;
	newSize( M.numRows );
	MatrixOperations.times(M.re,vRe,re);
    }

    /** this = M * V. */
    public void assignTimes( final RealMatrix M, final IntegerVector V ) {
	V.checkLength(M.numCols);
	newSize( M.numRows );
	MatrixOperations.times( M.re, V.re, re);
    }

     /** this = M * V. */
    public void assignTimes( final IntegerMatrix M, final IntegerVector V ) {
	V.checkLength(M.numCols);
	newSize( M.numRows );
	MatrixOperations.times( M.re, V.re, re);
    }

    /** this = M * V. */
    public void assignTimes( final IntegerMatrix M, final RealVector V ) {
	V.checkLength(M.numCols);
	newSize( M.numRows );
	MatrixOperations.times( M.re, V.re, re);
    }

    /** this = V * M. */
    public void assignTimes( final RealVector V, final RealMatrix M ) {
	V.checkLength(M.numRows);
	double[] vRe = (this != V ) ? V.re : (double[])V.re.clone();

	newSize( M.numCols );
	MatrixOperations.times(vRe,M.re,re);
    }
 
    /** this = V * M. */
    public void assignTimes( final RealVector V, final IntegerMatrix M ) {
	V.checkLength(M.numRows);
	double[] vRe = (this != V ) ? V.re : (double[])V.re.clone();

	newSize( M.numCols );
	MatrixOperations.times(vRe,M.re,re);
    }

    /** this = V * M. */
    public void assignTimes( final IntegerVector V, final RealMatrix M ) {
	V.checkLength(M.numRows);
	newSize( M.numCols );
	MatrixOperations.times(V.re,M.re,re);
    }

    /** this = V * M. */
    public void assignTimes( final IntegerVector V, final IntegerMatrix M ) {
	V.checkLength(M.numRows);
	newSize( M.numCols );
	MatrixOperations.times(V.re,M.re,re);
    }

    /* return norm of <i>V</i> */
    public static double norm( RealVector V ) {
	return Math.sqrt(VectorOperations.normSqr(V.re));
    }

    /* return norm of <i>this</i> vector */
    public double norm() {
	return Math.sqrt(VectorOperations.normSqr(re));
    }

    /* return square of norm of <i>V</i> */
    public static double normSqr( final RealVector V ){
	return VectorOperations.normSqr(V.re);
    }

    /* return square of norm of <i>this</i> vector */
    public double normSqr() {	
	return VectorOperations.normSqr(re);
    }

    /** return <tt>true</tt> if all entries are zero */
    public static boolean isZero(RealVector V) {
      return V.isZero();
    }

    /** return <tt>true</tt> if <i>this</i> vector is equal to <i>v</i> */
    public boolean equals( final RealVector v ){
	final double [] vRe = v.re;
	
	if(v.length!=length) return false;
	for(int i=0; i<length; i++)
	    if( (vRe[i]-re[i])*(vRe[i]-re[i]) > EPSILON )
		return false;
	return true;
    }
    
    /** return <tt>true</tt> if <i>this</i> vector is equal to <i>v</i>.
     * <i>diff</i> is maximum diffrence allowed between them 
     */
    public boolean equals( final RealVector v, final double diff ){
	final double [] vRe = v.re;
	
	if(v.length!=length) return false;
	for(int i=0; i<length; i++)
	    if( (vRe[i]-re[i])*(vRe[i]-re[i]) > diff )
		return false;
	return true;
    }

    /** return <tt>true</tt> if <i>this</i> vector is equal to <i>o</i> */
    public boolean equals(Object o){
	try
	    {
		return equals( (RealVector)o );
	    }
	catch(ClassCastException ex)
	    {
		return false;
	    }
    }
   
    /** creates a <b>String</b> containing this' matrix' entries */
    public String toString(){
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

    /** normalize <i>this</i> vector */
    public void assignNormalize(){
	double n = norm();
	if (n<Double.MIN_VALUE) throw new IllegalArgumentException("Vector has length zero");
	assignDivide(n);
    }

    /** normalize <i>V</i> and store the result in <i>this</i> vector */
    public void assignNormalize( RealVector V ){
        assign(V);
	assignNormalize();
    }

    /** normalize <i>this</i> vector and create and return a vector containing the result */
    public RealVector normalize(){
	return normalize( this );
    }
	
    /** normalize <i>V</i> and create and return a vector containing the result */
    public static RealVector normalize( RealVector V ){
	RealVector result = new RealVector ( V );
	result.assignNormalize();
	return result;
    }
	
	
    /** assign all entries to zero */
    public void assignZero() {	
	VectorOperations.assign( re, 0 );
    }

    /** creates int vector and assigns it with round of <code>v</code> */
    public final static RealVector round( final RealVector v ) {
	RealVector result = new RealVector ( v.length );
	VectorOperations.round( v.re, result.re );
	return result;
    }

    /** assign round of <code>v</code> */
    public final void assignRound( final RealVector v ) {
	newSize( v.length );
	VectorOperations.round( v.re, re );
    }

    /** assign this to round of <code>this</code> */
    public final void assignRound() {
	VectorOperations.round( re, re );
    }

    /** creates int vector and assigns it with round of <code>v</code> */
    public final static RealVector floor( final RealVector v ) {
	RealVector result = new RealVector ( v.length );
	VectorOperations.floor( v.re, result.re );
	return result;
    }

    /** assign round of <code>v</code> */
    public final void assignFloor( final RealVector v ) {
	newSize( v.length );
	VectorOperations.floor( v.re, re );
    }

    /** assign this to round of <code>this</code> */
    public final void assignFloor() {
	VectorOperations.floor( re, re );
    }

    /** creates real vector and assigns it with neg of <code>v</code> */
    public final static RealVector neg( final RealVector v ) {
	RealVector result = new RealVector ( v.length );
	VectorOperations.neg( v.re, result.re );
	return result;
    }

    /** assign neg of <code>v</code> */
    public final void assignNeg( final RealVector v ) {
	newSize( v.length );
	VectorOperations.neg( v.re, re );
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
	
    /** assign all entries randomly */
    public void assignRandom() {
      for( int i = 0; i<length; i++ ) re[i] = Math.random(); 
      
    }

}

   









