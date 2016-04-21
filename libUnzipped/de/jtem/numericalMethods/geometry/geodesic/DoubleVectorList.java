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

package de.jtem.numericalMethods.geometry.geodesic;

/**
 * Simple list of double array which have all the same length.
 * @author schmies
 */
public class DoubleVectorList {


	double [][] data;
	int size;
	final int dim;
	
	public DoubleVectorList( int dim, int capacity ) {
		this.data =new double[capacity][dim];
		this.size=0;
		this.dim =dim;
	}

	public DoubleVectorList( int dim ) {
		this( dim, 100 );
	}
	
	public DoubleVectorList( double [][] data ) {
		size=data.length;
		dim=data[0].length;
		this.data = new double[ data.length + (data.length+10 )/2][dim];
		System.arraycopy( data, 0, this.data, 0, size );
	}
	
	public int size() {
		return size;
	}
	
	public int capacity() {
		return data.length;
	}
	
	public double[] get( int index ) {
		if( index<0 || index>=size) 
			throw new IllegalArgumentException("index out of range");
		return data[index];
	}
	
	/**
	 * 
	 * @return new last double-vector, which might not be all zeros 
	 */
	public double [] getNewLast() {
		if( size==data.length ) {
			setCapacity( data.length + (data.length+10 )/2);
		}
		return data[size++];
	}
	
	public void set( int index, double [] value ) {
		if( index<0 || index>size || value.length != dim ) 
			throw new IllegalArgumentException("index out of range");
		if( value.length != dim ) 
			throw new IllegalArgumentException("vector has wrong dimension");
		if( size==index) {
			if( size==data.length ) {
				setCapacity( data.length + (data.length+10 )/2);
			}
			data[size++] =value;
		}
		
		data[index] = value;
	}
	
	public void add( double [] value ) {
		if( size==data.length ) {
			setCapacity( data.length + (data.length+10 )/2);
		}
		data[size++] = value;
	}
	
	public void remove( int index ) {
		if( index<0 || index>=size) 
			throw new IllegalArgumentException("index out of range");
		System.arraycopy( data, index+1, data, index, size-index-1 );
		data[data.length-1] = new double[dim];
		size--;
	}
	
	public void insert( int index, double [] value ) {
		if( index<0 || index>size ) 
			throw new IllegalArgumentException("index out of range");
		if( value.length != dim ) 
			throw new IllegalArgumentException("vector has wrong dimension");
		if( size==data.length ) {
			setCapacity( data.length + (data.length+10 )/2);
		}
		System.arraycopy( data, index, data, index+1, size-index );
		data[index]=value;
		size++;
	}
	
	public void setCapacity( int capacity ) {
		if( capacity<size)
			throw new IllegalArgumentException( "capacity must exeed size of list");
		double [][] newData = new double[capacity][];
		System.arraycopy( data, 0, newData, 0, Math.min( newData.length, data.length ) );
		for( int i=data.length; i<newData.length; i++)
			newData[i] = new double[dim];
		data = newData;
		
	}
	
	public void removeAll() {
		size=0;
	}
	
	void appendVector( double [] vector, StringBuffer sb ) {
		sb.append("(");
		if(size>0)
		{
			sb.append(data[0]);
			for(int j=1; j<size; j++)
			{
				sb.append(", ");
				sb.append(data[j]);
			}
		}
		sb.append(')');
	}
	
	public String toString()  {
		StringBuffer sb=new StringBuffer(300);
		sb.append("(");
		if(size>0)
		{
			sb.append(data[0]);
			for(int j=1; j<size; j++)
			{
				sb.append(", ");
				appendVector(data[j], sb );
			}
		}
		sb.append(')');
		return sb.toString();
	}
}
