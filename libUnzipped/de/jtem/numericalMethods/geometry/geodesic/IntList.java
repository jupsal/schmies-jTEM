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
 * Simple list of integers.
 * @author schmies
 *
 */
public final class IntList {

	int [] data;
	int size;

	public IntList( int capacity ) {
		this.data =new int[capacity];
		this.size=0;
	}

	public IntList() {
		this( 100 );
	}
	
	public IntList( int [] data ) {
		this.data = new int[ data.length + (data.length+10 )/2];
		size=data.length;
		System.arraycopy( data, 0, this.data, 0, size );
	}
	
	public int size() {
		return size;
	}
	
	public int capacity() {
		return data.length;
	}
	
	public int get( int index ) {
		if( index<0 || index>=size) 
			throw new IllegalArgumentException("index out of range");
		return data[index];
	}
	
	public void set( int index, int value ) {
		if( index<0 || index>size) 
			throw new IllegalArgumentException("index out of range");
		if( size==index) {
			if( size==data.length ) {
				setCapacity( data.length + (data.length+10 )/2);
			}
			data[size++] =value;
		}
			
		data[index] = value;
	}
	
	public void add( int value ) {
		if( size==data.length ) {
			setCapacity( data.length + (data.length+10 )/2);
		}
		data[size++] = value;
	}
	
	public void remove( int index ) {
		if( index<0 || index>=size) 
			throw new IllegalArgumentException("index out of range");
		System.arraycopy( data, index+1, data, index, size-index-1 );
		size--;
	}
	
	public void insert( int index, int value ) {
		if( index<0 || index>size) 
			throw new IllegalArgumentException("index out of range");
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
		int [] newData = new int[capacity];
		System.arraycopy( data, 0, newData, 0, size );
		data = newData;
	}
	
	public void removeAll() {
		size=0;
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
					sb.append(data[j]);
				}
			}
			sb.append(')');
			return sb.toString();
		}
}
