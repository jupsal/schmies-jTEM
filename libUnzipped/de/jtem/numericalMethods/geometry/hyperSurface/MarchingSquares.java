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

package de.jtem.numericalMethods.geometry.hyperSurface;

/**
 * Computes isolines on a rectangular domain.
 * 
 * @author schmies
 */
public class MarchingSquares {

	public static class IndexedLineSet {
		double [] points;
		int [] index;
		int [] length;
	
		IndexedLineSet( double [] points, int[] index, int [] length ) {
			this.points = points;
			this.index = index;
			this.length =length;
		}
		
		public int[] getIndex() {
			return index;
		}

		public int[] getLength() {
			return length;
		}

		public double[] getPoints() {
			return points;
		}
		
		
	}
	
	static abstract class RealFunctionOnTupleIndex {
		int xSize, ySize;

		abstract double valueAt( int i, int j );
	}
	
	public static interface  RealFunctionOnReal2 {
		public double valueAt( double x, double y );
	}
	

	public static RealFunctionOnTupleIndex  scanner( final RealFunctionOnReal2 f,
			final double xMin, final double xMax,
			final int xDiscr, final double yMin, 
			final double yMax, final int yDiscr ) {

		return new RealFunctionOnTupleIndex() {
			{
				xSize = xDiscr;
				ySize = yDiscr;
			}

			final double xs = 1./(xSize-1);
			final double ys = 1./(ySize-1);
					
			double valueAt( int i, int j ) {
				final double x = xMin * ( 1 - i * xs ) + xMax * i * xs;
				final double y = yMin * ( 1 - j * ys ) + yMax * j * ys;

				
				return f.valueAt(x,y);
			}
		};
	}
	
	public static double [][] scan( final RealFunctionOnReal2 f,
			final double xMin, final double xMax,
			final int xDiscr, final double yMin, 
			final double yMax, final int yDiscr ) {
		
		double [][] scan = new double[xDiscr][yDiscr];
		
		scan( scan, f, xMin, xMax, yMin, yMax );
		
		return scan;
	}
	
	public static void scan(  double [][] scan, final RealFunctionOnReal2 f,
			final double xMin, final double xMax, 
			final double yMin, final double yMax ) {
		
		final int xDiscr =scan.length;
		final int yDiscr = scan[0].length;
		
		RealFunctionOnTupleIndex scanner = scanner( f, xMin, xMax, xDiscr, yMin, yMax, yDiscr );
		
		for( int i=0; i<xDiscr; i++ ) {
			for( int j=0; j<yDiscr; j++ ) {
				scan[i][j] = scanner.valueAt(i,j);
			}
		}
	}
	
	public static IndexedLineSet compute( final RealFunctionOnReal2 f, final double level, final double xMin, 
									final double xMax,
									final int xDiscr, final double yMin, 
									final double yMax, final int yDiscr ) {

		RealFunctionOnTupleIndex scanner = scanner( f, xMin, xMax, xDiscr, yMin, yMax, yDiscr );
		
		return compute( scanner, level, xMin, xMax, yMin, yMax );
	}
	
	
	public static IndexedLineSet compute( RealFunctionOnTupleIndex f, double level,
			double xMin, double xMax, double yMin, double yMax ) {

		final int xSize = f.xSize;
		final int ySize = f.ySize;

		final int[][] xIndex = new int[xSize][ySize - 1];
		final int[][] yIndex = new int[xSize - 1][ySize];

		int numOfPoints = 0;

		double[] point = new double[100];

		// computing x-lines-intersections
		for (int i = 0; i < xSize; i++) {
			double p = f.valueAt(i,0) - level;

			for (int j = 1; j < ySize; j++) {
				double q = f.valueAt(i,j)-level;

				if (p * q < 0) {

					final double t = Math.abs(p) / (Math.abs(p) + Math.abs(q));

					final double y = (1 - t) * (j - 1) + t * j;

					point = setPoint(point, numOfPoints, i, y);

					xIndex[i][j - 1] = numOfPoints++;
				}
				else {
					xIndex[i][j - 1] = -1;
				}

				p = q;
			}
		}

		// computing y-lines-intersections
		for (int j = 0; j < ySize; j++) {
			double p = f.valueAt(0,j)-level;

			for (int i = 1; i < xSize; i++) {
				double q = f.valueAt(i,j) -level;

				if (p * q < 0) {

					final double t = Math.abs(p) / (Math.abs(p) + Math.abs(q));

					final double x = (1 - t) * (i - 1) + t * i;

					point = setPoint(point, numOfPoints, x, j);

					yIndex[i - 1][j] = numOfPoints++;
				}
				else {
					yIndex[i - 1][j] = -1;
				}

				p = q;
			}
		}

		point = resize(point, 2 * numOfPoints);

		//if( xMin!=0&&xMax!=1&&yMin!=0&&yMax!=1) {
			for( int i=0, j=0; i<numOfPoints; i++) {
				point[j] = xMin*(1-point[j]/(xSize-1)) + xMax*(point[j]/(xSize-1)); j++;
				point[j] = yMin*(1-point[j]/(ySize-1)) + yMax*(point[j]/(ySize-1)); j++;	
			}
		//}
		// compute edges
		int[] index = new int[100];

		int[] intersection = new int[4];

		int numOfEdges = 0;

		// march all cells
		for (int i = 0; i < xSize - 1; i++) {
			for (int j = 0; j < ySize - 1; j++) {

				intersection[0] = xIndex[i][j];
				intersection[1] = xIndex[i + 1][j];
				intersection[2] = yIndex[i][j];
				intersection[3] = yIndex[i][j + 1];

				// connect all to oneanother
				for (int ii = 0; ii < 4; ii++) {
					if (intersection[ii] != -1) {
						for (int jj = ii + 1; jj < 4; jj++) {
							if (intersection[jj] != -1) {
								index = setEdge(index, numOfEdges++,
										intersection[ii], intersection[jj]);
							}
						}
					}
				}
			}
		}

		index = resize(index, 2 * numOfEdges);

		// init the IndexedLineSet
		int[] count = new int[numOfEdges];

		for (int i = 0; i < numOfEdges; i++) {
			count[i] = 2;
		}

		return new IndexedLineSet(point, index, count);
	}


	public static IndexedLineSet compute( final double[][] field, double level,
			double xMin, double xMax, double yMin, double yMax ) {
	
		RealFunctionOnTupleIndex f = new RealFunctionOnTupleIndex() {
			{
				xSize = field.length;
				ySize = field[0].length;
			}

			double valueAt( int i, int j ) {
				return field[i][j];
			}
		};
		
		return compute( f, level, xMin, xMax, yMin, yMax );
	}

	private static double[] setPoint(double[] point, int pos, double x, double y) {

		if (point.length <= 2 * pos) {

			double[] newPoint = new double[2 * point.length + 10];

			System.arraycopy(point, 0, newPoint, 0, point.length);

			point = newPoint;
		}

		point[2 * pos] = x;
		point[2 * pos + 1] = y;

		return point;
	}

	private static double[] resize(double[] point, int size) {

		if (point.length == size) {
			return point;
		}

		double[] newPoint = new double[size];

		System.arraycopy(point, 0, newPoint, 0, Math.min(size, point.length));

		return newPoint;
	}

	private static int[] setEdge(int[] index, int pos, int p, int q) {

		if (index.length <= 2 * pos + 2) {

			int[] newIndex = new int[2 * index.length + 10];

			System.arraycopy(index, 0, newIndex, 0, index.length);

			index = newIndex;
		}

		index[2 * pos] = p;
		index[2 * pos + 1] = q;

		return index;
	}

	private static int[] resize(int[] index, int size) {

		if (index.length == size) {
			return index;
		}

		int[] newIndex = new int[size];

		System.arraycopy(index, 0, newIndex, 0, Math.min(size, index.length));

		return newIndex;
	}

	public static double[][] removeDirt(double[][] field, int x, int y,
			double min, double max, double ratio,
			double value) {

		final int xSize = field.length;
		final int ySize = field[0].length;

		final double[][] out = (double[][]) field.clone();

		final int xSizeMinusX = xSize - x;
		final int ySizeMinusY = ySize - y;

		final int count = (int) (x * y * ratio);

		for (int i = 0; i < xSizeMinusX; i++) {
			for (int j = 0; j < ySizeMinusY; j++) {

				if (getNumOfDirtyValues(field, i, j, x, y, min, max) > count) {
					removeDirt(out, i, j, x, y, value);
				}
			}
		}

		return out;
	}

	public static int getNumOfDirtyValues(double[][] field,
			int i, int j, int x, int y,
			double min, double max) {

		int numOfDirtyValues = 0;

		for (int I = 0, ii = i; I < x; I++, ii++) {
			for (int J = 0, jj = j; J < y; J++, jj++) {

				double v = field[ii][jj];

				if (v < min || v > max) {
					numOfDirtyValues++;
				}
			}
		}

		return numOfDirtyValues;
	}

	public static void removeDirt(double[][] field,
			int i, int j, int x, int y, double value) {

		for (int I = 0, ii = i; I < x; I++, ii++) {
			for (int J = 0, jj = j; J < y; J++, jj++) {
				field[ii][jj] = value;
			}
		}
	}

	public static double[][] removeDirt(double[][] field, int l,
			int maxNumOfSignChanges, double value) {

		final double[][] out = (double[][]) field.clone();

		removeDirt(field, out, l, maxNumOfSignChanges, value);

		return out;
	}

	public static void removeDirt(double[][] field, double[][] out, int l,
			int maxNumOfSignChanges, double value) {

		final int xSize = field.length;
		final int ySize = field[0].length;

		for (int i = 0; i < xSize; i++) {
			System.arraycopy(field[i], 0, out[i], 0, ySize);

		}
		final int xSizeMinusL = xSize - l;
		final int ySizeMinusL = ySize - l;

		for (int i = 0; i < xSize; i++) {
			for (int j = 0; j < ySizeMinusL; j++) {

				int count = 0;

				for (int kk = 1, k = j + 1; kk < l; kk++, k++) {

					if ( (field[i][k - 1] - value) * (field[i][k] - value) < 0) {
						count++;
					}
				}

				if (count > maxNumOfSignChanges) {
					removeDirt(out, i, j, 1, l, value);
				}
			}
		}

		for (int j = 0; j < ySize; j++) {
			for (int i = 0; i < xSizeMinusL; i++) {
				int count = 0;

				for (int kk = 1, k = i + 1; kk < l; kk++, k++) {

					if (field[k - 1][j] * field[k][j] < 0) {
						count++;
					}
				}

				if (count > maxNumOfSignChanges) {
					removeDirt(out, i, j, l, 1, value);
				}
			}
		}
	}

}
