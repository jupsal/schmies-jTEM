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
 * @author schmies
*/
public class PolygonOnTriangulation {

	final static double eps = 1e-12;

	IntList strip = new IntList();
	DoubleVectorList baryCoords = new DoubleVectorList(3);

	PolygonOnTriangulation() {
		super();
	}

	public void removeAll() {
		strip.removeAll();
		baryCoords.removeAll();
	}

	public void add(int element, double[] bary) {
		strip.add(element);
		baryCoords.add(bary);
	}

	public void assignGeodesic(
		int element,
		double[] bary,
		double[] direction,
		int[][] elementData,
		int[][] neighbourData,
		double[][] length) {

		double factor = 1;

		if (element < 0 || element >= elementData[0].length)
			throw new IllegalArgumentException("element not in element list of triangulation");

		if (!BaryCoords3d.isPoint(bary, eps)
			|| !BaryCoords3d.isInsideElement(bary, eps))
			throw new IllegalArgumentException("no valid barycentric coords");

		if (!BaryCoords3d.isVector(bary, eps))
			throw new IllegalArgumentException("no valid barycentric vector");

		removeAll();

		double[] newDirection = (double[]) direction.clone();
		direction = (double[]) direction.clone(); // do not want to change input

		add(element, (double[]) bary.clone());

		while (true) {

			double[] newBary = baryCoords.getNewLast();

			strip.add(element);

			double t = BaryCoords3d.add(bary, newBary, newDirection);

			if (1. - eps < t)
				break;

			int edge = BaryCoords3d.getEdge(newBary);

			if (edge != -1) {
				int newElement = neighbourData[element][edge];
				if (newElement == -1)
					throw new RuntimeException("hit boundary edge");
				BaryCoords3d.changeToOtherElementAtEdge(
					element,
					newElement,
					edge,
					direction,
					newDirection,
					elementData,
					length[element],
					length[newElement]);
				
				element=newElement;
				bary=newBary;
				double[] tmp=newDirection; newDirection=direction; direction=tmp;

			} else {

				int localVertex = BaryCoords3d.getVertex(bary);

				if (localVertex != -1)
					throw new RuntimeException("geodesic hit vertex; this case is not implemented yet");
				else
					throw new RuntimeException("fatal, error in code");
			}

			BaryCoords3d.vecScl(newDirection, newDirection, factor *= 1. - t);
		}
	}

	public double [][] getSpaceCurve( int [][] elementData, double [][] pointData ) {
		final int dim = pointData[0].length;
		final int nop = strip.size();
		
		final double [][] curveData = new double[nop][dim];
		
		for( int i=0; i<nop; i++) {
			
			final double [] bary = baryCoords.get(i);
			final double a = bary[0];
			final double b = bary[1];
			final double c = bary[2];
			
			final int [] element = elementData[strip.get(i)];
			final double []  A = pointData[ element[0] ];
			final double []  B = pointData[ element[1] ];
			final double []  C = pointData[ element[2] ];
			
			final double [] curve = curveData[i];
			
			for( int j=0; j<dim; j++ ) {
				curve[j] = a * A[j] + b * B[j] + c * C[j];
			}	
		}
		
		return curveData;
	}
	
	public double [] getSpaceCurve( int [][] elementData, double [] pointData, int dim ) {;
		final int nop = strip.size();
		
		final double [] curveData = new double[nop*dim];
		
		for( int i=0, J=0; i<nop; i++) {
			
			final double [] bary = baryCoords.get(i);
			final double a = bary[0];
			final double b = bary[1];
			final double c = bary[2];
			
			final int [] element = elementData[strip.get(i)];
			final int offsetA = element[0] * dim;
			final int offsetB = element[1] * dim;
			final int offsetC = element[2] * dim;
			
			for( int j=0; j<dim; j++, J++ ) {
				curveData[J] = a * pointData[ offsetA + j ] + b * pointData[ offsetB +j ] + c * pointData[ offsetC + j ];
			}	
		}
		
		return curveData;
	}
}
