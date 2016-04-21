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
 * Numerical utility function dealing with barycentric coordinates. This class
 * is a port from c code source taken form the objective-C oorange.
 * 
 * @author schmies
 */
public class BaryCoords3d {

	static double eps = 1.e-6;
	static int FALSE = 0; //TODO: remove
	static int TRUE = 0; //TODO: remove

	BaryCoords3d() {
	}

	static void correct(double[] bary) {

		double sum = 0;
		int j, k, i = 0;

		for (i = 0; i < 3; i++) {
			if (Math.abs(bary[i]) < eps)
				bary[i] = 0;
			sum += bary[i];
		}

		for (i = 0; bary[i] == 0. && i < 3; i++);

		bary[(j = (i + 1) % 3)] /= sum;
		bary[(k = (i + 2) % 3)] /= sum;
		bary[i] = 1. - bary[j] - bary[k];
	}

	public static double sqr(
		double a,
		double b,
		double c,
		double[] v) {
		return (
			- (
				a * a * v[1] * v[2]
					+ b * b * v[0] * v[2]
					+ c * c * v[1] * v[0]));
	}

	public static double abs(
		double a,
		double b,
		double c,
		double[] v) {
		return (
			Math.sqrt(
				- (
					a * a * v[1] * v[2]
						+ b * b * v[0] * v[2]
						+ c * c * v[1] * v[0])));
	}

	public static double dot(
		double a,
		double b,
		double c,
		double[] bV,
		double[] bW) {
		return (
			-.5
				* (a * a * (bV[1] * bW[2] + bV[2] * bW[1])
					+ b * b * (bV[0] * bW[2] + bV[2] * bW[0])
					+ c * c * (bV[1] * bW[0] + bV[0] * bW[1])));
	}

	public static double dist(
		double a,
		double b,
		double c,
		double[] A,
		double[] B) {
		double dist0 = A[0] - B[0];
		double dist1 = A[1] - B[1];
		double dist2 = A[2] - B[2];

		return (
			Math.sqrt(
				- (
					a * a * dist1 * dist2
						+ b * b * dist0 * dist2
						+ c * c * dist1 * dist0)));
	}

	static void norm(
		double a,
		double b,
		double c,
		double[] baryNorm,
		double[] bary) {

		vecScl(
			baryNorm,
			bary,
			1.
				/ Math.sqrt(
					- (
						a * a * bary[1] * bary[2]
							+ b * b * bary[0] * bary[2]
							+ c * c * bary[1] * bary[0])));
	}

	static void rotate90(
		double a,
		double b,
		double c,
		double[] normal,
		double[] bary) {

		double sp = sqr(a, b, c, bary);

		normal[0] = bary[2];
		normal[1] = bary[0];
		normal[2] = bary[1];

		vecLinComb(
			normal,
			normal,
			1.,
			bary,
			-dot(a, b, c, bary, normal) / sp);

		vecScl(
			normal,
			normal,
			Math.sqrt(sp) / abs(a, b, c, normal));
	}

	static void normal(
		double a,
		double b,
		double c,
		double[] normal,
		double[] bary) {
		normal[0] = bary[2];
		normal[1] = bary[0];
		normal[2] = bary[1];

		vecLinComb(
			normal,
			normal,
			1.,
			bary,
			-dot(a, b, c, bary, normal)
				/ sqr(a, b, c, bary));

		vecScl(normal, normal, 1. / abs(a, b, c, normal));
	}

	static void orthonormalize(
		double a,
		double b,
		double c,
		double[] bary,
		double[] normal) {

		norm(a, b, c, bary, bary);

		vecLinComb(
			normal,
			normal,
			1.,
			bary,
			-dot(a, b, c, bary, normal));

		norm(a, b, c, normal, normal);
	}

	public static double angle(
		double a,
		double b,
		double c,
		double[] bV,
		double[] bW) {
		double lV =
			Math.sqrt(
				- (
					a * a * bV[1] * bV[2]
						+ b * b * bV[0] * bV[2]
						+ c * c * bV[1] * bV[0]));
		double lW =
			Math.sqrt(
				- (
					a * a * bW[1] * bW[2]
						+ b * b * bW[0] * bW[2]
						+ c * c * bW[1] * bW[0]));

		double sp = dot(a, b, c, bV, bW) / lV / lW;

		sp = Math.min(1., Math.max(-1., sp));

		return Math.acos(sp);
	}

	public static double orientedAngle(
		double a,
		double b,
		double c,
		double[] bV,
		double[] bW) {
		double[] nV = new double[3];

		double lV =
			Math.sqrt(
				- (
					a * a * bV[1] * bV[2]
						+ b * b * bV[0] * bV[2]
						+ c * c * bV[1] * bV[0]));
		double lW =
			Math.sqrt(
				- (
					a * a * bW[1] * bW[2]
						+ b * b * bW[0] * bW[2]
						+ c * c * bW[1] * bW[0]));

		double sp = dot(a, b, c, bV, bW) / lV / lW;

		sp = Math.min(1., Math.max(-1., sp));

		normal(a, b, c, nV, bV);

		if (dot(a, b, c, nV, bW) < 0)
			return -Math.acos(sp);
		else
			return Math.acos(sp);
	}

	static void rotate(
		double a,
		double b,
		double c,
		double[] bary,
		double[] rotated,
		double alpha) {
		double[] normal = new double[3];

		normal(a, b, c, normal, bary);

		vecScl(normal, normal, abs(a, b, c, bary));

		vecLinComb(
			rotated,
			bary,
			Math.cos(alpha),
			normal,
			Math.sin(alpha));
	}

	/*
	 * calc point x with for barycentric coordinates (el, bary) on t2
	 */
	static void convertToVec3(
		double[] x,
		double[] bary,
		double[] A,
		double[] B,
		double[] C) {
		x[0] = bary[0] * A[0] + bary[1] * B[0] + bary[2] * C[0];
		x[1] = bary[0] * A[1] + bary[1] * B[1] + bary[2] * C[1];
		x[2] = bary[0] * A[2] + bary[1] * B[2] + bary[2] * C[2];
	}

	/*
	 * calc barycentric coordinates bary for point x in anElement el, not
	 * necessarily 0 <= b[i] <= 1
	 */
	int gmBaryInTriangleConvertFromVec3(
		double[] bary,
		double[] x0,
		double[] x1,
		double[] x2,
		double[] x) {
		//TODO: use exception handling for degenerate situations
		//TODO: the return type should be void or double which should give
		// distance from plane

		int i0 = 0, i1 = 1, i2 = 2;

		double det;

		/* find two linear independent rows */
		for (;;) {
			det =
				x1[i0] * x2[i1]
					- x1[i1] * x2[i0]
					- (x0[i0] * x2[i1] - x0[i1] * x2[i0])
					+ x0[i0] * x1[i1]
					- x0[i1] * x1[i0];

			if (Math.abs(det) > eps)
				break;

			if (i1 == 1) {
				i1 = 2;
				i2 = 1;
			} else if (i0 == 0) {
				i0 = 1;
				i2 = 0;
			} else {

				System.out.println(
					"mBaryInElementConvertFromVec3: triangle degenerated?");
				//				
				//				fprintf(stderr, "gmBaryInElementConvertFromVec3: triangle
				// degenerated?\n");
				//				fprintf(stderr, "\tcan't compute barycentric
				// coordinates.\n");
				//				fprintf(stderr, "vertex0: (%f, %f, %f)\n", x0[0], x0[1],
				// x0[2]);
				//				fprintf(stderr, "vertex1: (%f, %f, %f)\n", x1[0], x1[1],
				// x1[2]);
				//				fprintf(stderr, "vertex2: (%f, %f, %f)\n", x2[0], x2[1],
				// x2[2]);
				//				fprintf(stderr, "point : (%f, %f, %f)\n", x[0], x[1], x[2]);
				//				fflush(stderr);

				return FALSE;
			}
		}

		/* calculate barycentric coordinates */
		bary[0] =
			(x1[i0] * x2[i1]
				- x1[i1] * x2[i0]
				- (x[i0] * x2[i1] - x[i1] * x2[i0])
				+ x[i0] * x1[i1]
				- x[i1] * x1[i0])
				/ det;
		bary[1] =
			(x[i0] * x2[i1]
				- x[i1] * x2[i0]
				- (x0[i0] * x2[i1] - x0[i1] * x2[i0])
				+ x0[i0] * x[i1]
				- x0[i1] * x[i0])
				/ det;
		bary[2] = 1.0 - bary[0] - bary[1];

		/* test third row */
		if (Math
			.abs(x0[i2] * bary[0] + x1[i2] * bary[1] + x2[i2] * bary[2] - x[i2])
			> 1.e-3) {
			System.out.println(
				"gmBaryInElementConvertFromVec3: test for third row failed.");
			//			fprintf(stderr, "gmBaryInElementConvertFromVec3: test for third
			// row failed.\n");
			//			fprintf(stderr, "\tpoint not in triangle plane?\n");
			//			fflush(stderr);

			return FALSE;
		}

		correct(bary);

		return TRUE;
	}

	/*
	 * clip the line from bary to nbary at the edge of the triangle, bary is
	 * expected to be inside the triangle, nbary is clipped if necessary, the
	 * function returns if a clipping was performed
	 */
	int clip(double[] bary, double[] nbary) {
		double[] v = new double[3];
		double fac;
		int i;
		int clipped = FALSE;
		double fac_min = 1.e98;
		/* minimal scaling factor for v to reach first edge */

		for (i = 0; i < 3; i++)
			if (!((0.0 <= nbary[i]) && (nbary[i] <= 1.0)) && nbary[i] < 0.) {
				fac = bary[i] / (bary[i] - nbary[i]);
				if (fac < fac_min)
					fac_min = fac;
			}

		if (fac_min != 1.e98) {
			clipped = TRUE;

			DoubleTriple.sub(v, nbary, bary);
			DoubleTriple.scale(v, v, fac_min);
			DoubleTriple.add(nbary, bary, v);
			correct(nbary);
		}

		return clipped;
	}

	static void interpolate(
		double[] bary,
		double[] position,
		double[][] triple) {
		bary[0] =
			position[0] * triple[0][0]
				+ position[1] * triple[1][0]
				+ position[2] * triple[2][0];
		bary[1] =
			position[0] * triple[0][1]
				+ position[1] * triple[1][1]
				+ position[2] * triple[2][1];
		bary[2] =
			position[0] * triple[0][2]
				+ position[1] * triple[1][2]
				+ position[2] * triple[2][2];
	}

	static void interpolateVec3(
		double[] out,
		double[] bary,
		double[] A,
		double[] B,
		double[] C) {
		out[0] = bary[0] * A[0] + bary[1] * B[0] + bary[2] * C[0];
		out[1] = bary[0] * A[1] + bary[1] * B[1] + bary[2] * C[1];
		out[2] = bary[0] * A[2] + bary[1] * B[2] + bary[2] * C[2];

	}

	static void interpolateVec2f(
		float[] out,
		double[] bary,
		float[] A,
		float[] B,
		float[] C) {
		throw new RuntimeException("check behavior of this function. ");
		//		out[0] = (float)(bary[0] * A[0] + bary[1] * B[0] + bary[2] * C[0]);
		//		out[1] = (float)(bary[0] * A[1] + bary[1] * B[1] + bary[2] * C[1]);
		//		out[2] = (float)(bary[0] * A[2] + bary[1] * B[2] + bary[2] * C[2]);

	}

	static double add(double[] new_pos, double[] old_pos, double[] dir) {
		//TODO: check behavior of function; add exception handling; question
		// return type
		double[] h = new double[3];
		double t, sum, abs_sum;
		int i, j, k;

		/*
		 * gmBaryCorrect( old_pos );
		 */

		for (abs_sum = 0, sum = 0, i = 0; i < 3; i++) {
			sum += dir[i];
			abs_sum += Math.abs(dir[i]);
		}

		if (Math.abs(sum) > abs_sum * eps) {
			System.out.println("not a barycentric vector");
			//			printf("gmBaryAdd:");
			//			printf("not a barycentric vector\n" );
			//
			//			printf(" %f %f %f\n", dir[0], dir[1], dir[2] );
			return (-1);
		}

		for (i = 0; i < 3; i++)
			if (Math.abs(dir[i]) < abs_sum * eps) {
				dir[(i + 1) % 3] += dir[i] / 2.;
				dir[(i + 2) % 3] = -dir[(i + 1) % 3];
				dir[i] = 0;
			}

		addVec(h, old_pos, dir);

		if (isInsideElement(h)) {
			DoubleTriple.copy(new_pos, h);
			return (1.);
		}

		for (i = 0; i < 3; i++)
			if (dir[i] != 0)
				if ((t = -old_pos[i] / dir[i]) > 0
					&& (h[(j = (i + 1) % 3)] = old_pos[j] + t * dir[j]) >= 0
					&& h[j] <= 1
					&& (h[(k = (i + 2) % 3)] = old_pos[k] + t * dir[k]) >= 0
					&& h[k] <= 1) {

					new_pos[i] = 0;
					new_pos[j] = h[j];
					new_pos[k] = 1. - h[j];

					return (t);
				}

		/*
		 * in this situation at least one of the barycentric coord. is zero and
		 * the direction vector points out of the triangle.
		 */

		for (i = 0; i < 3; i++)
			if (dir[i] != 0)
				if ((t = -old_pos[i] / dir[i]) > eps
					&& (h[(j = (i + 1) % 3)] = old_pos[j] + t * dir[j]) > -eps
					&& h[j] < 1 + eps
					&& (h[(k = (i + 2) % 3)] = old_pos[k] + t * dir[k]) > -eps
					&& h[k] < 1 + eps) {

					new_pos[i] = 0;
					new_pos[j] = Math.min(Math.max(h[j], 0), 1);
					new_pos[k] = 1. - h[j];

					return (t);
				}

		/* Last attempt failed again, so exit with doing nothing */

		DoubleTriple.copy(new_pos, old_pos);

		return (0);
	}

	/*
	 * test for barycentric coordinates
	 */

	public static int getVertex(double[] bary) {
		if (bary[0] == 1. && bary[1] == 0. && bary[2] == 0.)
			return (0);
		if (bary[0] == 0. && bary[1] == 1. && bary[2] == 0.)
			return (1);
		if (bary[0] == 0. && bary[1] == 0. && bary[2] == 1.)
			return (2);

		return (-1);
	}

	public static int getEdge(double[] bary) {
		if (bary[0] == 0. && bary[1] != 0. && bary[2] != 0.)
			return (0);
		if (bary[0] != 0. && bary[1] == 0. && bary[2] != 0.)
			return (1);
		if (bary[0] != 0. && bary[1] != 0. && bary[2] == 0.)
			return (2);

		return (-1);
	}

	/*
	 * Functions for handling with barycentric coordinats, they sometimes need
	 * element or/and neighbour information given as arrays of int [].
	 */

	public static int onEdgeSwitchElement(
		double[] bary,
		int[] anElement,
		int[][] elementData,
		int[][] neighbourData) {
		//TODO: seperate this into 2 steps: 1st getting the neighbour element;
		// 2nd recompute the bary coords
		double[] temp = new double[3];
		int edge = getEdge(bary);
		int T = neighbourData[anElement[0]][edge];

		if (T == -1);

		temp[IntTriple.getLocInd(
			elementData[T],
			elementData[anElement[0]][(edge + 1) % 3])] =
			bary[(edge + 1) % 3];
		temp[IntTriple.getLocInd(
			elementData[T],
			elementData[anElement[0]][(edge + 2) % 3])] =
			bary[(edge + 2) % 3];

		DoubleTriple.copy(bary, temp);

		anElement[0] = T;

		return (1);
	}

	public static void outsidePointingEdgeNormal( double [] l, int element, int edge, double [] normal  ) {
		double [] n = normal;
		double [] v = new double[3];
		
		n[ edge ] = -1;
		n[ (edge + 1) % 3 ] = 0.5;
		n[ (edge + 1) % 3 ] = 0.5;
		
		v[ (edge + 1) % 3 ] = -1/l[edge];
		v[ (edge + 2) % 3 ] =  1/l[edge];
		
		vecLinComb( n, n, 1, v, -dot(l[0], l[1], l[2], n, v) );
		
		vecScl(n, n, 1. / abs(l[0], l[1], l[2],n));
	}
	
	public static void changeToOtherElementAtEdge(
			int element1, int element2, int edge1,
			double[] vec1,
			double[] vec2,
			int[][] elementData,
			double [] l1, double [] l2 ) {
		
		double [] n = new double[3];
		double [] v = new double[3];
		
		int localIndexOfV1InElement1 = (edge1 + 1) % 3;
		int localIndexOfV2InElement1 = (edge1 + 2) % 3;
		
		int V1 = elementData[element1][localIndexOfV1InElement1];
		int V2 = elementData[element1][localIndexOfV2InElement1];
		
		v[ localIndexOfV1InElement1 ] = -1 / l1[edge1];
		v[ localIndexOfV2InElement1 ] =  1 / l1[edge1];
		
		outsidePointingEdgeNormal( l1, element1, edge1, n);
		
		double rCosAlpha = dot( l1[0], l1[1], l1[2], n, vec1 );
		double rSinAlpha = dot( l1[0], l1[1], l1[2], v, vec1 );
	
		int localIndexOfV1InElement2 = IntTriple.getLocInd( elementData[element2], V1 );
		int localIndexOfV2InElement2 = IntTriple.getLocInd( elementData[element2], V2 );
		
		int edge2 = 2 *(localIndexOfV1InElement2+localIndexOfV2InElement2) % 3;
		
		v[0]=v[1]=v[2] = 0;
		v[ localIndexOfV1InElement2 ] = -1 / l2[edge2];
		v[ localIndexOfV2InElement2 ] =  1 / l2[edge2];
		
		outsidePointingEdgeNormal( l2, element2, edge2, n);
		
		vecLinComb( vec2, n, -rCosAlpha, v, rSinAlpha );
	}
	
	/**
	 * @deprecated
	 * @param bary
	 * @param anElement
	 * @param elementData
	 * @param neighbourData
	 * @return
	 */
	public static int getOtherElementAtEdge(
		double[] bary,
		int anElement,
		int[][] elementData,
		int[][] neighbourData) {
		int i;

		for (i = 0; i < 3; i++)
			if (bary[i] == 0.)
				break;

		if (i == 3) {
			throw new IllegalArgumentException("bary0d is not on an edge.");
		}

		return neighbourData[anElement][i];
	}

	public static int setNewElement(
		double[] bary,
		int[] anElement,
		int[][] elementData,
		int[][] neighbourData,
		int newElement) {
		//TODO: check usage; add exeption handle; question return type
		int local;

		if (anElement[0] == newElement)
			return (1);

		if (isInnerPoint(bary))
			return (newElement == anElement[0] ? TRUE : FALSE);

		if (isOnEdge(bary))
			if (onEdgeSwitchElement(bary,
				anElement,
				elementData,
				neighbourData)
				== TRUE)
				if (anElement[0] == newElement)
					return (1);
				else {
					onEdgeSwitchElement(
						bary,
						anElement,
						elementData,
						neighbourData);
					return (0);
				}
			else
				return (0);

		else {

			local =
				IntTriple.getLocInd(
					elementData[newElement],
					elementData[anElement[0]][getVertex(bary)]);
			if (local == -1)
				return (0);

			anElement[0] = newElement;

			bary[local] = 1.;
			bary[(local + 1) % 3] = bary[(local + 2) % 3] = 0.;

			return (1);
		}
	}

	/** @deprecated */
	public static int getVerticesOfEdge(
		double[] bary,
		int anElement,
		int[][] elementData,
		int[] left,
		int[] right) {
		int edge = getEdge(bary);

		if (edge == -1)
			return (0);

		left[0] = elementData[anElement][(edge + 2) % 3];
		right[0] = elementData[anElement][(edge + 1) % 3];

		return (1);
	}

	public static boolean vertexIsOnEdge(
		double[] bary,
		int anElement,
		int[][] elementData,
		int aVertex) {
		int edge = getVertex(bary);

		if (edge == -1)
			return (false);

		return (
			elementData[anElement][(edge + 2) % 3] == aVertex
				|| elementData[anElement][(edge + 1) % 3] == aVertex);
	}

	public static int leftVertexOfEdge(
		double[] bary,
		int anElement,
		int[][] elementData) {
		int i;

		for (i = 0; i < 3; i++)
			if (bary[i] == 0.)
				break;

		if (i == 3)
			throw new RuntimeException("bary0d is not on an edge");

		return elementData[anElement][(i + 2) % 3];
	}

	public static int rightVertexOfEdge(
		double[] bary,
		int anElement,
		int[][] elementData) {
		int i;

		for (i = 0; i < 3; i++)
			if (bary[i] == 0.)
				break;

		if (i == 3)
			throw new RuntimeException("bary0d is not on an edge.");

		return (elementData[anElement][(i + 1) % 3]);
	}

	public static int changeBase(
		double[] bary,
		int[] oldBase,
		int[] newBase) {
		int i, j;
		double[] tmp = new double[3];

		DoubleTriple.copy(tmp, bary);
		DoubleTriple.zero(bary);

		for (i = 0; i < 3; i++)
			if (tmp[i] != 0.) {
				for (j = 0; j < 3; j++)
					if (oldBase[i] == newBase[j]) {
						bary[j] = tmp[i];
						break;
					}
				if (j == 3) {
					DoubleTriple.copy(bary, tmp);
					return (0);
				}
			}

		return 1;
	}

	public static int canChangeBase(
		double[] bary,
		int[] oldBase,
		int[] newBase) {
		int i, j;

		for (i = 0; i < 3; i++)
			if (bary[i] != 0.) {
				for (j = 0; j < 3; j++)
					if (oldBase[i] == newBase[j]) {
						break;
					}
				if (j == 3) {
					return (0);
				}
			}

		return 1;
	}

	public static void vecScl(
		double[] baryVec1,
		double[] baryVec2,
		double f) {
		double factor1 = (f);
		(baryVec1)[0] = (baryVec2)[0] * factor1;
		(baryVec1)[1] = (baryVec2)[1] * factor1;
		(baryVec1)[2] = - (baryVec1)[0] - (baryVec1)[1];
	}

	public static void vecLinComb(
		double[] baryVec0,
		double[] baryVec1,
		double f1,
		double[] baryVec2,
		double f2) {
		double factor1 = (f1);
		double factor2 = (f2);
		(baryVec0)[0] = (baryVec1)[0] * factor1 + (baryVec2)[0] * factor2;
		(baryVec0)[1] = (baryVec1)[1] * factor1 + (baryVec2)[1] * factor2;
		(baryVec0)[2] = - (baryVec0)[0] - (baryVec0)[1];
	}

	public static void linComb(
		double[] bary0,
		double[] bary1,
		double f1,
		double[] bary2,
		double f2) {
		double factor1 = (f1);
		double factor2 = (f2);
		(bary0)[0] = (bary1)[0] * factor1 + (bary2)[0] * factor2;
		(bary0)[1] = (bary1)[1] * factor1 + (bary2)[1] * factor2;
		(bary0)[2] = 1. - (bary0)[0] - (bary0)[1];
	}

	public static void addVec(
		double[] bary1,
		double[] bary2,
		double[] baryVec) {
		(bary1)[0] = (bary2)[0] + (baryVec)[0];
		(bary1)[1] = (bary2)[1] + (baryVec)[1];
		(bary1)[2] = 1. - (bary1)[0] - (bary1)[1];
	}

	public static void subVec(
		double[] bary1,
		double[] bary2,
		double[] baryVec) {
		(bary1)[0] = (bary2)[0] - (baryVec)[0];
		(bary1)[1] = (bary2)[1] - (baryVec)[1];
		(bary1)[2] = 1. - (bary1)[0] - (bary1)[1];
	}

	public static void sub(
		double[] baryVec,
		double[] bary1,
		double[] bary2) {
		(baryVec)[0] = (bary1)[0] - (bary2)[0];
		(baryVec)[1] = (bary1)[1] - (bary2)[1];
		(baryVec)[2] = - (baryVec)[0] - (baryVec)[1];
	}

	public static void vecAdd(
		double[] baryVec0,
		double[] baryVec1,
		double[] baryVec2) {
		(baryVec0)[0] = (baryVec1)[0] + (baryVec2)[0];
		(baryVec0)[1] = (baryVec1)[1] + (baryVec2)[1];
		(baryVec0)[2] = - (baryVec0)[0] - (baryVec0)[1];
	}

	public static void vecSub(
		double[] baryVec0,
		double[] baryVec1,
		double[] baryVec2) {
		(baryVec0)[0] = (baryVec2)[0] - (baryVec1)[0];
		(baryVec0)[1] = (baryVec2)[1] - (baryVec1)[1];
		(baryVec0)[2] = - (baryVec0)[0] - (baryVec0)[1];
	}

	public static boolean isPoint( double [] bary ) {
		return isPoint( bary, eps );
	}
	
	public static boolean isPoint( double [] bary, double eps ) {
		return Math.abs( bary[0] + bary[1] + bary[2] - 1 ) < eps;
	}
	
	public static boolean isVector( double [] bary ) {
		return isVector( bary, eps );
	}
	
	public static boolean isVector( double [] bary, double eps ) {
		return Math.abs( bary[0] + bary[1] + bary[2] ) < eps;
	}
	
	public static boolean isInsideElement(double[] bary ) {
		return isInsideElement( bary, eps );
	}

	public static boolean isInsideElement(double[] bary, double eps) {
		return (
			(bary)[0] > - (eps)
				&& (bary)[0] < 1 + eps
				&& (bary)[1] > - (eps)
				&& (bary)[1] < 1 + eps
				&& (bary)[2] > - (eps)
				&& (bary)[2] < 1 + eps);
	}

	public static boolean isInnerPoint(double[] bary) {
		return (
			(bary)[0] > 0.
				&& (bary)[0] < 1.
				&& (bary)[1] > 0.
				&& (bary)[1] < 1.
				&& (bary)[2] > 0.
				&& (bary)[2] < 1.);
	}

	public static boolean isOnEdge(double[] bary) {
		return (
			((bary)[0] == 0. && (bary)[1] != 0. && (bary)[2] != 0.)
				|| ((bary)[0] != 0. && (bary)[1] == 0. && (bary)[2] != 0.)
				|| ((bary)[0] != 0. && (bary)[1] != 0. && (bary)[2] == 0.));
	}

	public static boolean isOnVertex(double[] bary) {
		return (
			((bary)[0] == 1. && (bary)[1] == 0. && (bary)[2] == 0.)
				|| ((bary)[0] == 0. && (bary)[1] == 1. && (bary)[2] == 0.)
				|| ((bary)[0] == 0. && (bary)[1] == 0. && (bary)[2] == 1.));
	}

}
