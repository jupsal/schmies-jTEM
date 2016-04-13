// JTEM - Java Tools for Experimental Mathematics
// Copyright (C) 2001 JEM-Group
//
// This program is free software; you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation; either version 2 of the License, or
// any later version.
//
// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with this program; if not, write to the Free Software
// Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307 USA

package de.jtem.riemann.schottky;


import java.beans.PropertyChangeListener;
import java.beans.PropertyChangeSupport;
import java.util.Vector;

import de.jtem.mfc.field.Complex;
import de.jtem.moebiusViewer.Attributes;
import de.jtem.moebiusViewer.MoebiusGraphics;
import de.jtem.moebiusViewer.shape.AbstractShape;

public class SchottkyShape extends AbstractShape {

	private static final long serialVersionUID = 1L;

	SchottkyData schottky;

	double[] alpha = new double[0];

	Complex C = new Complex();
	Complex D = new Complex();

	public SchottkyShape() {
	}

	public SchottkyShape(SchottkyData schottky) {
		setSchottky(schottky);
	}

	public SchottkyData getSchottky() {
		return schottky;
	}

	public void setSchottky(SchottkyData v) {
		schottky = v;
		if (alpha.length != schottky.getNumGenerators()) {
			alpha = new double[schottky.getNumGenerators()];

		}
		update();
	}

	protected PropertyChangeSupport propertyChangeSupport =
		new PropertyChangeSupport(this);

	private boolean showFixpoints = true;
	private boolean showCircles = true;
	private boolean showRadiusPoints = true;
	private boolean showInnerCircles = false;
	
	public void addPropertyChangeListener(PropertyChangeListener listener) {

		propertyChangeSupport.addPropertyChangeListener(listener);
	}

	public void removePropertyChangeListener(PropertyChangeListener listener) {
		propertyChangeSupport.removePropertyChangeListener(listener);
	}

	public void firePropertyChange(String propertyName) {
		propertyChangeSupport.firePropertyChange(propertyName, null, null);
	}

	public void update() {
	}

	public void draw(MoebiusGraphics context) {

		super.draw(context);

		context.setPointRadius(2);

		context.setVerticalTextLayout(Attributes.TOP);
		context.setHorizontalTextLayout(Attributes.RIGHT);

		for (int i = 0; i < schottky.numGenerators; i++) {

			draw(context, 0, i);
			draw(context, 1, i);

			if (showRadiusPoints) {

				computeCAndD(i, C, D);

				context.point(C.re, C.im);
				context.point(D.re, D.im);
			}

			if ( showInnerCircles && schottky instanceof Schottky) {

				for (int j = 0; j < 2; j++) {
					Complex[] innerCenter =
						((Schottky) schottky).innerCenter[i][j];
					double[] innerRadius =
						((Schottky) schottky).innerRadius[i][j];

					for (int k = 0; k < innerCenter.length; k++) {
						context.circle(
							innerCenter[k].re,
							innerCenter[k].im,
							innerRadius[k]);
					}
				}
			}
		}
	}

	void draw(MoebiusGraphics context, int j, int n) {

		if (showFixpoints) {

			Complex fp = schottky.fixpoint[n][j];
			context.point(fp.re, fp.im);

			if( getShowLabel()) { 
			context.text(fp.re, fp.im, (j == 0 ? "A(" : "B(") + (n + 1) + ") ");
			}
		}

		if (showCircles) {
			Complex center = schottky.center[n][j];
			context.circle(center.re, center.im, schottky.radius[n]);
		}
	}

	void compute(
		int i,
		boolean isPrime,
		Complex Z,
		double alpha,
		double radius) {

		Complex e = new Complex(Math.cos(alpha), Math.sin(alpha));
		Complex d =
			schottky.fixpoint[i][isPrime
				? 0
				: 1].minus(schottky.center[i][isPrime ? 0 : 1]);

		Complex kappa = d.times(e.conjugate()).plus(d.conjugate().times(e));

		double delta = d.absSqr() - radius * radius;

		double t = -kappa.re / 2 + Math.sqrt(kappa.re * kappa.re / 4 - delta);

		Z.assign(Math.cos(alpha) * t, Math.sin(alpha) * t);

		Z.assignPlus(schottky.fixpoint[i][isPrime ? 0 : 1]);
	}

	void computeCAndD(int i, Complex C, Complex D) {

		compute(i, true, C, alpha[i], schottky.radius[i]);

		//  C.assign( Math.cos( alpha[i] ) * schottky.radius[i],
		//  		  Math.sin( alpha[i] ) * schottky.radius[i] );

		//  	C.assignPlus( schottky.center[i][0] );

		schottky.generator[i].applyTo(C, D);

	}

	void computeD(Complex D) {

	}

	public Vector getTools() {

		Vector toolList = new Vector();

		toolList.addElement(new DragSchottkyTool());

		return toolList;
	}

	public boolean isShowFixpoints() {
		return showFixpoints;
	}

	public void setShowFixpoints(boolean showFixpoints) {
		this.showFixpoints = showFixpoints;
		this.firePropertyChange("showFixPoints");
	}

	public boolean isShowCircles() {
		return showCircles;
	}

	public void setShowCircles(boolean showCircles) {
		this.showCircles = showCircles;
		this.firePropertyChange("showCircles");
	}

	public boolean isShowRadiusPoints() {
		return showRadiusPoints;
	}

	public void setShowRadiusPoints(boolean showRadiusPoints) {
		this.showRadiusPoints = showRadiusPoints;
		this.firePropertyChange("showRadiusPoints");
	}
	
	public boolean isShowInnerCircles() {
		return showInnerCircles;
	}

	public void setShowInnerCircles(boolean showInnerCircles) {
		this.showInnerCircles = showInnerCircles;
		this.firePropertyChange("showInnerCircles");
	}
	

	class SchottkyCircleShape extends AbstractShape {

		int j, n;

		SchottkyCircleShape(int j, int n) {
			this.j = j;
			this.n = n;
		}

		public void draw(MoebiusGraphics context) {

			super.draw(context);

			context.setPointRadius(2);

			context.setVerticalTextLayout(Attributes.TOP);
			context.setHorizontalTextLayout(Attributes.RIGHT);

			SchottkyShape.this.draw(context, j, n);

			if (showRadiusPoints) {

				computeCAndD(n, C, D);

				if (j == 0)
					context.point(C.re, C.im);
				else
					context.point(D.re, D.im);
			}
		}
	}

	public AbstractShape getSchottkyCircleShape(int j, int m) {
		return new SchottkyCircleShape(j, m);
	}

}