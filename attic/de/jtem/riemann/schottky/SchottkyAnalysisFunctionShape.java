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
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with this program; if not, write to the Free Software
// Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA


package de.jtem.riemann.schottky;

import de.jtem.mfc.field.Complex;
import de.jtem.moebiusViewer.MoebiusGraphics;
import de.jtem.moebiusViewer.shape.AbstractShape;
import de.jtem.moebiusViewer.shape.ColorGradientBox;
import de.jtem.moebiusViewer.util.ColorGradient;

public class SchottkyAnalysisFunctionShape extends AbstractShape {

    public static final long serialVersionUID = 1L;

    SchottkyAnalysisFunction function;

    SchottkyShape schottkyShape = new SchottkyShape();

    double [] field, point;

    static final double logOf10 = Math.log( 10 );

    double pointDensity          = 0.1;
    double pointDensityOnCircles = 0.1;

    ColorGradient colorGradient = new ColorGradient();// Color.green, Color.red,  Color.yellow, Color.blue, 1 );

    double xDist = 1, yDist = 1;

    double durationOfLastUpdate;

    double minVal, maxVal;

    Complex min  = new Complex();
    Complex max  = new Complex();

    double[] bound ;

    boolean drawSchottkyShape    = true;
    boolean drawColorGradientBox = true;
    boolean computeExtrems       = true;
    boolean applyLog             = false;

    ColorGradientBox colorGradientBox = new ColorGradientBox( 0, 0 , 1 , 1, colorGradient );

    public SchottkyAnalysisFunctionShape() {
    }

    public SchottkyAnalysisFunctionShape( SchottkyAnalysisFunction function ) {
      this.function = function;
      update();
    }

    public SchottkyShape getSchottkyShape() {
      return schottkyShape;
    }

    public double getDurationOfLastUpdate() {
      return durationOfLastUpdate;
    }

    public boolean getApplyLog() {
        return applyLog;
    }

    public void setApplyLog(boolean  v) {
      if( applyLog == v )
        return;

      applyLog = v;

      firePropertyChange( "applyLog" );
    }

    public double [] getPoints() {
      return point;
    }

    public boolean getComputeExtrems() {
        return computeExtrems;
    }

    public void setComputeExtrems( boolean aBool ) {
        if( aBool == computeExtrems )
            return;

        computeExtrems = aBool;

        if( computeExtrems )
          computeExtrems();

        firePropertyChange( null );
    }

    public boolean getDrawSchottkyShape() {
        return drawSchottkyShape;
    }

    public void setDrawSchottkyShape( boolean aBool ) {
        drawSchottkyShape = aBool;

        firePropertyChange( "drawSchottkyShape" );
        fireMoebiusShapeChange();
    }

    public double getMinVal() {
        return minVal;
    }

    public void setMinVal( double aMin ) {
      minVal = aMin;
      setComputeExtrems( false );

      updateColorGradientBox();

      firePropertyChange( "minVal" );
      fireMoebiusShapeChange();
    }


    public double getMaxVal() {
        return maxVal;
    }

    public void setMaxVal( double aMax ) {
        maxVal = aMax;
        setComputeExtrems( false );

        updateColorGradientBox();

        firePropertyChange( "maxVal" );
        fireMoebiusShapeChange();
    }

    public Complex getMin() {
        return min;
    }

    public void setMin( Complex newMin ) {
      min.assign( newMin );

      minimize();
    }

    public Complex getMax() {
        return max;
    }

    public void setMax( Complex newMax ) {

        max.assign( newMax );

        maximize();
    }

    public double[] getBound() {
        return bound;
    }

    public SchottkyAnalysisFunction getSchottkyAnalysisFunction() {
        return function;
    }

    public void setSchottkyAnalysisFunction(SchottkyAnalysisFunction  v) {
        if( function == v )
            return;
        function = v;
        update();
    }

    public double getPointDensity() {
        return pointDensity;
    }

    public void setPointDensity(double  v) {
        pointDensity = v;
        update();
    }

    public double getPointDensityOnCircles() {
        return pointDensityOnCircles;
    }

    public void setPointDensityOnCircles(double  v) {
        pointDensityOnCircles = v;
        update();
    }

    public ColorGradient getColorGradient() {
        return colorGradient;
    }

    public void setColorGradient( ColorGradient aColorGradient ) {

        if( aColorGradient == null )
            throw new IllegalArgumentException("can not set color map to null" );

        colorGradient = aColorGradient;

        colorGradientBox.setColorMap( colorGradient );
        firePropertyChange( "colorGradient" );
        fireMoebiusShapeChange();
    }


    public double getXDist() {
        return xDist;
    }

    public void setXDist( double aDouble ) {
        xDist = aDouble;
        bound = getDefaultBound();
        update();
        firePropertyChange( "xDist" );
        fireMoebiusShapeChange();
    }

    public double getYDist() {
        return yDist;
    }

    public void setYDist( double aDouble ) {
        yDist = aDouble;
        bound = getDefaultBound();
        update();
        firePropertyChange( "yDist" );
        fireMoebiusShapeChange();
    }

    public double [] getDefaultBound() {

        double [] defaultBound = SchottkyDomainSampler.getBound( function.getSchottky() );

        defaultBound[0] -= xDist;
        defaultBound[1] -= yDist;
        defaultBound[2] += xDist;
        defaultBound[3] += yDist;

        return defaultBound;
    }

    public double [] getFocusedBound( int j, int n, double factor ) {

      double [] bound = new double[4];

      Complex center = function.getSchottky().center[n][j];
      double radius = function.getSchottky().radius[n];

      bound[0] = center.re - radius * factor;
      bound[1] = center.im - radius * factor;
      bound[2] = center.re + radius * factor;
      bound[3] = center.im + radius * factor;

      return bound;
    }

    public void setBound( double [] bound ) {
      this.bound = bound;
      update();
      firePropertyChange( "bound" );
      fireMoebiusShapeChange();
    }

    final double altLog( double a ) {

        if( applyLog )
            return Math.log( a ) / logOf10;
        else
            return a;
    }


    double evaluate( Complex z ) {
      return altLog( function.eval( z ) );
    }

    public void computeExtrems() {

        minVal = Double.MAX_VALUE;
        maxVal =-Double.MAX_VALUE;

        Complex z = new Complex();

        for( int i=0, k=0; k<point.length; i++, k+=2 ) {

            if( minVal > field[i] ) {
                minVal = field[i];
              z.assign( point[k], point[k+1] );
            min.assign( z );
          }

            if( maxVal < field[i] ) {
                maxVal = field[i];
              z.assign( point[k], point[k+1] );
            max.assign( z );
          }
        }

    }

    public void update() {

        if( function == null )
            return;

        schottkyShape.setSchottky( function.getSchottky() );

        if( bound == null )
            bound = getDefaultBound();

        point = SchottkyDomainSampler.getCover( function.getSchottky(),
                                          pointDensity, pointDensityOnCircles, bound );

        long timeMillis = System.currentTimeMillis();

        if( field == null || point.length != field.length * 2 )
            field = function.eval( point );
          else
            function.eval( point, field );

        durationOfLastUpdate = (System.currentTimeMillis()-timeMillis)/1000.;

        if( computeExtrems )
            computeExtrems();

          updateColorGradientBox();

        firePropertyChange(null);
        fireMoebiusShapeChange();
    }

    void updateColorGradientBox() {

      colorGradientBox.setX(bound[0]);
      colorGradientBox.setY(bound[1] - Math.abs(bound[1] - bound[3]) / 20);
      colorGradientBox.setW(bound[2] - bound[0]);
      colorGradientBox.setH(Math.abs(bound[1] - bound[3]) / 10);

      double minValue = altLog( minVal );
      double maxValue = altLog( maxVal );

      colorGradientBox.setShowLabel( true );

      if( (float)minValue - (int)minValue == 0 &&
          (float)maxValue - (int)maxValue == 0 ) {

        colorGradientBox.setLabel
            ("[ " + (int)minValue + " , " + (int)maxValue + " ]");
      } else {
        colorGradientBox.setLabel
            ("[ " + (float)minValue + " , " + (float)maxValue + " ]");
      }
    }

    public void draw( MoebiusGraphics context ) {

        super.draw( context );

        context.setPointOutline(0);

        final double minValue = altLog(minVal);
        final double maxValue = altLog(maxVal);

        for( int i=0, k=0; k<point.length; i++, k+=2 ) {
          final double t = ( altLog( field[i] ) - minValue ) / ( maxValue - minValue );
            context.setColor( colorGradient.getColor( t ) );
            context.point( point[k], point[k+1] );
        }

        if( drawSchottkyShape ) {
            schottkyShape.draw( context );
        }

        if( drawColorGradientBox )
          colorGradientBox.draw( context );
    }


    public void maximize() {
        maximize( 50, 0 );
    }
    public void maximizeStep() {
        maximize( 1, 1e-10 );
    }

    public void maximize( int n, double ftol ) {
      maxVal = function.maximize( max, ftol, n );

      updateColorGradientBox();

      firePropertyChange( null );
      fireMoebiusShapeChange();
    }

    public void minimize() {
        minimize( 50, 0 );
    }

    public void minimizeStep() {
        minimize( 1, 1e-10 );
    }

    public void minimize( int n, double ftol ) {
      minVal = function.minimize( min, ftol, n );

      updateColorGradientBox();

      firePropertyChange( null );
      fireMoebiusShapeChange();
    }

  public boolean isDrawColorGradientBox() {
    return drawColorGradientBox;
  }
  public void setDrawColorGradientBox(boolean drawColorGradientBox) {
    this.drawColorGradientBox = drawColorGradientBox;
  }

}
