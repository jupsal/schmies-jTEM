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

package de.jtem.numericalMethods.geometry.meshGeneration.ruppert;

import java.awt.Button;
import java.awt.Color;
import java.awt.Frame;
import java.awt.Graphics;
import java.awt.GridLayout;
import java.awt.Label;
import java.awt.Panel;
import java.awt.Scrollbar;
import java.awt.TextField;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
import java.awt.event.AdjustmentEvent;
import java.awt.event.AdjustmentListener;
import java.awt.event.WindowAdapter;
import java.awt.event.WindowEvent;

public class DelaunayViewer implements java.io.Serializable 
{
  private static final long serialVersionUID = 1L;

  double xmin,xmax,ymin,ymax;
  int FX_SIZE = 500;
  int FY_SIZE = 500;
  int X_SIZE = 500;
  int Y_SIZE;
  double[] point;
  int[] face;
  int nof;

  public DelaunayViewer(Ruppert d) { setTriangle(d); init(d); }

  public DelaunayViewer(int size,Ruppert d) { X_SIZE=size; setTriangle(d); init(d); }

  public void setTriangle(Ruppert d) { 
    point=d.getPoints(); face=d.getIndices(); nof=face.length/3; 
    int l=point.length; xmin=xmax=point[0]; ymin=ymax=point[1]; 
    for (int i=2;i<l;) {
      double x=point[i++]; if (x<xmin) xmin=x; if (x>xmax) xmax=x;
      double y=point[i++]; if (y<ymin) ymin=y; if (y>ymax) ymax=y;
    }
    Y_SIZE = (int)(((double)X_SIZE)*(ymax-ymin)/(xmax-xmin));
  }

  final void init(Ruppert d) {
    Frame f = new Frame("Triangulation");
    f.setSize(FX_SIZE,FY_SIZE);
    DrawingPanel p = new DrawingPanel();
    f.add(p);
    f.setVisible(true);
    f.addWindowListener(new WindowAdapter() { public void windowClosing(WindowEvent e) { System.exit(0); } });

    Frame f2 = new Frame("Ruppertment constraints");
    DelaunayPanel p2=new DelaunayPanel(d,this,p);
    f2.add(p2);
    f2.setSize(X_SIZE,X_SIZE/3);
    f2.setVisible(true);
  }
  
  int[] xx=new int[3], yy=new int[3];
  class DrawingPanel extends Panel {
    public void paint(Graphics g) {
      g.setColor(Color.lightGray); g.fillRect(0,0,X_SIZE,Y_SIZE);
      double xdif=(4*X_SIZE)/(5*(xmax-xmin)),ydif=(4*Y_SIZE)/(5*(ymax-ymin));
      int x_border=X_SIZE/10, y_border=Y_SIZE/10;
      for (int i=0,fp=0;i<nof;i++) {
	for (int j=0;j<3;j++,fp++) {
	  xx[j]=x_border+(int)((point[2*face[fp]]-xmin)*xdif); 
	  yy[j]=Y_SIZE-y_border-(int)((point[2*face[fp]+1]-ymin)*ydif);
	}
	g.setColor(Color.green); 
	g.fillPolygon(xx,yy,3); 
	g.setColor(Color.blue); 
	g.drawPolygon(xx,yy,3); 
      }
    }
  }

  class DelaunayPanel extends Panel {

    int area,angle;
    Scrollbar angleS,areaS;
    TextField angleT,areaT;
    Button quitB,refineB;
    
    DelaunayPanel(final Ruppert d,final DelaunayViewer dv,final DrawingPanel dp) {
      angle = (int)(d.getAngleConstraint());
      area  = (int)(100*d.getAreaConstraint());
      setLayout(new GridLayout(3,3));
      angleS = new Scrollbar(Scrollbar.HORIZONTAL, 0, 1, 0, 35);
      angleS.addAdjustmentListener(new AdjustmentListener() { 
	public void adjustmentValueChanged(AdjustmentEvent e) { angle=e.getValue(); repaint(); }});
      angleT = new TextField();
      angleT.addActionListener(new ActionListener() { 
	public void actionPerformed(ActionEvent e) { 
	  angle=Integer.valueOf(angleT.getText()).intValue(); repaint(); }});
      areaS = new Scrollbar(Scrollbar.HORIZONTAL, 0, 10, 0, 1000);
      areaS.addAdjustmentListener(new AdjustmentListener() { 
	public void adjustmentValueChanged(AdjustmentEvent e) { area=1+e.getValue(); repaint(); }});
      areaT = new TextField();
      areaT.addActionListener(new ActionListener() { 
	public void actionPerformed(ActionEvent e) { 
	  area=(int)(100*Double.valueOf(areaT.getText()).doubleValue()); repaint(); }});
      quitB = new Button("quit");
      quitB.addActionListener(new ActionListener() { 
	public void actionPerformed(ActionEvent e) { System.exit(0); }});
      refineB = new Button("refine");
      refineB.addActionListener(new ActionListener() { 
	public void actionPerformed(ActionEvent e) { 
	  d.setAngleConstraint(angle); d.setAreaConstraint(0.01*(double)area); d.refine();
	  dv.setTriangle(d); dp.repaint(); }});
      add(new Label ("angle")); add(angleS); add(angleT);
      add(new Label ("area ")); add(areaS);  add(areaT);
      add(new Label (""));      add(quitB);  add(refineB); 
    }

    public void paint(Graphics g) {
      if (angle<0) angle=0; if (angle>34) angle=34; if (area<1) area=1; if (area>1000) area=1000; 
      angleT.setText(Integer.toString(angle)); angleS.setValue(angle);
      areaT.setText(Double.toString(0.01*(double)area)); areaS.setValue(area-1);
    }
  }


  public static void main(String[] arg) 
  {
    double[] h={-2,0,1, 2,0,1};
    int[] ref={20,20};
    new DelaunayViewer(Initializer.init(-4,-2,4,2,2,2,h,ref));
  }

}
