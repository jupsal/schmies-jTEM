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

/*--- formatted by Jindent 2.1, (www.c-lab.de/~jindent) ---*/

package de.jtem.riemann.schottky;

import java.awt.Graphics;
import java.beans.Customizer;

import de.jtem.moebiusViewer.moebius.JMoebiusViewer;

public class SchottkyCustomizer extends JMoebiusViewer
    implements Customizer {

  private static final long     serialVersionUID = 1L;

  public SchottkyCustomizer () {}

  Object        theObject;
  boolean       firstPaint;

  public void setObject (Object o) {
    if (o != theObject) {
      theObject = o;

      removeAll();
      //add ( new CoordinateSystem() );
      add ( new SchottkyShape( (SchottkyData) o) );

      firstPaint = true;
    }
  }

  public void paint (Graphics g) {
    if (firstPaint) {
      encompass ();

      firstPaint = false;
    }

    super.paint (g);
  }

}

