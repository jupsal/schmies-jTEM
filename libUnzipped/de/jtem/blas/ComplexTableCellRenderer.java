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

import java.awt.Color;
import java.awt.Component;
import java.awt.FontMetrics;
import java.text.DecimalFormat;
import java.util.StringTokenizer;

import javax.swing.BorderFactory;
import javax.swing.JLabel;
import javax.swing.JTable;
import javax.swing.SwingUtilities;
import javax.swing.border.Border;
import javax.swing.table.TableCellRenderer;

import de.jtem.mfc.field.ComplexValue;
import de.jtem.mfc.field.Field;

/**
 * This special TableCellRenderer made for ComplexMatrix and ComplexVector
 * changes especially those values, which are too long to fit in the cell.
 * The decimal places of both parts (the real part and the imaginary part) the
 * complex number consits of were cut as much as necessary but as few as
 * possible. Trunks and decimal powers are even displayed, to give the user an
 * imagine of the range the value has. To show the user that a value was
 * changed in its precision, the text color changes to green.
 */
public class ComplexTableCellRenderer extends JLabel
  implements TableCellRenderer
{
  private Color lightBlue, darkBlue, darkGreen, lightGreen;
  Border b;
  String complexStr;
  Field.Complex complexVal;
  StringBuffer buffer;
  StringTokenizer tokenizer;
  FontMetrics fontMetrics;
  int columnWidth;
  int stringWidth;
  boolean isCut=false;
  boolean noMoreDigs=false;
  /**
   * Creates a new instance of a ComplexTableCellRenderer.
   */
  public ComplexTableCellRenderer()
  {
    lightBlue = new Color(160, 160, 255);
    darkBlue = new Color( 64,  64, 128);
    lightGreen = new Color(192, 255, 192);
    darkGreen = new Color( 32,  96,  32);
    b = BorderFactory.createEmptyBorder(1, 1, 1, 1);
    setOpaque(true);
    setBorder(b);
  }
  /**
   * Returns the Component which renders the cells of the specified table.
   * @return <code>this</code>
   */
  public Component getTableCellRendererComponent(JTable table,
    Object value, boolean isSelected, boolean hasFocus, int row, int column)
  {
    isCut=false;
    noMoreDigs=false;
    complexVal = (Field.Complex)value;
    DecimalFormat format = new DecimalFormat("0.");
    String re=format.format(complexVal.getRe());
    String im=format.format(complexVal.getIm());
    complexStr = ComplexValue.toString((Field.Complex)value, false);
    fontMetrics = table.getFontMetrics(table.getFont());
    stringWidth = SwingUtilities.computeStringWidth(fontMetrics, complexStr);
    columnWidth = table.getColumnModel().getColumn(column).getWidth()-2;
    setFont(table.getFont());
    setForeground(table.getForeground());
    setBackground(table.getBackground());
    if(stringWidth > columnWidth)
    {
      setForeground(darkGreen);
      isCut=true;
      tokenizer = new StringTokenizer(complexStr, ".+-Ei", true);
      String t;
      buffer = new StringBuffer("");
      int beginIndex=0, endIndex=-1;
      int cutIndex1=-1, cutEnd1=-1;
      int cutIndex2=-1, cutEnd2=-1;
      while(tokenizer.hasMoreTokens())
      {
	t = tokenizer.nextToken();
        beginIndex=endIndex+1;
        endIndex=beginIndex+(t.length()-1);
	if(t.equals("."))
	{
          buffer.append(t);
          t = tokenizer.nextToken();
          beginIndex=endIndex+1;
          endIndex=beginIndex+(t.length()-1);
          if(cutIndex1<0 && cutEnd1<0)
          {
            cutEnd1=beginIndex;
            cutIndex1=endIndex;
          }
          else
          {
            cutEnd2=beginIndex;
            cutIndex2=endIndex;
          }
        }
	buffer.append(t);
      }
      if(cutIndex1>=0)
      {
        int possibleCut1=(cutIndex1-cutEnd1)+1;
        if(cutIndex2>=0)
        {
          int possibleCut2=(cutIndex2-cutEnd2)+1;
          while(stringWidth>=columnWidth && noMoreDigs==false)
          {
            if(possibleCut1>=possibleCut2)
            {
              buffer.deleteCharAt(cutIndex1--);
              cutIndex2--;
              cutEnd2--;
              possibleCut1--;
            }
            else
            {
              buffer.deleteCharAt(cutIndex2--);
              possibleCut2--;
            }
            if(possibleCut1<1 && possibleCut2<1)
              noMoreDigs=true;
            stringWidth=SwingUtilities.
                          computeStringWidth(fontMetrics, buffer.toString());
          }
        }
        else
        {
          while(stringWidth>=columnWidth && possibleCut1>0)
          {
            buffer.deleteCharAt(cutIndex1--);
            possibleCut1--;
            stringWidth=SwingUtilities.
                          computeStringWidth(fontMetrics, buffer.toString());
          }
        }
      }
      complexStr=buffer.toString();
    }
    setText(complexStr);
    if (hasFocus)
    {
      if(isCut)
        setForeground(lightGreen);
      else
        setForeground(Color.white);
      setBackground(darkBlue);
    }
    else if (isSelected)
    {
      setBackground(lightBlue);
    }
    return this;
  }
}


