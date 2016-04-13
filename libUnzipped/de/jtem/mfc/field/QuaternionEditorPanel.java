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

package de.jtem.mfc.field;

import java.awt.GridBagConstraints;
import java.awt.GridBagLayout;

import javax.swing.JPanel;
import javax.swing.event.ChangeEvent;
import javax.swing.event.ChangeListener;

/**
 * An instance of this class is used by the QuaternionEditor class to display 4
 * MinMaxPanels representing a Quaternion.
 */
public class QuaternionEditorPanel extends JPanel
{
  static final public String QUATERNION_PROPERTY="quaternion";
  MinMaxPanel rScrol, iScrol, jScrol, kScrol;
  Quaternion value;
  boolean changing;
  /**
   * Constructs a new QuaternionEditorPanel initializing 4 MinMaxPanels titled
   * r, i, j and k. The MinMaxPanel r represents the real part of a Quaternion,
   * i, j and k represent the 3 imaginary parts of a Quaternion.
   */
  public QuaternionEditorPanel()
  {
    value=new Quaternion();
    setLayout(new GridBagLayout());
    GridBagConstraints gbc=new GridBagConstraints();
    gbc.gridwidth=GridBagConstraints.REMAINDER;
    gbc.fill=GridBagConstraints.HORIZONTAL;
    gbc.weightx=1f;
    gbc.insets.top=10;
    add(rScrol=new MinMaxPanel("r"), gbc);
    rScrol.getModel().addChangeListener(new QuaternionChangeListener(0));
    gbc.insets.top=0;
    add(iScrol=new MinMaxPanel("i"), gbc);
    iScrol.getModel().addChangeListener(new QuaternionChangeListener(1));
    add(jScrol=new MinMaxPanel("j"), gbc);
    jScrol.getModel().addChangeListener(new QuaternionChangeListener(2));
    gbc.insets.bottom=10;
    add(kScrol=new MinMaxPanel("k"), gbc);
    kScrol.getModel().addChangeListener(new QuaternionChangeListener(3));
  }
  /**
   * Assigns quat to the Quaternion of this QuaternionEditorPanel.
   * Fires a PropertyChangeEvent with property name {@link
   * #QUATERNION_PROPERTY} and null for oldValue and newValue.
   * Get the updated value by calling one of the getQuaternion methods.
   */
  public void setQuaternion(Quaternion quat)
  {
    if(!value.equals(quat))
    try
    {
      changing=true;
      value.assign(quat);
      rScrol.setValue(value.getRe());
      iScrol.setValue(value.getX());
      jScrol.setValue(value.getY());
      kScrol.setValue(value.getZ());

      firePropertyChange(QUATERNION_PROPERTY, null, null);
    }finally
    {
      changing=false;
    }
  }
  public Quaternion getQuaternion()
  {
    return new Quaternion(value);
  }
  public void getQuaternion(Quaternion target)
  {
    target.assign(value);
  }
  private class QuaternionChangeListener  implements ChangeListener
  {
    final int index;
    final double[] quaternion=new double[4];
    QuaternionChangeListener(int i)
    {
      index=i;
    }
    public void stateChanged(ChangeEvent ev)
    {
      if(changing) return;
      final DoubleBoundedRangeModel src=(DoubleBoundedRangeModel)ev.getSource();
      value.get(quaternion, 0);
      quaternion[index]=src.getDoubleValue();
      value.assign(quaternion, 0);
      firePropertyChange(QUATERNION_PROPERTY, null, null);
    }
  }
}
