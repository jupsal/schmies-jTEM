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
import java.awt.Dimension;
import java.awt.GridBagConstraints;
import java.awt.GridBagLayout;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;

import javax.swing.JLabel;
import javax.swing.JPanel;
import javax.swing.JScrollPane;
import javax.swing.JTable;
import javax.swing.JTextField;
import javax.swing.event.TableModelEvent;
import javax.swing.event.TableModelListener;

import de.jtem.mfc.field.Complex;

/**
 * An instance of this class is a JPanel containing a one columned JTable
 * representing a Vector. A VectorEditorPanel instance is used for every
 * kind of VectorEditor.
 */
public class VectorEditorPanel extends JPanel
  implements TableModelListener, ActionListener, java.io.Serializable
{
  private JTable table, rowHeaderTable;
  private JScrollPane sPane;
  private JTextField elementsTF;
  private AbstractVectorTableModel tableModel;
  private LineNumberModel lineNumberModel;
  private AbstractVector vector;
  private GridBagConstraints gbc=new GridBagConstraints();

  /**
   * Constructs a new VectorEditorPanel.
   */
  public VectorEditorPanel()
  {
    setLayout(new GridBagLayout());
    gbc.fill=GridBagConstraints.HORIZONTAL;
    gbc.insets.top=10;
    gbc.insets.bottom=gbc.insets.right=gbc.insets.left=4;
    gbc.gridwidth=1;
    gbc.weightx=1.f;
    gbc.anchor=GridBagConstraints.SOUTHEAST;
    add(new JLabel("number of elements", JLabel.RIGHT), gbc);
    gbc.gridwidth=GridBagConstraints.REMAINDER;
    gbc.weightx=0.f;
    add(elementsTF=new JTextField(), gbc);
    elementsTF.addActionListener(this);
    elementsTF.setColumns(8);
    gbc.insets.top=0;

    table=new JTable();
    table.setDefaultRenderer(Complex.class, new ComplexTableCellRenderer());
    table.setDefaultEditor(Complex.class, new ComplexTableCellEditor());
    //table.setAutoResizeMode(JTable.AUTO_RESIZE_OFF);
    table.setCellSelectionEnabled(true);
    sPane=new JScrollPane(table,
		          JScrollPane.VERTICAL_SCROLLBAR_AS_NEEDED,
                          JScrollPane.HORIZONTAL_SCROLLBAR_AS_NEEDED);
    lineNumberModel=new LineNumberModel();
    rowHeaderTable=new JTable(lineNumberModel);
    Color color=table.getTableHeader().getBackground();
    rowHeaderTable.setBackground(color);
    rowHeaderTable.setIntercellSpacing(new Dimension(20, 0));
    //rowHeaderTable-cells should not be selectable
    rowHeaderTable.setEnabled(false);
    rowHeaderTable.setPreferredScrollableViewportSize(new Dimension(60,20));
    sPane.setRowHeaderView(rowHeaderTable);
    //sPane.setCorner(JScrollPane.UPPER_LEFT_CORNER, new JButton("element"));
    add(sPane, gbc);
  }
  /**
   * Sets the relation between this object and the used Vector object.
   * The specified AbstractVector object defines the kind of vector working
   * together and provides the needed reference to this class. The
   * AbstractVectorTableModel object transfers the right TableModel storing
   * all data for the vector.
   */
  public void setVector(AbstractVector v, AbstractVectorTableModel model)
  {
    vector=v;

    tableModel=model;    
    tableModel.setVector(vector);
    tableModel.addTableModelListener(this);
    table.setModel(tableModel);

    int elements=tableModel.getRowCount();
       
    elementsTF.setText(String.valueOf(elements));
   
    //maximal 20 Elemente sichtbar
    if(elements*16+21>341)	   
      sPane.setPreferredSize(new Dimension(100, 341));
    else
      sPane.setPreferredSize(new Dimension(100, elements*16+21));
    lineNumberModel.setRowCount(elements);
  }
  /**
   * Sets the relation between this object and the specified IntegerVector
   * object.
   */
  public void setVector(IntegerVector v)
  {
    setVector(v, new IntegerVectorTableModel());
  }
  /**
   * Sets the relation between this object and the specified RealVector
   * object.
   */
  public void setVector(RealVector v)
  {
    setVector(v, new RealVectorTableModel());
  }
  /**
   * Sets the relation between this object and the specified ComplexVector
   * object.
   */
  public void setVector(ComplexVector v)
  {
    setVector(v, new ComplexVectorTableModel());
  }
  /**
   * Updates the related Vector object whenever TableModel datas were changed.
   */
  public void tableChanged(TableModelEvent ev)
  {
    vector=tableModel.getVector();
  }
  /**
   * Updates the count of rows of the displayed table and the related vector
   * whenever the responsible textfield was edited.
   */
  public void actionPerformed(ActionEvent ev)
  {
    JTextField source=(JTextField)ev.getSource();
    int n=Integer.parseInt(source.getText());

    if(n==tableModel.getRowCount())
      return; 
    if(n*16+21>341)	   
      sPane.setPreferredSize(new Dimension(100, 341));
    else
      sPane.setPreferredSize(new Dimension(100, n*16+21));
    vector.setNumEntries(n);
    lineNumberModel.setRowCount(n);
    tableModel.setVector(vector);
  }
}

