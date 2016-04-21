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
import java.awt.Insets;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
import java.beans.PropertyChangeEvent;
import java.beans.PropertyChangeListener;

import javax.swing.JLabel;
import javax.swing.JPanel;
import javax.swing.JScrollPane;
import javax.swing.JTable;
import javax.swing.JTextField;
import javax.swing.event.TableModelEvent;
import javax.swing.event.TableModelListener;
import javax.swing.table.DefaultTableColumnModel;
import javax.swing.table.TableColumn;

import de.jtem.mfc.field.Complex;

/**
 * An instance of this class is a JPanel containing a JTable representing a
 * Matrix. A MatrixEditorPanel instance is used for every kind of MatrixEditor.
 */
public class MatrixEditorPanel extends JPanel
  implements TableModelListener, ActionListener, java.io.Serializable
{
  private JTable table, rowHeaderTable;
  private IntegerMatrix rowHeaderMatrix;
  private JScrollPane sPane;
  private JTextField rowTF, colTF;
  private AbstractMatrixTableModel tableModel;
  private LineNumberModel lineNumberModel;
  private AbstractMatrix matrix;
  private GridBagConstraints gbc=new GridBagConstraints();
  //finals
  private final int ROW_HEIGHT=16;
  private final int REST_HEIGHT=36;
  private final int MAX_VISIBLE_HEIGHT=356;
  /**
   * Constructs a new MatrixEditorPanel.
   */
  public MatrixEditorPanel()
  {
    setLayout(new GridBagLayout());
    gbc.insets.top=10;
    gbc.insets.bottom=gbc.insets.right=gbc.insets.left=4;
    gbc.gridwidth=1;
    gbc.weightx=0.f;
    gbc.anchor=GridBagConstraints.SOUTHEAST;
    add(new JLabel("rows"), gbc);
    gbc.gridwidth=2;
    add(rowTF=new JTextField(), gbc);
    rowTF.addActionListener(this);
    rowTF.setColumns(8);
    gbc.gridwidth=3;
    gbc.weightx=1.f;
    add(new JLabel("colums", JLabel.RIGHT), gbc);
    gbc.gridwidth=GridBagConstraints.REMAINDER;
    gbc.weightx=0.f;
    add(colTF=new JTextField(), gbc);
    colTF.addActionListener(this);
    colTF.setColumns(8);
    gbc.insets.top=0;

    table=new JTable()
    {
      public void validate()
      {
        tableFitsInPane();
        super.validate();
      }
    };
    table.setDefaultRenderer(Complex.class, new ComplexTableCellRenderer());
    table.setDefaultEditor(Complex.class, new ComplexTableCellEditor());
    table.setCellSelectionEnabled(true);
    sPane=new JScrollPane(table, JScrollPane.VERTICAL_SCROLLBAR_AS_NEEDED,
      JScrollPane.HORIZONTAL_SCROLLBAR_ALWAYS);
    table.setAutoResizeMode(JTable.AUTO_RESIZE_OFF);
    columnModel=(DefaultTableColumnModel)table.getColumnModel();

    lineNumberModel=new LineNumberModel();
    rowHeaderTable=new JTable(lineNumberModel);
    Color color=table.getTableHeader().getBackground();
    rowHeaderTable.setBackground(color);
    rowHeaderTable.setIntercellSpacing(new Dimension(20, 0));
    rowHeaderTable.setEnabled(false);
    rowHeaderTable.setPreferredScrollableViewportSize(new Dimension(60,20));
    sPane.setRowHeaderView(rowHeaderTable);
    gbc.fill=GridBagConstraints.HORIZONTAL;
    add(sPane, gbc);
  }
  /**
   * Sets the relation between this object and the used Matrix object.
   * The specified AbstractMatrix object defines the kind of matrix working
   * together and provides the needed reference to this class. The
   * AbstractMatrixTableModel object transfers the right TableModel storing
   * all data for the matrix.
   */
  public void setMatrix(AbstractMatrix m, AbstractMatrixTableModel model)
  {
    matrix=m;

    tableModel=model;    
    tableModel.setMatrix(matrix);
    tableModel.addTableModelListener(this);
    table.setModel(tableModel);

    int rows=tableModel.getRowCount();
    int cols=tableModel.getColumnCount();
    
    rowTF.setText(String.valueOf(rows));
    colTF.setText(String.valueOf(cols));

    if(rows*ROW_HEIGHT+REST_HEIGHT>MAX_VISIBLE_HEIGHT)	   
      sPane.setPreferredSize(new Dimension(100, MAX_VISIBLE_HEIGHT));
    else
      sPane.setPreferredSize(new Dimension(100, rows*ROW_HEIGHT+REST_HEIGHT));
    lineNumberModel.setRowCount(rows);
    noCells();
    tableFitsInPane();
  }
  private void noCells()
  {
    if(rowTF.getText().equals("0") || colTF.getText().equals("0"))
    {
      table.getTableHeader().setVisible(false);
      rowHeaderTable.setVisible(false);
      sPane.setHorizontalScrollBarPolicy(
        JScrollPane.HORIZONTAL_SCROLLBAR_NEVER);
      sPane.setPreferredSize(new Dimension(100, REST_HEIGHT));
    }
    else
    {
      table.getTableHeader().setVisible(true);
      rowHeaderTable.setVisible(true);
      sPane.setHorizontalScrollBarPolicy(
        JScrollPane.HORIZONTAL_SCROLLBAR_ALWAYS);
      if(tableModel.getRowCount()*ROW_HEIGHT+REST_HEIGHT>MAX_VISIBLE_HEIGHT)	   
        sPane.setPreferredSize(new Dimension(100, MAX_VISIBLE_HEIGHT));
      else
	sPane.setPreferredSize(
	  new Dimension(100, tableModel.getRowCount()*ROW_HEIGHT+REST_HEIGHT));
    }	    
  }
  private DefaultTableColumnModel columnModel;
  private int columnWidth;
  PropertyChangeListener columnListener=new PropertyChangeListener()
  {
    public void propertyChange(PropertyChangeEvent ev)
    {
      if(!"width".equals(ev.getPropertyName())) return;
      TableColumn c=(TableColumn)ev.getSource();
      int tabWidth=table.getPreferredSize().width;
      Insets ins=sPane.getViewport().getInsets();
      int panWidth=sPane.getViewport().getWidth()-ins.left-ins.right;
      if(tabWidth<panWidth)
      {
        int delta=panWidth-tabWidth;
        TableColumn lastColumn=
          columnModel.getColumn(table.getColumnCount()-1);
        if(c!=lastColumn)
          lastColumn.setPreferredWidth(
            lastColumn.getPreferredWidth()+delta);
        else
        {
          TableColumn firstColumn=columnModel.getColumn(0);
          firstColumn.setPreferredWidth(
            firstColumn.getPreferredWidth()+delta);
        }
      }
    }
  };
  private void tableFitsInPane()
  {
    final int columnCount=tableModel.getColumnCount();
    if(columnCount>0)
    {
      int tableWidth = table.getPreferredSize().width;
      Insets ins=sPane.getViewport().getInsets();
      int paneWidth=sPane.getViewport().getWidth()-ins.left-ins.right;
      if(tableWidth<paneWidth)
      {
        columnWidth=paneWidth/columnCount;
        int rest=paneWidth-columnWidth*columnCount;
        for(int i=0; i<columnCount-1; i++)
          columnModel.getColumn(i).setPreferredWidth(columnWidth);
        columnModel.
          getColumn(columnCount-1).setPreferredWidth(columnWidth+rest);
      }
      for(int i=0; i<columnCount; i++)
        columnModel.getColumn(i).addPropertyChangeListener(columnListener);
    }
  }
  public void validate()
  {
    super.validate();
    tableFitsInPane();
  }
  /**
   * Sets the relation between this object and the specified IntegerMatrix
   * object.
   */
  public void setMatrix(IntegerMatrix m)
  {
    setMatrix(m, new IntegerMatrixTableModel());
  }
  /**
   * Sets the relation between this object and the specified RealMatrix
   * object.
   */
  public void setMatrix(RealMatrix m)
  {
    setMatrix(m, new RealMatrixTableModel());
  }
  /**
   * Sets the relation between this object and the specified ComplexMatrix
   * object.
   */
  public void setMatrix(ComplexMatrix m)
  {
    setMatrix(m, new ComplexMatrixTableModel());
  }
  /**
   * Updates the related Matrix object whenever TableModel datas were changed.
   */
  public void tableChanged(TableModelEvent ev)
  {
    matrix=tableModel.getMatrix();
  }
  /**
   * Updates the count of rows respectively columns of the displayed table and
   * the related matrix whenever the responsible textfield was edited.
   */
  public void actionPerformed(ActionEvent ev)
  {
    JTextField source=(JTextField)ev.getSource();
    int n=Integer.parseInt(source.getText());

    if(source==rowTF)
    {
      if(n==tableModel.getRowCount())
        return; 
      matrix.setNumRows(n);
      lineNumberModel.setRowCount(n);
      tableModel.setMatrix(matrix);
    }
    else if(source==colTF)
    {
      if(n==tableModel.getColumnCount())
        return;
      matrix.setNumCols(n);
      tableModel.setMatrix(matrix);
    }
    noCells();
    tableFitsInPane();
  }
}
