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

import java.awt.GridBagConstraints;
import java.awt.GridBagLayout;
import java.beans.Customizer;

import javax.swing.JLabel;
import javax.swing.JPanel;
import javax.swing.JScrollPane;
import javax.swing.JTable;
import javax.swing.JTextField;
import javax.swing.event.TableModelEvent;
import javax.swing.event.TableModelListener;
import javax.swing.table.DefaultTableColumnModel;
import javax.swing.table.TableColumn;

abstract class AbstractVectorCustomizer extends JPanel implements Customizer, TableModelListener 
{
    protected AbstractVectorTableModel tableModel;
    protected JTable table;
    protected JScrollPane sPane;
    protected JTextField  numRows;
    protected boolean loop;

    public void tableChanged(TableModelEvent e) {
	if(!loop) try {
	    loop=true;
	    firePropertyChange( null, null, null );
	}finally{
	    loop=false;
	}
    }

    public AbstractVectorCustomizer( AbstractVectorTableModel aTableModel ) {
	super(new GridBagLayout());
	GridBagConstraints gbc=new GridBagConstraints();
	gbc.fill=gbc.HORIZONTAL;
	tableModel = aTableModel; 
	aTableModel.addTableModelListener( this );
	table=new JTable(tableModel);
 // table.setOpaque(true);
	table.setDoubleBuffered(false);
	sPane=new JScrollPane(table,JScrollPane.VERTICAL_SCROLLBAR_AS_NEEDED,
			      JScrollPane.HORIZONTAL_SCROLLBAR_NEVER);
	add(new JLabel("Number of Elements: "),gbc);
	gbc.weightx=1.0f;
	gbc.insets.left=10;
	numRows=new JTextField();
	numRows.setEditable(false);
	add(numRows,gbc);
	gbc.weightx=1.5f;
	gbc.gridwidth=gbc.REMAINDER;
	add(new JLabel(""),gbc);
	gbc.insets.left=0;
	gbc.fill=gbc.BOTH;
	gbc.weighty=1.0f;
	add(sPane,gbc);
	sPane.setDoubleBuffered(false);
	DefaultTableColumnModel mod=(DefaultTableColumnModel)table.getColumnModel();
	TableColumn tc=mod.getColumn(0);
	/*    tc.setMinWidth(80);
	      tc.setPreferredWidth(130);
	      tc.setMaxWidth(Integer.MAX_VALUE);
	      tc.setResizable(false); */
    }

 public void setObject(Object o){
     try{
	 loop=true;
	 AbstractVector m=(AbstractVector)o;
	 tableModel.setVector(m);
	 numRows.setText(String.valueOf(m.size()));
     } finally{
	 loop=false;
     }
  }

}
















