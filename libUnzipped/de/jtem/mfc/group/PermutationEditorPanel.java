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

package de.jtem.mfc.group;

import java.awt.GridBagConstraints;
import java.awt.GridBagLayout;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;

import javax.swing.JButton;
import javax.swing.JLabel;
import javax.swing.JPanel;
import javax.swing.JTextField;


/**
 * An instance of this class is a JPanel containing two JTextField's, one for
 * the permutation, the other for its cycles, and two buttons, next and
 * nextDerangement, the latter can be greyed. A PermutationEditorPanel
 * instance is used for every kind of PermutationEditor.
 */
public class PermutationEditorPanel extends JPanel
  implements  ActionListener, java.io.Serializable
{
    private JTextField permutation;
    private JTextField label;
    private JTextField cycle;
    private JButton next;
    private JButton nextD;
    private GridBagConstraints gbc=new GridBagConstraints();

    private int[] value;

    private PermutationEditor pE;
    /**
     * Constructs a new PermutationEditorPanel.
     */
    public PermutationEditorPanel(PermutationEditor pE)
    {
	this.pE = pE;

	setLayout(new GridBagLayout());

	gbc.fill=GridBagConstraints.NONE;
	gbc.insets.top=gbc.insets.bottom=0;
	gbc.insets.right=gbc.insets.left=4;
	gbc.gridwidth=1;
	gbc.weightx=0.f;
	add(new JLabel("", JLabel.RIGHT), gbc);
	
	gbc.fill=GridBagConstraints.HORIZONTAL;
	gbc.gridwidth=GridBagConstraints.REMAINDER;
	gbc.weightx=1.f;
	add(label=new JTextField("( 0 1 2 3 4 )",JLabel.LEFT), gbc);
	label.setEditable(false);

	gbc.insets.bottom=4;

	gbc.gridwidth=1;
	gbc.weightx=0.f;
	add(new JLabel("permutation", JLabel.RIGHT), gbc);


	gbc.gridwidth=GridBagConstraints.REMAINDER;
	gbc.weightx=1.f;
	add(permutation=new JTextField(), gbc);
	permutation.addActionListener(this);


	gbc.insets.top=10;



	gbc.gridwidth=1;
	gbc.weightx=0.f;
	add(new JLabel("cycles", JLabel.RIGHT), gbc);

	gbc.gridwidth=GridBagConstraints.REMAINDER;
	gbc.weightx=1.f;
	add(cycle=new JTextField(), gbc);
	cycle.addActionListener(this);


	gbc.gridwidth=1;
	gbc.weightx=0.;
	gbc.fill=GridBagConstraints.NONE;

	add(next = new JButton("next"), gbc);
	next.setActionCommand("next");
	next.setToolTipText("Click to jump to the next permutation.");
	next.addActionListener(this);

	gbc.gridwidth=GridBagConstraints.REMAINDER;
	gbc.weightx=1.f;
	add(nextD = new JButton("next Derangement"), gbc);
	nextD.setActionCommand("nextD");
	nextD.setToolTipText("Click to jump to the next derangement.");
	nextD.setEnabled(false);
	nextD.addActionListener(this);
	gbc.insets.top=0;

    }
    /**
     * Sets the relation between this object and the used Permutation object.
     */
    public void setValue(Permutation v)
    {
	setValue(v.getPermutation());
    }

    public void setValue(int[] v)
    {

	value = v;
	pE.setValue(v);
	    
	label.setText(de.jtem.numericalMethods.algebra.group.Permutation.toString
		      (de.jtem.numericalMethods.algebra.group.Permutation.identity
		       (v.length)));

	permutation.setText
	    (de.jtem.numericalMethods.algebra.group.Permutation.toString(v));
	cycle.setText
	    (de.jtem.numericalMethods.algebra.group.Permutation.cyclesToString
	     (de.jtem.numericalMethods.algebra.group.Permutation.cycles(v)));
	nextD.setEnabled
	    (de.jtem.numericalMethods.algebra.group.Permutation.isDerangement(v));
    }
    /**
     * Updates the value whenever a button was pushed or something was edited.
     */
    public void actionPerformed(ActionEvent ev)
    {
	if ("next".equals(ev.getActionCommand())) {
	    de.jtem.numericalMethods.algebra.group.Permutation.next(value);
	    setValue(value);
	    return;
	}
	if ("nextD".equals(ev.getActionCommand())) {

	    de.jtem.numericalMethods.algebra.group.Permutation.nextDerangement
		     (value);
	    setValue(value);
	    return;
	}
	// Otherwise it comes from the JTextField's.
	JTextField source=(JTextField)ev.getSource();
	String s = source.getText();

	if (source==permutation)  
	    setValue
		(de.jtem.numericalMethods.algebra.group.Permutation.stringToIntArray(s));  

	if (source==cycle)  
	    setValue
		(de.jtem.numericalMethods.algebra.group.Permutation.fromCycles
		 (de.jtem.numericalMethods.algebra.group.Permutation.stringToCycles(s)));  
	
    }
}

