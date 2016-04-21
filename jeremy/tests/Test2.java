package jeremy.tests;

import de.jtem.mfc.field.Complex;

public class Test2 {
    
    public static void main(String args[]) {

        GroupData P = new GroupData();
        Complex h = new Complex(3,3);
        System.out.println(h);
        System.out.println(P);
        System.out.println(P.c);


        Complex cc1 = new Complex(3,1);
        Complex cc2 = new Complex(3,2);
        P.updateC(cc1);
        System.out.println("PrintC1");
        P.printC();
        P.c = cc2;
        System.out.println("PrintC2");
        P.printC();
        System.out.println("PrintC3");
        P.updateC(cc1);
        System.out.println("PrintC4");
        P.printC();

    }
}

