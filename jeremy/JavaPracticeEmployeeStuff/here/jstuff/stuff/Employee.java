package jstuff.stuff;

import java.io.*;

public class Employee{

    String name;
    int age;

    public Employee(String name) {
        this.name = name;
    }

    public void empAge(int empAge) {
        age = empAge;
    }

    public void printEmployee() {
        System.out.println("Name:"+name );
        System.out.println("Age:" + age);
    }
}
