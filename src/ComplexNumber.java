import java.lang.Math.*;

public class ComplexNumber {
    private double re = 0;
    private double im = 0;

    public ComplexNumber(double re, double im) {
        this.re = re;
        this.im = im;
//        System.out.println("Complex number " + this.toString() + " was created");
    }

    public ComplexNumber(double re) {
        this.re = re;
        this.im = 0;
//        System.out.println("Complex number " + this.toString() + " was created copying another Real Number");
    }

    public ComplexNumber(ComplexNumber that) {
        this.re = that.re;
        this.im = that.im;
//        System.out.println("Complex number " + this.toString() + " was created copying another Complex Number");
    }

    public ComplexNumber() {
//        System.out.println("Complex number " + this.toString() + " was created with default value");
    }

//    // construct a complex number from polar coordinates
//    public static ComplexNumber fromPolar(double r, double theta) {
//        return new ComplexNumber(r*Math.cos(theta), r*Math.sin(theta));
//    }

    public double reOutput() {
        return this.re;
    }

    public double imOutput() {
        return this.im;
    }

    public double absSquared() {
        return this.re*this.re + this.im*this.im;
    }

    public double abs() {
        return Math.sqrt(this.absSquared());
    }

    public ComplexNumber complexConjugate() {
        return new ComplexNumber(this.re, -this.im);
    }

    public ComplexNumber add(ComplexNumber that) {
        return new ComplexNumber(this.re+that.re, this.im+that.im);
    }

    public ComplexNumber add(double that) {
        return new ComplexNumber(this.re+that, this.im);
    }

    public ComplexNumber subtract(ComplexNumber that) {
        return new ComplexNumber(this.re-that.re, this.im-that.im);
    }

    public ComplexNumber subtract(double that) {
        return new ComplexNumber(this.re-that, this.im);
    }

    public ComplexNumber multiply(ComplexNumber that) {
        double reNew = this.re * that.re - this.im * that.im;
        double imNew = this.re * that.im + this.im * that.re;
        return new ComplexNumber(reNew, imNew);
    }

    public ComplexNumber multiply(double that) {
        return new ComplexNumber(this.re*that, this.im*that);
    }

    public String toString() {
        if (this.re == 0) {
            return "" + this.im + "i";
        }else if (this.im == 0) {
            return "" + this.re;
        }else {
            return "" + this.re + "+" + this.im + "i";
        }
    }

    public boolean equals(Object obj) {
        if (obj == null) {
            System.out.println("target obj null");
            return false;
        }else if (!(obj instanceof ComplexNumber)) {
            System.out.println("target obj not a ComplexNumber");
            return false;
        }else {
            ComplexNumber that = (ComplexNumber) obj;
            if (this.re != that.re) {
                return false;
            }else if (this.im != that.im) {
                return false;
            }
        }
        return true;
    }

    public void show() {
        System.out.println(this.toString());
    }
}
