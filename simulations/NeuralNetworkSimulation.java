import java.util.*;
import java.io.*;

/**
 * This program simulates neural network interactions between excitatory and inhibitory neurons.
 * It initializes large arrays to represent these neurons and processes their interactions over a simulated time.
 * The neural activity, including voltage changes and spike transmissions, is logged to an output file.
 * Written as a playground/bad code example by Athanasios Theocharis (~2014), with feedback by Asterios Arampatzis.
 */

public class NeuralNetworkSimulation {

    public static void main(String[] args) throws IOException{
        PrintStream out = new PrintStream(new FileOutputStream("output.txt"));
        System.setOut(out);
        Random r = new Random();
        double[][] ex = new double[10000][130];
        double[][] in = new double[2500][130];
        int i, j, k, n;
        double t, noise, y;
        boolean same = true;
        
        //Excitatory and Inhibitory arrays
        for (i = 0; i < 10000; i++) {
            if (i < 2500) {
                in[i][0] = 10;
                in[i][129] = in[i][0];
                for (j = 1; j < 101; j++) {
                    in[i][j] = r.nextInt(10000) + 1;
                }
                for (j = 101; j < 126; j++) {
                    in[i][j] = r.nextInt(2500) + 1;
                }
            }
            ex[i][0] = 10;
            ex[i][129] = ex[i][0];
            for (j = 1; j < 101; j++) {
                ex[i][j] = r.nextInt(10000) + 1;
            }
            for (j = 101; j < 126; j++) {
                ex[i][j] = r.nextInt(2500) + 1;
            }
        }
        
        //Arrays created
        
        //Checking for double excitatory-excitatory neurons
        for (i = 0; i < 10000; i++) {
            j = 1;
            while (j < 101) {
                k = 1;
                while (k < j && same == true) {
                    if (ex[i][k] == ex[i][j]) {
                        ex[i][j] = r.nextInt(10000) + 1;
                        same = false;
                    }
                    k++;
                }
                if (same) { //Reason why the while was broken
                    j++;
                }
                same = true;
            }
        }

        //Checking for double excitatory-inhibitory neurons
        for (i = 0; i < 10000; i++) {
            j = 101;
            while (j < 126) {
                k = 101;
                while (k < j && same == true) {
                    if (ex[i][k] == ex[i][j]) {
                        ex[i][j] = r.nextInt(2500) + 1;
                        same = false;
                    }
                    k++;
                }
                if (same) { //Reason why the while was broken
                    j++;
                }
                same = true;
            }
        }

        //Checking for double inhibitory-excitatory neurons
        for (i = 0; i < 2500; i++) {
            j = 1;
            while (j < 101) {
                k = 1;
                while (k < j && same == true) {
                    if (in[i][k] == in[i][j]) {
                        in[i][j] = r.nextInt(10000) + 1;
                        same = false;
                    }
                    k++;
                }
                if (same) { //Reason why the while was broken
                    j++;
                }
                same = true;
            }
        }

        //Checking for double inhibitory-inhibitory neurons
        for (i = 0; i < 2500; i++) {
            j = 101;
            while (j < 126) {
                k = 101;
                while (k < j && same == true) {
                    if (in[i][k] == in[i][j]) {
                        in[i][j] = r.nextInt(2500) + 1;
                        same = false;
                    }
                    k++;
                }
                if (same) { //Reason why the while was broken
                    j++;
                }
                same = true;
            }
        }
        
        //BY THIS POINT THERE IS A NEURAL NETWORK WITH 10000 EX AND 2500 IN NEURONS CONNECTED WITH EACH OTHER
        //[0] CONTAINS VOLTS, [1-100] CONTAIN WITH WHICH EX NEURON THEY'RE CONNECTED AND [101-125] CONTAIN WITH WHICH IN NEURON THEY'RE CONNECTED
        for (i = 0; i < 10000; i++) { //For every neuron
            if (i < 2500) {
                for (j = 1; j < 101; j++) {
                    n = (int) in[i][j] - 1;
                    if (ex[n][126] == 0) {
                        ex[n][128]++;
                    }
                }
                for (j = 101; j < 126; j++) {
                    n = (int) in[i][j] - 1;
                    if (in[n][126] == 0) {
                        in[n][128]++;
                    }
                }
                in[i][126] = 0.5;
            }
            for (j = 1; j < 101; j++) {
                n = (int) ex[i][j] - 1;
                if (ex[n][126] == 0) {
                    ex[n][127]++;
                }
            }
            for (j = 101; j < 126; j++) {
                n = (int) ex[i][j] - 1;
                if (ex[n][126] == 0) {
                    ex[n][127]++;
                }
            }
            ex[i][126] = 0.5;
        }
        
        //Running the simulation for 0.5ms in every loop	
        for (t = 1; t < 200; t = t + 0.5) { //For every 0.5ms
            System.out.print(t + " ");
            noise = gaussian();
            for (i = 0; i < 10000; i++) { //For every neuron
                if (i < 2500) { //If i < 2500 then inhibitories included
                    if (t >= in[i][126] + 2) { //Could be wrong
                        y = y(in[i][127], in[i][128], noise);
                        in[i][0] = Volt(in[i][127] + in[i][128], y, t, t - 1.5, in[i][129]); //(int k,double y, double t, double c, a)
                        if (in[i][0] >= 20) { //Sends spike to every neuron it's connected to
                            in[i][129] = 10;
                            for (j = 1; j < 101; j++) {
                                n = (int) in[i][j] - 1;
                                if ((t + 1.5) > (ex[n][126] + 2)) {
                                    ex[n][128]++;
                                }
                            }
                            for (j = 101; j < 126; j++) {
                                n = (int) in[i][j] - 1;
                                if ((t + 1.5) > (in[n][126] + 2)) {
                                    in[n][128]++;
                                }
                            }
                            in[i][126] = t;
                        } else {
                            in[i][129] = in[i][0];
                        }
                    }
                }
                if (t > ex[i][126] + 2) {
                    y = y(ex[i][127], ex[i][128], noise);
                    ex[i][0] = Volt(ex[i][127] + ex[i][128], y, t, t - 1.5, ex[i][129]);
                    if (ex[i][0] > 20) { //Every excitatory sends spike to every neuron it's connected to
                        ex[i][129] = 10;
                        if (i < 50) {
                            System.out.print(i + " ");
                        }
                        for (j = 1; j < 101; j++) {
                            n = (int) ex[i][j] - 1;
                            if ((t + 1.5) > (ex[n][126] + 2)) {
                                ex[n][127]++;
                            }
                        }
                        for (j = 101; j < 126; j++) {
                            n = (int) ex[i][j] - 1;
                            if ((t + 1.5) > (ex[n][126] + 2)) {
                                in[n][127]++;
                            }
                        }
                        ex[i][126] = t;
                    } else {
                        ex[i][129] = ex[i][0];
                    }
                }
            }
            System.out.println();
        }
    }

    public static double gaussian() {
        Random r = new Random();
        return Math.round(Math.abs(r.nextGaussian()));
    }

    public static double Volt(double k, double y, double t, double c, double a) {
        double V;
        double D = 1.5;
        int tau = 20;
        V = (2 * k * y * theta(t - 1) * theta(c + D - t * theta(1 - t) - theta(t - 1)) * theta(-c - D + theta(1 - t) + t * theta(t - 1)) * Math.exp((c + D - t) / tau)) / tau - (k * y * theta(c + D - t * theta(1 - t) - theta(t - 1)) * theta(-c - D + theta(1 - t) + t * theta(t - 1)) * Math.exp((c + D - t) / tau)) / tau + a * Math.exp(-t / tau);
        return V;
    }

    public static int theta(double x) {
        if (x >= 0) {
            return 1;
        } else {
            return 0;
        }
    }

    public static double y(double x, double z, double noise) {
        int tau = 20;
        return tau * (x * 0.1 - z * 0.3 + noise * 0.1);
    }
}