import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.IOException;

public class tester {
    public static void main(String[] args) {
        double[][] rho_eq_data = new double[1000][1000];
        try {
            String path = System.getProperty("user.dir");
            File file = new File(path + "\\src\\rand.txt");
            FileReader fileReader = new FileReader(file);
            BufferedReader bufferedReader = new BufferedReader(fileReader);
            String line;
            int idx = 0;
            while ((line = bufferedReader.readLine()) != null) {
                String[] lineArray = line.split("  ");
                for(int idxx = 1; idxx < lineArray.length; idxx++) {
                    rho_eq_data[idx][idxx-1] = Double.parseDouble(lineArray[idxx]);
                }
                idx++;
            }
            fileReader.close();
        } catch (IOException e) {
            e.printStackTrace();
        }
        Matrix rho_eq = new Matrix(rho_eq_data);


        double tstart = 0;
        tstart = System.currentTimeMillis();
        Matrix x = rho_eq.multiply(rho_eq);
        double tend = 0;
        tend = System.currentTimeMillis();
        x.show();
        System.out.println("Calculation took : " + (tend - tstart) + " mills");

    }
}
