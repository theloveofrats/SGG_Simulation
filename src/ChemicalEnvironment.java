import javax.swing.*;
import java.awt.image.BufferedImage;
import java.awt.image.Raster;
import java.util.ArrayList;

/**
 * Created with IntelliJ IDEA.
 * User: luke
 * Date: 19/02/2014
 * Time: 08:23
 * To change this template use File | Settings | File Templates.
 */
public class ChemicalEnvironment {

    public static double grain = 20.0;
    public double[][] profile;
    public double[][] oldprofile;

    public static double production  = 0.0;
    public static double degradation = 0.00; // Equilibrium values of 20.


    public ChemicalEnvironment(double dMin, double dMax){

        int iMax = ((int) (MelaMigration.dimensions[2]/grain))-1;
        profile = new double[iMax][];

        for(int i = 0; i<iMax; i++){

            double[] dNew = new double[((int) (MelaMigration.dimensions[3]/grain))];
            for(int j = 0; j<dNew.length; j++){
                if(i<=(MigrationSimulation.padding/grain)) dNew[j] = 0.0;
                else dNew[j] = dMin+ (((double) i)/iMax)*(dMax-dMin);
            }
            profile[i] = dNew;
        }

        for(int i = 0; i<iMax; i++){

        double[] dNew = new double[((int) (MelaMigration.dimensions[3]/grain))];
        for(int j = 0; j<dNew.length; j++){
            if(i<=(MigrationSimulation.padding/grain)) dNew[j] = 0.0;
            //else if(i==iMax-1) dNew[j] = dMax;
            else dNew[j] = dMin+ (((double) i)/iMax)*(dMax-dMin);
        }
        profile[i] = dNew;
    }

        oldprofile = profile.clone();
    }

    private int[] ConvertCoordinates(double cx, double cy){

        int x = (int) Math.floor(cx/grain);
        int y = (int) Math.floor(cy/grain);

        x = Math.max(0,Math.min(x,profile.length-1));
        y = Math.max(0,Math.min(y,profile[x].length-1));

        return new int[]{x,y};
    }

    public void TakeAtLocation(double dx, double dy, double amount){
        int[] ic = ConvertCoordinates(dx,dy);

        profile[ic[0]][ic[1]] -=amount;
    }

    public void add(double cn){
        for(int i=0; i<profile.length; i++){
            for(int j=0; j<profile[i].length; j++){
                profile[i][j]+=cn;
            }
        }
    }

    public double GetLocationConcentration(double dx, double dy){
        int[] ic = ConvertCoordinates(dx,dy);
        return profile[ic[0]][ic[1]];
    }

    public double[] GetGradientAtLocation(double dx, double dy){
        int[] ic = ConvertCoordinates(dx,dy);
        double xDiff = profile[Math.min(ic[0]+1, profile.length-1)][ic[1]] - profile[(Math.max(ic[0]-1, 0))][ic[1]];
        double yDiff = profile[ic[0]][Math.min(ic[1]+1, profile[ic[0]].length-1)] - profile[ic[0]][Math.max(ic[1]-1, 0)];

        return new double[]{xDiff/(2*grain),yDiff/(2*grain)};
    }


    public double GetConcentrationAtLocation(double dx, double dy){
        int[]  ic = ConvertCoordinates(dx,dy);

        double c0 = profile[ic[0]][ic[1]];

        /*c0       += profile.get(Math.min(ic[0]+1, profile.size()-1))[ic[1]];
        c0       += profile.get(Math.max(ic[0]-1, 0))[ic[1]];
        c0       += profile.get(ic[0])[Math.min(ic[1]+1, profile.get(ic[0]).length-1)]
        c0       += profile.get(ic[0])[Math.max(ic[1]-1, 0)];        */

        return c0/5;
    }

    public void diffuse(double D, double dt,double dx){
        int iB = MelaMigration.pinned ? 1 : 0;
        /*if(MelaMigration.CNM){
            //Crank-Nicolson,
            for(int i = 0; i<profile.length; i++){
                for(int j = 0; j<profile[i].length; j++){
                    int yLength = oldprofile[i].length-1;
                    double ip1,im1,jp1,jm1;
                    ip1 = im1 = oldprofile[i][j];
                    jp1 =       oldprofile[i][0];
                    jm1 = oldprofile[i][yLength];


                    if(i > 0) im1 = oldprofile[i-1][j];
                    if(i < oldprofile.length-1) ip1 = oldprofile[i+1][j];
                    if(j > 0) jm1 = oldprofile[i][j-1];
                    if(j < yLength) jp1 = oldprofile[i][j+1];

                    double r = (D*dt/(dx*dx));


                    profile[i][j] += (D*dt/dx)*(ip1+im1+jp1+jm1-4*oldprofile[i][j]);
                }
            }
        }*/
        //else{
            for(int i = iB; i<profile.length-iB; i++){
                for(int j = 0; j<profile[i].length; j++){
                    int yLength = oldprofile[i].length-1;
                    double ip1,im1,jp1,jm1;
                    ip1 = im1 = oldprofile[i][j];
                    jp1 =       oldprofile[i][0];
                    jm1 = oldprofile[i][yLength];


                    if(i > 0) im1 = oldprofile[i-1]   [j];
                    if(i < oldprofile.length-1) ip1 = oldprofile[i+1][j];
                    if(j > 0) jm1 = oldprofile[i][j-1];
                    if(j < yLength) jp1 = oldprofile[i][j+1];

                    profile[i][j] += (D*dt/(dx*dx))*(ip1+im1+jp1+jm1-4*oldprofile[i][j])+(ChemicalEnvironment.production-ChemicalEnvironment.degradation*profile[i][j])*dt;
                }
            }
        //}
        oldprofile = profile.clone();
    }

    public double[] meanProfile(){
        double[] dOut = new double[profile.length];
        for(int i = 0; i<profile.length; i++){
            double dMean = 0;
            for(int j = 0; j<profile[i].length; j++){
                dMean+=profile[i][j];
            }
            dMean/=profile[i].length;
            dOut[i] = dMean;
        }
        return dOut;
    }
}
