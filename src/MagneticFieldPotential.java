public class MagneticFieldPotential {
    private double[][] phi;
    private int N_inside;
    private int N_outside;
    private int N_surface;
    private int N_sum;
    private double[] S;
    private double Chi = 4.2;
    private int N = 100;
    private double[] localR;
    private double[] localZ;
    private double[] Integral;
    private double[] c;
    private double[] cForK;
    private double[] cUnderK;
    private double[] cAboveK;
    private double epsilon = 0.0001;
    private double q = 0.6;
    private int indexK = 3;
    private
    private Pair[][] G;

    MagneticFieldPotential(Pair[][] G, int N_inside, int N_outside, int N_surface, int indexK){
        //this.indexK = (int)Math.floor(N_outside / 2);
        this.N_inside = N_inside;
        this.N_outside = N_outside;
        this.N_surface = N_surface;
        this.N_sum = N_inside + N_outside;
        this.G = new Pair[G.length][G.length];
        for(int i = 0; i <  G.length; i++){
            for(int j = 0; j < G.length; j++){
                this.G[i][j] = new Pair(G[i][j].getFirst(), G[i][j].getSecond());
            }
        }
        phi = new double[N_surface + 1][N_sum + 1];
        for(int j = 0; j <= N_sum; j++){
            phi[0][j] = 0.0;
        }
        for(int i = 1; i <= N_surface; i++){
            phi[i][N_sum] = G[i][N_inside].getSecond();
        }
        for(int i = 1; i <= N_surface; i++){
            for(int j = 0; j <= N_sum - 1; j++){
                phi[i][j] = 0.0;
            }
        }
        this.localR = new double[N + 1];
        this.localZ = new double[N + 1];
        for(int i = 0; i <= N; i++){
            localZ[i] = Math.cos(Math.PI / 2.0 * i / N);
            localR[i] = Math.sin(Math.PI / 2.0 * i / N);
        }
        calculateS();
        calculateIntegral();
        calculateInsideC();
        calculateCforK();
        calculateCunderK();
        calculatePotential();
    }


    private void calculateS(){
        S = new double[7];
        for(int i = 1; i <= 6; i++){
            S[i] = (localR[0] - localR[i]) * (localZ[i+1] - localZ[i]) - (localR[i+1] - localR[i]) * (localZ[0] - localZ[i]);
        }
    }

    private void calculateIntegral(){
        Integral = new double[7];
        for(int i = 1; i <= 6; i++){
            Integral[i] = S[i] * (1 + Chi) * (localR[0] + localR[i] + localR[i+1]) / 6.0;
        }
    }
//
//    private Pair point(int index){
//        switch (index){
//            case 0 :{
//                return G[i][j];
//            }
//        }
//    }

    private void calculateInsideC(){
        c = new double[7];
        c[0] = (Math.pow((localZ[2] - localZ[1]), 2) + Math.pow((localR[2] - localR[1]), 2)) / Math.pow(S[1], 2) * Integral[1]
                + (Math.pow((localZ[3] - localZ[2]), 2) + Math.pow((localR[3] - localR[2]), 2)) / Math.pow(S[2], 2) * Integral[2]
                + (Math.pow((localZ[4] - localZ[3]), 2) + Math.pow((localR[4] - localR[3]), 2)) / Math.pow(S[3], 2) * Integral[3]
                + (Math.pow((localZ[5] - localZ[4]), 2) + Math.pow((localR[5] - localR[4]), 2)) / Math.pow(S[4], 2) * Integral[4]
                + (Math.pow((localZ[6] - localZ[5]), 2) + Math.pow((localR[6] - localR[5]), 2)) / Math.pow(S[5], 2) * Integral[5]
                + (Math.pow((localZ[1] - localZ[6]), 2) + Math.pow((localR[1] - localR[6]), 2)) / Math.pow(S[6], 2) * Integral[6];
        c[1] = ((localZ[1] - localZ[6]) * (localZ[6] - localZ[0]) + (localR[1] - localR[6]) * (localR[6] - localR[0])) / Math.pow(S[6], 2) * Integral[6]
                + ((localZ[2] - localZ[1]) * (localZ[0] - localZ[2]) + (localR[2] - localR[1]) * (localR[0] - localR[2])) / Math.pow(S[1], 2) * Integral[1];
        c[2] = ((localZ[2] - localZ[1]) * (localZ[1] - localZ[0]) + (localR[2] - localR[1]) * (localR[1] - localR[0])) / Math.pow(S[1], 2) * Integral[1]
                + ((localZ[3] - localZ[2]) * (localZ[0] - localZ[3]) + (localR[3] - localR[2]) * (localR[0] - localR[3])) / Math.pow(S[2], 2) * Integral[2];
        c[3] = ((localZ[3] - localZ[2]) * (localZ[2] - localZ[0]) + (localR[3] - localR[2]) * (localR[2] - localR[0])) / Math.pow(S[2], 2) * Integral[2]
                + ((localZ[4] - localZ[3]) * (localZ[0] - localZ[4]) + (localR[4] - localR[3]) * (localR[0] - localR[4])) / Math.pow(S[3], 2) * Integral[3];
        c[4] = ((localZ[4] - localZ[3]) * (localZ[3] - localZ[0]) + (localR[4] - localR[3]) * (localR[3] - localR[0])) / Math.pow(S[3], 2) * Integral[3]
                + ((localZ[5] - localZ[4]) * (localZ[0] - localZ[5]) + (localR[5] - localR[4]) * (localR[0] - localR[5])) / Math.pow(S[4], 2) * Integral[4];
        c[5] = ((localZ[5] - localZ[4]) * (localZ[4] - localZ[0]) + (localR[5] - localR[4]) * (localR[4] - localR[0])) / Math.pow(S[4], 2) * Integral[4]
                + ((localZ[6] - localZ[5]) * (localZ[0] - localZ[6]) + (localR[6] - localR[5]) * (localR[0] - localR[6])) / Math.pow(S[5], 2) * Integral[5];
        c[6] = ((localZ[6] - localZ[5]) * (localZ[5] - localZ[0]) + (localR[6] - localR[5]) * (localR[5] - localR[0])) / Math.pow(S[5], 2) * Integral[5]
                + ((localZ[1] - localZ[6]) * (localZ[0] - localZ[1]) + (localR[1] - localR[6]) * (localR[0] - localR[1])) / Math.pow(S[6], 2) * Integral[6];
    }

    private void calculateCforK(){
        cForK = new double[7];
        cForK[0] = (Math.pow((localZ[6] - localZ[5]), 2) + Math.pow((localR[6] - localR[5]), 2)) / Math.pow(S[5], 2) * Integral[5]
                + (Math.pow((localZ[1] - localZ[6]), 2) + Math.pow((localR[1] - localR[6]), 2)) / Math.pow(S[6], 2) * Integral[6];
        cForK[1] = ((localZ[1] - localZ[6]) * (localZ[6] - localZ[0]) + (localR[1] - localR[6]) * (localR[6] - localR[0])) / Math.pow(S[6], 2) * Integral[6];
        cForK[2] = cForK[3] = cForK[4] = 0;
        cForK[5] = ((localZ[6] - localZ[5]) * (localZ[0] - localZ[6]) + (localR[6] - localR[5]) * (localR[0] - localR[6])) / Math.pow(S[5], 2) * Integral[5];
        cForK[6] = ((localZ[6] - localZ[5]) * (localZ[5] - localZ[0]) + (localR[6] - localR[5]) * (localR[5] - localR[0])) / Math.pow(S[5], 2) * Integral[5]
                + ((localZ[1] - localZ[6]) * (localZ[0] - localZ[1]) + (localR[1] - localR[6]) * (localR[0] - localR[1])) / Math.pow(S[6], 2) * Integral[6];
    }

    private void calculateCunderK(){
        cUnderK = new double[7];
        cUnderK[0] = (Math.pow((localZ[2] - localZ[1]), 2) + Math.pow((localR[2] - localR[1]), 2)) / Math.pow(S[1], 2) * Integral[1]
                + (Math.pow((localZ[6] - localZ[5]), 2) + Math.pow((localR[6] - localR[5]), 2)) / Math.pow(S[5], 2) * Integral[5]
                + (Math.pow((localZ[1] - localZ[6]), 2) + Math.pow((localR[1] - localR[6]), 2)) / Math.pow(S[6], 2) * Integral[6];
        cUnderK[1] = ((localZ[1] - localZ[6]) * (localZ[6] - localZ[0]) + (localR[1] - localR[6]) * (localR[6] - localR[0])) / Math.pow(S[6], 2) * Integral[6]
                + ((localZ[2] - localZ[1]) * (localZ[0] - localZ[2]) + (localR[2] - localR[1]) * (localR[0] - localR[2])) / Math.pow(S[1], 2) * Integral[1];
        cUnderK[2] = ((localZ[2] - localZ[1]) * (localZ[1] - localZ[0]) + (localR[2] - localR[1]) * (localR[1] - localR[0])) / Math.pow(S[1], 2) * Integral[1];
        cUnderK[3] = cUnderK[4] = 0;
        cUnderK[5] = ((localZ[6] - localZ[5]) * (localZ[0] - localZ[6]) + (localR[6] - localR[5]) * (localR[0] - localR[6])) / Math.pow(S[5], 2) * Integral[5];
        cUnderK[6] = ((localZ[6] - localZ[5]) * (localZ[5] - localZ[0]) + (localR[6] - localR[5]) * (localR[5] - localR[0])) / Math.pow(S[5], 2) * Integral[5]
                + ((localZ[1] - localZ[6]) * (localZ[0] - localZ[1]) + (localR[1] - localR[6]) * (localR[0] - localR[1])) / Math.pow(S[6], 2) * Integral[6];
    }

    public void calculatePotential(){
        for(int i = 1; i <= N_surface - 1; i++){
            double prevPhi = 0;
            for(int j = N_sum - 1; j >= 1; j--){
                while(true){
                    prevPhi = phi[i][j];
                    phi[i][j] = (
                            c[2] * phi[i+1][j] + c[3] * phi[i+1][j-1] + c[4] * phi[i][j-1]
                                    + c[5] * phi[i-1][j] + c[6] * phi[i-1][j+1]
                                    + c[1] * phi[i][j+1]
                    ) / c[0];
                    phi[i][j] = (1 - q) * prevPhi + q * phi[i][j];
                    if(Math.abs(phi[i][j] - prevPhi) <= epsilon){
                        break;
                    }
                }
            }
        }
        while(true){
            double prevPhi = phi[indexK][0];
            phi[indexK][0] = (
                    cForK[5] * phi[indexK-1][0] + cForK[6] * phi[indexK-1][1]
                            + cForK[1] * phi[indexK][1]
            ) / cForK[0];
            phi[indexK][0] = (1 - q) * prevPhi + q * phi[indexK][0];
            if(Math.abs(phi[indexK][0] - prevPhi) <= epsilon){
                break;
            }
        }
        for(int i = indexK + 1; i <= N_surface - 1; i ++){
            while(true){
                double prevPhi = phi[i][0];
                phi[i][0] = (
                        cUnderK[2] * phi[i+1][0] + cUnderK[5] * phi[i-1][0] + cUnderK[6] * phi[i-1][1]
                                + cUnderK[1] * phi[indexK][1]
                ) / cUnderK[0];
                phi[i][0] = (1 - q) * prevPhi + q * phi[i][0];
                if(Math.abs(phi[i][0] - prevPhi) <= epsilon){
                    break;
                }
            }
        }
    }

    public double[][] getPotential(){
        return phi;
    }
}

