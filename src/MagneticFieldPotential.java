public class MagneticFieldPotential {
    private double[][] phi;
    private double[][] prevPhi;
    private int N_surface;
    private int N_sum;
    private double[] S;
    private double Chi = 0.1;
    private double[] localR;
    private double[] localZ;
    private double[] Integral;
    private double[] c;
    private double[] cForK;
    private double[] cUnderK;
    private double[] cAboveK;
    private double epsilon = 0.0001;
    private double q = 0.6;
    private Pair[][] G;

    MagneticFieldPotential(Pair[][] G, int N_inside, int N_outside, int N_surface){
        this.N_surface = N_surface;
        this.N_sum = N_inside + N_outside;
        this.G = new Pair[N_surface + 1][N_sum + 1];
        for(int i = 0; i <= N_surface  ; i++){
            for(int j = 0; j <= N_sum; j++){
                this.G[i][j] = new Pair(G[i][j].getFirst(), G[i][j].getSecond());
            }
        }
        phi = new double[N_surface + 1][N_sum + 1];
        prevPhi = new double[N_surface + 1][N_sum + 1];
        for(int j = 0; j <= N_sum; j++){
            phi[0][j] = 0.0;
        }
        for(int i = 1; i <= N_surface; i++){
            phi[i][N_sum] = G[i][N_sum].getSecond();
        }
        for(int i = 1; i <= N_surface; i++){
            for(int j = 0; j <= N_sum - 1; j++){
                phi[i][j] = 0.0;
            }
        }
        for(int i = 0; i <= N_surface; i++){
            for(int j = 0; j <= N_sum; j++){
                prevPhi[i][j] = phi[i][j];
            }
        }
        this.localR = new double[7];
        this.localZ = new double[7];
        for(int i = 0; i < 7; i++){
            localZ[i] = 0;
            localR[i] = 0;
        }
        cForK = new double[7];
        S = new double[7];
        Integral = new double[7];
        c = new double[7];
        cUnderK = new double[7];
        cAboveK = new double[7];
    }

    public void calculatePotential(){
        calculatePotentialInsideAndAboveK();
        //calculatePotentialUnderK();
        //calculatePotentialInK();
    }

    private void calculateS(){
        for(int i = 1; i <= 5; i++){
            S[i] = Math.abs((localR[0] - localR[i]) * (localZ[i+1] - localZ[i]) - (localR[i+1] - localR[i]) * (localZ[0] - localZ[i]));
        }
        S[6] = Math.abs((localR[0] - localR[6]) * (localZ[1] - localZ[6]) - (localR[1] - localR[6]) * (localZ[0] - localZ[6]));
    }

    private void calculateIntegral(){
        for(int i = 1; i <= 5; i++){
            Integral[i] = S[i] * (1.0 + Chi) * (localR[0] + localR[i] + localR[i+1]) / 6.0;
        }
        Integral[6] = S[6] * (1.0 + Chi) * (localR[0] + localR[6] + localR[1]) / 6.0;
    }

    private void calculateInsideC(){
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
        for(int i = 0; i < 7; i++){
            c[i] = Math.abs(c[i]);
        }
    }

    private void calculateCAboveK(){
        cAboveK[0] = (Math.pow((localZ[5] - localZ[4]), 2) + Math.pow((localR[5] - localR[4]), 2)) / Math.pow(S[4], 2) * Integral[4]
                + (Math.pow((localZ[6] - localZ[5]), 2) + Math.pow((localR[6] - localR[5]), 2)) / Math.pow(S[5], 2) * Integral[5]
                + (Math.pow((localZ[1] - localZ[6]), 2) + Math.pow((localR[1] - localR[6]), 2)) / Math.pow(S[6], 2) * Integral[6];
        cAboveK[1] = ((localZ[1] - localZ[6]) * (localZ[6] - localZ[0]) + (localR[1] - localR[6]) * (localR[6] - localR[0])) / Math.pow(S[6], 2) * Integral[6];
        cAboveK[2] = cAboveK[3] = 0;
        cAboveK[4] = ((localZ[5] - localZ[4]) * (localZ[0] - localZ[5]) + (localR[5] - localR[4]) * (localR[0] - localR[5])) / Math.pow(S[4], 2) * Integral[4];
        cAboveK[5] = ((localZ[5] - localZ[4]) * (localZ[4] - localZ[0]) + (localR[5] - localR[4]) * (localR[4] - localR[0])) / Math.pow(S[4], 2) * Integral[4]
                + ((localZ[6] - localZ[5]) * (localZ[0] - localZ[6]) + (localR[6] - localR[5]) * (localR[0] - localR[6])) / Math.pow(S[5], 2) * Integral[5];
        cAboveK[6] = ((localZ[6] - localZ[5]) * (localZ[5] - localZ[0]) + (localR[6] - localR[5]) * (localR[5] - localR[0])) / Math.pow(S[5], 2) * Integral[5]
                + ((localZ[1] - localZ[6]) * (localZ[0] - localZ[1]) + (localR[1] - localR[6]) * (localR[0] - localR[1])) / Math.pow(S[6], 2) * Integral[6];
        for(int i = 0; i < 7; i++){
            cAboveK[i] = Math.abs(cAboveK[i]);
        }
    }


    private void calculateCUnderK(){
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
        for(int i = 0; i < 7; i++){
            cUnderK[i] = Math.abs(cUnderK[i]);
        }
    }

    private void calculateCforK(){
        cForK[0] = (Math.pow((localZ[6] - localZ[5]), 2) + Math.pow((localR[6] - localR[5]), 2)) / Math.pow(S[5], 2) * Integral[5]
                + (Math.pow((localZ[1] - localZ[6]), 2) + Math.pow((localR[1] - localR[6]), 2)) / Math.pow(S[6], 2) * Integral[6];
        cForK[1] = ((localZ[1] - localZ[6]) * (localZ[6] - localZ[0]) + (localR[1] - localR[6]) * (localR[6] - localR[0])) / Math.pow(S[6], 2) * Integral[6];
        cForK[2] = cForK[3] = cForK[4] = 0;
        cForK[5] = ((localZ[6] - localZ[5]) * (localZ[0] - localZ[6]) + (localR[6] - localR[5]) * (localR[0] - localR[6])) / Math.pow(S[5], 2) * Integral[5];
        cForK[6] = ((localZ[6] - localZ[5]) * (localZ[5] - localZ[0]) + (localR[6] - localR[5]) * (localR[5] - localR[0])) / Math.pow(S[5], 2) * Integral[5]
                + ((localZ[1] - localZ[6]) * (localZ[0] - localZ[1]) + (localR[1] - localR[6]) * (localR[0] - localR[1])) / Math.pow(S[6], 2) * Integral[6];
        for(int i = 0; i < 7; i++){
            cForK[i] = Math.abs(cForK[i]);
        }
    }


    public double[][] getPotential(){
        return phi;
    }
//
//    private void calculatePotentialInK(){
//        int indexK = N_surface;
//        newLocalCoord(indexK, 0);
//        calculateS();
//        calculateIntegral();
//        calculateCforK();
//        double prevPhi;
//        while(true){
//            prevPhi = phi[indexK][0];
//            phi[indexK][0] = (
//                    cForK[5] * phi[indexK-1][0] + cForK[6] * phi[indexK-1][1]
//                            + cForK[1] * phi[indexK][1]
//            ) / cForK[0];
//            phi[indexK][0] = (1 - q) * prevPhi + q * phi[indexK][0];
//            if(Math.abs(phi[indexK][0] - prevPhi) <= epsilon){
//                break;
//            }
//        }
//    }


    private void newLocalCoord(int i, int j){
        localR[0] = G[i][j].getFirst();
        localZ[0] = G[i][j].getSecond();
        if(j + 1 <= N_sum) {
            localR[1] = G[i][j + 1].getFirst();
            localZ[1] = G[i][j + 1].getSecond();
        }
        if(i + 1 <= N_surface) {
            localR[2] = G[i + 1][j].getFirst();
            localZ[2] = G[i + 1][j].getSecond();
        }
        if(j - 1 >= 0 && i + 1 <= N_surface) {
            localR[3] = G[i + 1][j - 1].getFirst();
            localZ[3] = G[i + 1][j - 1].getSecond();
        }
        if(j - 1 >= 0){
            localR[4] = G[i][j - 1].getFirst();
            localZ[4] = G[i][j - 1].getSecond();
        }
        if(i - 1 >= 0) {
            localR[5] = G[i - 1][j].getFirst();
            localZ[5] = G[i - 1][j].getSecond();
        }
        if(i - 1 >= 0 && j+1 <= N_sum){
            localR[6] = G[i - 1][j + 1].getFirst();
            localZ[6] = G[i - 1][j + 1].getSecond();
        }
    }

//    private void calculatePotentialUnderK(){
//        for(int i = N_surface - 1; i >= 1; i--){
//            newLocalCoord(i, 0);
//            calculateS();
//            calculateIntegral();
//            calculateCUnderK();
//            while(true){
//                double prevPhi = phi[i][0];
//                phi[i][0] = (
//                        cUnderK[2] * phi[i+1][0] + cUnderK[5] * phi[i-1][0] + cUnderK[6] * phi[i-1][1]
//                                + cUnderK[1] * phi[i][1]
//                ) / cUnderK[0];
//                phi[i][0] = (1 - q) * prevPhi + q * phi[i][0];
//                if(Math.abs(phi[i][0] - prevPhi) <= epsilon){
//                    break;
//                }
//            }
//        }
//    }

    private void calculatePotentialInsideAndAboveK(){
        boolean isCorrectResult = false;
        int count = 0;
        while(!isCorrectResult) {
            for (int j = N_sum - 1; j >= 1; j--) {
                for (int i = 1; i <= N_surface - 1; i++) {
                    newLocalCoord(i, j);
                    calculateS();
                    calculateIntegral();
                    calculateInsideC();
                    prevPhi[i][j] = phi[i][j];
                    phi[i][j] = (
                            c[2] * phi[i + 1][j] + c[3] * phi[i + 1][j - 1] + c[4] * phi[i][j - 1]
                                    + c[5] * phi[i - 1][j] + c[6] * phi[i - 1][j + 1]
                                    + c[1] * phi[i][j + 1]
                    ) / c[0];
                    //phi[i][j] = (1 - q) * prevPhi[i][j] + q * phi[i][j];
                }
                newLocalCoord(N_surface, j);
                calculateS();
                calculateIntegral();
                calculateCAboveK();
                prevPhi[N_surface][j] = phi[N_surface][j];
                phi[N_surface][j] = (
                        cAboveK[4] * phi[N_surface][j - 1]
                                + cAboveK[5] * phi[N_surface - 1][j]
                                + cAboveK[6] * phi[N_surface - 1][j + 1]
                                + cAboveK[1] * phi[N_surface][j + 1]
                ) / cAboveK[0];
                //phi[N_surface][j] = (1 - q) * prevPhi[N_surface][j] + q * phi[N_surface][j];
            }
            for (int i = N_surface - 1; i >= 1; i--) {
                newLocalCoord(i, 0);
                calculateS();
                calculateIntegral();
                calculateCUnderK();
                prevPhi[i][0] = phi[i][0];
                phi[i][0] = (
                        cUnderK[2] * phi[i + 1][0] + cUnderK[5] * phi[i - 1][0] + cUnderK[6] * phi[i - 1][1]
                                + cUnderK[1] * phi[i][1]
                ) / cUnderK[0];
                //phi[i][0] = (1 - q) * prevPhi[i][0] + q * phi[i][0];
            }
            int indexK = N_surface;
            newLocalCoord(indexK, 0);
            calculateS();
            calculateIntegral();
            calculateCforK();
            prevPhi[indexK][0] = phi[indexK][0];
            phi[indexK][0] = (
                    cForK[5] * phi[indexK - 1][0] + cForK[6] * phi[indexK - 1][1]
                            + cForK[1] * phi[indexK][1]
            ) / cForK[0];
            //phi[indexK][0] = (1 - q) * prevPhi[indexK][0] + q * phi[indexK][0];
            isCorrectResult = true;
            count++;
            for(int i = 1; i <= N_surface; i++){
                for(int j = 0; j <= N_sum - 1; j++){
                    if(Math.abs(phi[i][j] - prevPhi[i][j]) >= epsilon * q){
                        isCorrectResult = false;
                        //break;
                    }
                }
//                if(!isCorrectResult){
//                    break;
//                }
            }
        }
        System.out.println("count = " + count);
    }
}

