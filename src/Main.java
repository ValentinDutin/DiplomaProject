public class Main {
    public static void main(String[] args) {
        int N_inside = 5;
        int N_outside = 5;
        int N_surface = 5;
        CreateGrid grid = new CreateGrid(N_inside, N_outside, N_surface);
        for(int i = 0; i <= N_surface; i++){
            for(int j = 0; j <= N_inside + N_outside; j++){
                System.out.format("(%f;%f) ",grid.getG()[i][j].getFirst(), grid.getG()[i][j].getSecond());
            }
            System.out.println();
        }
        System.out.println("\nPotential\n");
        MagneticFieldPotential mfp = new MagneticFieldPotential(grid.getG(), N_inside, N_outside, N_surface);
        mfp.calculatePotential();
        for(int i = 0; i <= N_surface; i++){
            for(int j = 0; j <= N_inside + N_outside; j++){
                System.out.format("(%f)   ", mfp.getPotential()[i][j]);
            }
            System.out.println();
        }
    }
}
