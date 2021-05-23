public class CreateGrid {
    private int N_inside;
    private int N_outside;
    private int N_surface;
    private int N   ;
    private int N_sum;
    private double[] r;
    private double[] z;
    private double K;

    private Pair[][] G;


    public CreateGrid(int N_inside, int N_outside, int N_surface){
        this.N_inside = N_inside;
        this.N_outside = N_outside;
        this.N_surface = N_surface;
        this.N_sum = N_inside + N_outside;
        this.N = 100;
        this.z = new double[N + 1];
        this.r = new double[N + 1];
        for(int i = 0; i <= N; i++){
            z[i] = Math.cos(Math.PI / 2.0 * i / N);
            r[i] = Math.sin(Math.PI / 2.0 * i / N);
        }
        this.K = N_surface * z[0] / (N_inside + N_surface);
        System.out.println("z0 = " + z[0]);
        System.out.println("K = " + K);
        initMatrG();
    }

    private double linInt(double A, double B, double t){
        return A + t * (B - A);
    }

    private Pair linInt(Pair A, Pair B, double t){
        return new Pair(
                A.getFirst() + t * (B.getFirst() - A.getFirst()),
                A.getSecond() + t * (B.getSecond() - A.getSecond())
        );
    }

    private Pair surf(double s){
        return new Pair(
                linInt(r[(int)Math.floor(s * N)], r[(int)Math.ceil(s * N)], s),
                linInt(z[(int)Math.floor(s * N)], z[(int)Math.ceil(s * N)], s)
        );
    }

    private void initMatrG(){
        G = new Pair[N_surface + 1][N_sum + 1];
        for(int i = 0; i <= N_surface; i++){
            G[i][0] = new Pair(0, 1.0 * i / N_surface * K);
            G[i][N_inside] = surf(1 - 1.0 * i / N_surface);
            G[i][N_sum] = new Pair(4 * surf(1 - 1.0 * i / N_surface).getFirst(), 4 * surf(1 - 1.0 * i / N_surface).getSecond());
        }
        for(int i = 1; i <= N_inside - 1; i++){
            G[0][i] = new Pair(1.0 * i / N_inside * r[N], 0);
            G[N_surface][i] = new Pair(0, linInt(K, z[0], 1.0 * i / N_inside));
        }
        for(int i = N_inside + 1; i <= N_sum; i++){
            G[0][i] = new Pair(linInt(r[N], 4 * r[N], (1.0 * i - N_inside)/ (N_sum - N_inside)), 0);
            G[N_surface][i] = new Pair(0, linInt(z[0], 4 * z[0], 1.0 * (i - N_inside) / (N_sum - N_inside)));
        }
        for(int i = 1; i <= N_surface - 1; i++){
            for(int j = 1; j <= N_sum; j++) {
                if (j <= N_inside - 1) {
                    G[i][j] = linInt(G[i][0], G[i][N_inside], 1.0 * j / N_inside);
                }
                else if(j >= N_inside + 1){
                    G[i][j] = linInt(G[i][N_inside], G[i][N_sum], 1.0 * (j - N_inside) / (N_sum - N_inside));
                }
            }
        }
    }

    public Pair[][] getG(){
        return G;
    }
    public int getIndexK(){
        for(int i = 0; i <= N_surface; i ++){
            if(G[i][0].getSecond() == K){
                return i;
            }
        }
        return -1;
    }
}
