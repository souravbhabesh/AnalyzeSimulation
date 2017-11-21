//#define KAPPA 10.0
//#define STEPS 84000000
#define PERIOD 10000
//#define FRAMES STEPS/PERIOD
//#define NX 101
//#define NY 151
//#define LEN NX*NY
#define NMAX 50000
#define NXMAX 2000
#define MAXFRAMES 20000
//#define EPSILON 7200.0
#define a 1.0
//#define RUN 10

extern int N,Nb,Nd,bondGroup[NMAX*2],dihedralGroup[NMAX*4];
//N:#particles, Nb:#bonds, Nd:#dihedrals
extern float position[NMAX*3];
extern uint32_t particleID[NMAX];
extern char particleType[5][2];
extern double bendingEner[NMAX];
extern double bondHarmonicEner[NMAX];
extern double total_DHE,total_BHE;
extern double hgt_fluctuation[NMAX];
extern double h_avg_node[NMAX];
extern double hgt_fluctuation[NMAX];
extern double h_width[MAXFRAMES][NXMAX];
extern double h_bb[MAXFRAMES][NXMAX];


extern int nx,NY,RUN,STEPS,LEN,FRAMES;
extern double KAPPA,EPSILON;

