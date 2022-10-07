#define DEGREES 0.05
#define NODATA -99.0
#define PI 3.141592653589793
#define DEG2RAD (PI/180.0)
#define RAD2DEG 1/(PI/180.0)
#define MINDBZ 0.0

/* function prototypes. */

double degrees2distance(double gridLatitude, double latitude, 
			double gridLongitude, double longitude);
double cressman(double d2gp, float RADIUS);
void searcher(double distance[], int indices[], int indices2[], 
	      int *count, int *count2, float RADIUS);
void searcher2(double *latitude_2d, double *longitude_2d, 
	       double *gridLatitude, double *gridLongitude,  
	       double *reflectivity_2d, int indices2[], 
	       int *count2, int size2D, int i);
void interp(int size2D, int gridColumns, int gridRows, int numAlts,
	    float RADIUS, char *interp_type, double *latitude, double *longitude,
	    double *inputArray, double *gridLatitude,double *gridLongitude, 
        double *gridNoSwath, double *outputArray);

